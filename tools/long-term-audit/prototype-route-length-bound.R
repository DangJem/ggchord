# Prune only graph vertices that cannot occur on an already permitted route.
# The original annotation caller accepts at most 3*old_length+5 physical mm.
script <- readLines('tools/long-term-audit/check-label-count-transition.R', warn = FALSE)
cut <- grep('^results <- list', script)[1L]
stopifnot(is.finite(cut))
eval(parse(text = script[seq_len(cut)]))
namespace <- asNamespace('ggchord')
bounded <- get('ggchord_route_shortest', envir = namespace)
formals(bounded) <- c(formals(bounded), alist(max_distance = Inf))
parts <- as.list(body(bounded))
at <- which(vapply(parts, function(x) grepl('nodes <- unique(rbind(start, end, corners))',
  paste(deparse(x), collapse = ' '), fixed = TRUE), logical(1)))
stopifnot(length(at) == 1L)
trace <- new.env(parent = emptyenv())
trace$calls <- 0L; trace$bounded_calls <- 0L
trace$maximum_nodes_before <- 0L; trace$maximum_nodes_after <- 0L
prune <- quote(if (is.finite(max_distance)) {
  route_bound_trace$maximum_nodes_before <- max(route_bound_trace$maximum_nodes_before, nrow(nodes))
  from_start <- sqrt(rowSums(sweep(nodes, 2, start, '-')^2))
  from_end <- sqrt(rowSums(sweep(nodes, 2, end, '-')^2))
  nodes <- nodes[from_start+from_end <= max_distance+1e-9, , drop = FALSE]
  route_bound_trace$maximum_nodes_after <- max(route_bound_trace$maximum_nodes_after, nrow(nodes))
})
body(bounded) <- as.call(append(parts, list(prune), after = at))
private <- new.env(parent = environment(bounded)); private$route_bound_trace <- trace
environment(bounded) <- private
wrapper <- function(start, end, boxes, canvas, clearance = .6, barriers = NULL) {
  trace$calls <- trace$calls+1L
  parent <- parent.frame()
  path <- get0('path', envir = parent, inherits = FALSE)
  unit <- get0('units', envir = parent, inherits = FALSE)
  maximum <- Inf
  if (is.data.frame(path) && nrow(path) &&
      all(c('x1', 'y1', 'x2', 'y2') %in% names(path)) && length(unit) == 1L) {
    maximum <- 3*sum(sqrt((path$x2-path$x1)^2+(path$y2-path$y1)^2))+5*unit/25.4
    trace$bounded_calls <- trace$bounded_calls+1L
  }
  bounded(start, end, boxes, canvas, clearance, barriers, max_distance = maximum)
}
unlockBinding('ggchord_route_shortest', namespace)
assign('ggchord_route_shortest', wrapper, envir = namespace)
lockBinding('ggchord_route_shortest', namespace)
prefix <- 'v0.27.0-label-transition-bound-prototype-'
remaining <- script[-seq_len(cut)]
remaining <- gsub('v0.27.0-label-transition-', 'v0.27.0-label-transition-bound-prototype-', remaining, fixed = TRUE)
eval(parse(text = remaining))
file <- paste0('tools/long-term-audit/', prefix, paste(counts, collapse = '-'), '-results.json')
report <- jsonlite::fromJSON(file, simplifyVector = FALSE)
report$baseline_runtime_sha256 <- report$runtime_sha256; report$runtime_sha256 <- NULL
report$production_override <- 'Graph nodes violating the existing accepted-length lower bound are pruned only for bounded annotation calls'
report$candidate_body_sha256 <- digest::digest(deparse(body(bounded)), algo = 'sha256')
report$search_observations <- as.list(trace)
report$scope <- paste('Isolated route-length candidate for', length(counts), 'selected count(s); production source unchanged.')
jsonlite::write_json(report, file, pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Route-length-bound candidate observations recorded\n')
