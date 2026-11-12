# Isolated packing experiment; never edits runtime sources or source records.
script <- readLines('tools/long-term-audit/check-flagship-label-quality.R', warn = FALSE)
cut <- grep('^results <- list', script)[1L]
eval(parse(text = script[seq_len(cut)]))
namespace <- asNamespace('ggchord')
candidate <- get('ggchord_place_dense_perimeter', namespace)
replacements <- 0L
events <- new.env(parent = emptyenv()); events$rows <- list()
rewrite <- function(expr) {
  if (!is.call(expr)) return(expr)
  if (identical(expr[[1L]], as.name('if')) &&
      grepl('^ggchord_perimeter_rail_crossings', paste(deparse(expr[[2L]]), collapse = ' '))) {
    replacements <<- replacements+1L
    return(quote({
      # Search a bounded lateral/vertical set jointly. A wider departure
      # boundary can avoid crossings without lifting the whole column.
      original_axis <- axis; original_boundary <- boundary
      step <- max(.008, .04*units_per_inch)
      shifts <- unique(c(0, -mean(range(axis)),
        as.vector(rbind(-seq_len(60L)*step, seq_len(60L)*step))))
      best <- NULL
      for (factor in c(1, 1.2, 1.5, 2, 2.5, 3)) {
        edge <- original_boundary*factor
        width <- 2*(abs(edge)+max(tasks$w[rows]))
        for (shift in shifts) {
          yy <- original_axis+shift
          span <- range(c(-backbone_outer, backbone_outer,
            yy-tasks$h[rows]/2, yy+tasks$h[rows]/2))
          score <- max(width/max(canvas[1]/canvas[2], 1e-8), diff(span)) +
            .01*abs(mean(span))
          if (!is.null(best) && score >= best$score) next
          if (ggchord_perimeter_rail_crossings(tasks, rows, edge, yy) != 0L) next
          best <- list(score = score, boundary = edge, axis = yy,
            factor = factor, shift = shift)
        }
      }
      if (!is.null(best)) {boundary <- best$boundary; axis <- best$axis}
      events$rows <- c(events$rows, list(list(side = side,
        count = length(rows), accepted = !is.null(best),
        factor = if (!is.null(best)) best$factor else NA_real_,
        shift = if (!is.null(best)) best$shift else NA_real_)))
    }))
  }
  as.call(lapply(as.list(expr), rewrite))
}
body(candidate) <- rewrite(body(candidate)); stopifnot(replacements == 1L)
private <- new.env(parent = namespace); private$events <- events
environment(candidate) <- private
unlockBinding('ggchord_place_dense_perimeter', namespace)
assign('ggchord_place_dense_perimeter', candidate, namespace)
lockBinding('ggchord_place_dense_perimeter', namespace)
prefix <- 'v0.27.0-balanced-perimeter-prototype-'
eval(parse(text = script[-seq_len(cut)]))
file <- paste0('tools/long-term-audit/', prefix, paste(names(profiles), collapse = '-'), '-results.json')
report <- jsonlite::fromJSON(file, simplifyVector = FALSE)
report$baseline_runtime_sha256 <- report$runtime_sha256; report$runtime_sha256 <- NULL
report$production_override <- 'Joint bounded lateral/vertical rail candidate search; no production certification'
report$candidate_body_sha256 <- digest::digest(deparse(body(candidate)), algo = 'sha256')
report$candidate_events <- events$rows
jsonlite::write_json(report, file, pretty = TRUE, auto_unbox = TRUE, na = 'null')
