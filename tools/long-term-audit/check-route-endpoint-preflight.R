# Isolated solver preflight equivalence on explicit geometric controls.
pkgload::load_all(quiet = TRUE)
production <- ggchord_route_shortest
original <- production
old_parts <- as.list(body(original))
old_parts <- Filter(function(x) !grepl('ggchord_route_segment_hits(start, start, padded)',
  paste(deparse(x), collapse = ' '), fixed = TRUE), old_parts)
body(original) <- as.call(old_parts)
candidate <- original
body(candidate) <- as.call(c(list(as.name('{')), as.list(body(original))[-1L]))
parts <- as.list(body(candidate))
at <- which(vapply(parts, function(x) grepl('padded <- ggchord_route_inflate', paste(deparse(x), collapse = ' '), fixed = TRUE), logical(1)))
stopifnot(length(at) == 1L)
preflight <- quote(if (any(ggchord_route_segment_hits(start, start, padded)) ||
    any(ggchord_route_segment_hits(end, end, padded))) return(NULL))
body(candidate) <- as.call(append(parts, list(preflight), after = at))
box <- data.frame(xmin = 4, xmax = 6, ymin = 4, ymax = 6)
barrier <- data.frame(x1 = 5, y1 = 2, x2 = 5, y2 = 8)
cases <- list(start_inside = list(start = c(5, 5), end = c(9, 9)),
  end_inside = list(start = c(1, 1), end = c(5, 5)),
  padded_endpoint = list(start = c(3.8, 5), end = c(9, 9)),
  boundary = list(start = c(3.75, 5), end = c(9, 9)),
  route_around = list(start = c(1, 5), end = c(9, 5)),
  clear = list(start = c(1, 1), end = c(9, 1)))
results <- lapply(cases, function(case) {
  args <- c(case, list(boxes = box, canvas = c(10, 10), clearance = .25, barriers = barrier))
  before <- do.call(original, args); after <- do.call(candidate, args)
  stopifnot(identical(before, after))
  list(same_result = TRUE, status = if (is.null(after)) 'no_route' else 'route', path = after)
})
same_production <- identical(parse(text = deparse(body(production)), keep.source = FALSE),
  parse(text = deparse(body(candidate)), keep.source = FALSE))
stopifnot(same_production)
jsonlite::write_json(list(status = 'passed production preflight equivalence', cases = results,
  production_has_preflight = same_production,
  source_sha256 = digest::digest(file = 'R/fusion-routing.R', algo = 'sha256'),
  scope = 'Six explicit geometry controls; dense figure observations are recorded separately'),
  'tools/long-term-audit/v0.27.0-route-endpoint-prototype-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Endpoint preflight production equivalence: passed\n')
