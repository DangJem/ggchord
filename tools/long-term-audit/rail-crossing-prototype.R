# Preserve the existing proper-interior crossing rule, batching later rails.
pkgload::load_all(quiet = TRUE)
candidate_rail_crossings <- function(tasks, rows, edge_x, edge_y) {
  n <- length(rows)
  if (n < 2L) return(0L)
  segments <- data.frame(x1 = tasks$route_x[rows], y1 = tasks$route_y[rows],
    x2 = rep(edge_x, n), y2 = edge_y)
  crossings <- 0L
  for (a in seq_len(n-1L)) {
    at <- seq.int(a+1L, n)
    at <- at[abs(segments$x1[at]-segments$x1[a]) >= 1e-8 |
      abs(segments$y1[at]-segments$y1[a]) >= 1e-8]
    if (length(at)) crossings <- crossings+sum(ggchord_route_segment_crossings(
      c(segments$x1[a], segments$y1[a]), c(segments$x2[a], segments$y2[a]),
      segments[at, , drop = FALSE]))
  }
  as.integer(crossings)
}
set.seed(270)
for (n in c(0L, 1L, 2L, 8L, 31L, 58L)) {
  tasks <- data.frame(route_x = runif(n, -1, 1), route_y = runif(n, -1, 1))
  if (n > 2L) {
    tasks[2, ] <- tasks[1, ]
    tasks[3, ] <- tasks[1, ]+1e-10
  }
  for (edge_x in c(-2, 0, 2)) {
    edge_y <- runif(n, -2, 2)
    expected <- ggchord_perimeter_rail_crossings(tasks, seq_len(n), edge_x, edge_y)
    observed <- candidate_rail_crossings(tasks, seq_len(n), edge_x, edge_y)
    stopifnot(identical(observed, expected))
  }
}
cat('Batched rail crossing counts: 18 explicit geometry controls equivalent\n')
