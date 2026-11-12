pkgload::load_all(quiet = TRUE)
# Explicit computational geometry, no biological observations.
blocked <- function(a, targets, boxes, barriers, tolerance = 1e-7) {
  n <- nrow(targets)
  if (!n) return(logical())
  result <- rep(FALSE, n)
  if (nrow(boxes)) {
    enter <- matrix(0, n, nrow(boxes)); leave <- matrix(1, n, nrow(boxes))
    for (axis in 1:2) {
      lo <- boxes[[if (axis == 1L) 'xmin' else 'ymin']] + tolerance
      hi <- boxes[[if (axis == 1L) 'xmax' else 'ymax']] - tolerance
      delta <- targets[, axis]-a[axis]
      moving <- abs(delta) >= 1e-12
      if (any(moving)) {
        first <- outer(rep(1, sum(moving)), lo-a[axis])/delta[moving]
        last <- outer(rep(1, sum(moving)), hi-a[axis])/delta[moving]
        enter[moving, ] <- pmax(enter[moving, , drop = FALSE], pmin(first, last))
        leave[moving, ] <- pmin(leave[moving, , drop = FALSE], pmax(first, last))
      }
      outside <- a[axis] <= lo | a[axis] >= hi
      if (any(!moving) && any(outside)) leave[!moving, outside] <- -Inf
    }
    result <- rowSums(enter < leave & leave > 0 & enter < 1) > 0
  }
  if (!is.null(barriers) && nrow(barriers)) {
    rx <- targets[, 1]-a[1]; ry <- targets[, 2]-a[2]
    sx <- barriers$x2-barriers$x1; sy <- barriers$y2-barriers$y1
    denominator <- outer(rx, sy)-outer(ry, sx)
    valid <- is.finite(denominator) & abs(denominator) > 1e-9
    qx <- barriers$x1-a[1]; qy <- barriers$y1-a[2]
    t <- matrix(rep(qx*sy-qy*sx, each = n), n)/denominator
    u <- (outer(ry, qx)-outer(rx, qy))/denominator
    result <- result | rowSums(valid & t > 1e-9 & t < 1-1e-9 & u > 1e-9 & u < 1-1e-9,
      na.rm = TRUE) > 0
  }
  unname(result)
}
set.seed(27)
boxes <- data.frame(xmin = runif(35, 1, 8), ymin = runif(35, 1, 8))
boxes$xmax <- boxes$xmin+runif(35, .05, 1); boxes$ymax <- boxes$ymin+runif(35, .05, 1)
barriers <- as.data.frame(matrix(runif(4*40, 0, 10), ncol = 4))
names(barriers) <- c('x1', 'y1', 'x2', 'y2')
targets <- rbind(matrix(runif(2*100, 0, 10), ncol = 2), c(0, 0), c(10, 10), c(5, 0), c(0, 5))
targets <- rbind(targets, as.matrix(boxes[1:3, c('xmin', 'ymin')]),
  as.matrix(boxes[1:3, c('xmax', 'ymax')]),
  as.matrix(boxes[1:3, c('xmin', 'ymin')])+1e-7)
origins <- rbind(matrix(runif(2*15, 0, 10), ncol = 2), c(0, 0), c(5, 5))
checks <- 0L
for (i in seq_len(nrow(origins))) for (use_boxes in c(FALSE, TRUE)) for (use_barriers in c(FALSE, TRUE)) {
  bs <- boxes[if (use_boxes) seq_len(nrow(boxes)) else integer(), ]
  br <- barriers[if (use_barriers) seq_len(nrow(barriers)) else integer(), ]
  expected <- vapply(seq_len(nrow(targets)), function(j)
    any(ggchord_route_segment_hits(origins[i, ], targets[j, ], bs)) ||
      any(ggchord_route_segment_crossings(origins[i, ], targets[j, ], br)), logical(1))
  actual <- ggchord_route_targets_blocked(origins[i, ], targets, bs, br)
  stopifnot(identical(actual, expected)); checks <- checks+length(actual)
}
cat('Vector obstruction decisions equivalent:', checks, '\n')
source('tools/long-term-audit/route-scalar-reference.R')
original <- reference_route_shortest
candidate <- ggchord_route_shortest
cases <- list(
  list(start = c(0, 0), end = c(10, 10), boxes = boxes[1:3, ], barriers = barriers[1:3, ]),
  list(start = c(0, 5), end = c(10, 5), boxes = boxes[1:8, ], barriers = barriers[1:8, ]),
  list(start = c(5, 0), end = c(5, 10), boxes = boxes[1:10, ], barriers = barriers[1:12, ]),
  list(start = c(1, 1), end = c(9, 9), boxes = boxes[1:2, ], barriers = barriers[1:2, ]))
cases$fixed_inside <- list(start = c(mean(c(boxes$xmin[1], boxes$xmax[1])),
  mean(c(boxes$ymin[1], boxes$ymax[1]))), end = c(10, 10), boxes = boxes[1, ], barriers = barriers[0, ])
route_status <- character()
for (case in cases) {
  args <- c(case, list(canvas = c(10, 10), clearance = .05))
  a <- do.call(original, args); b <- do.call(candidate, args)
  stopifnot(identical(a, b))
  route_status <- c(route_status, if (is.null(b)) 'no_route' else 'route')
}
stopifnot(any(route_status == 'route'), any(route_status == 'no_route'))
cat('Complete graph routes equivalent, including an explicit no-route case\n')

for (n in c(0L, 1L, 2L, 8L, 31L, 58L)) {
  tasks <- data.frame(route_x = runif(n, -1, 1), route_y = runif(n, -1, 1))
  if (n > 2L) { tasks[2, ] <- tasks[1, ]; tasks[3, ] <- tasks[1, ]+1e-10 }
  for (edge_x in c(-2, 0, 2)) {
    edge_y <- runif(n, -2, 2)
    stopifnot(identical(reference_rail_crossings(tasks, seq_len(n), edge_x, edge_y),
      ggchord_perimeter_rail_crossings(tasks, seq_len(n), edge_x, edge_y)))
  }
}
old <- data.frame(.component = 'segment', source_row = 1,
  x = c(0, .2), y = c(0, .1), xend = c(.2, 1), yend = c(.1, 0))
new <- data.frame(.component = 'segment', source_row = 1,
  x = c(0, .1, .5), y = c(0, -.2, -.1), xend = c(.1, .5, 1), yend = c(-.2, -.1, 0))
a <- ggchord_annotation_leader_roots(list(format = list(gene_label_repel = old)))
b <- ggchord_annotation_leader_roots(list(format = list(gene_label_repel = new)))
stopifnot(identical(a, b), nrow(a) == 1L, a$x0 == 0, a$y0 == 0)
cat('Native rail counts and fixed directed roots: passed\n')
jsonlite::write_json(list(status = 'passed native batching equivalence controls',
  decision_count = checks, full_route_status = route_status,
  candidate_body_sha256 = digest::digest(deparse(body(candidate)), algo = 'sha256'),
  helper_body_sha256 = digest::digest(deparse(body(ggchord_route_targets_blocked)), algo = 'sha256'),
  reference_commit = '665651c', runtime_source_sha256 = digest::digest(file = 'R/fusion-routing.R', algo = 'sha256'),
  production_changed = TRUE, scope = 'Explicit finite computational geometry; no full-figure or platform certification'),
  'tools/long-term-audit/v0.27.0-route-batching-results.json', pretty = TRUE, auto_unbox = TRUE)
