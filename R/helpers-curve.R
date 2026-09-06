# Curve coordinate mapping function
map_to_curve_many <- function(angle, radius, ref) {
  n <- length(angle)
  nref <- length(ref$angles)

  # Vectorised nearest-reference lookup (O(n log n) overall).
  fi <- findInterval(angle, ref$angles)
  fi <- pmax(1, pmin(fi, nref - 1))
  idx <- ifelse(abs(ref$angles[fi] - angle) <=
                  abs(ref$angles[fi + 1] - angle), fi, fi + 1)

  base <- ref$path[idx, , drop = FALSE]
  idx_next <- pmin(idx + 1, nref)
  idx_prev <- pmax(idx - 1, 1)
  dx <- ref$path$x[idx_next] - base$x
  dy <- ref$path$y[idx_next] - base$y
  last <- idx == nref
  if (any(last)) {
    dx[last] <- base$x[last] - ref$path$x[idx_prev[last]]
    dy[last] <- base$y[last] - ref$path$y[idx_prev[last]]
  }

  norm_x <- -dy
  norm_y <- dx
  nl <- sqrt(norm_x^2 + norm_y^2)
  ok <- nl > 0
  norm_x[ok] <- norm_x[ok] / nl[ok]
  norm_y[ok] <- norm_y[ok] / nl[ok]

  offset <- radius - ref$r0
  cbind(x = base$x + norm_x * offset,
        y = base$y + norm_y * offset)
}

map_to_curve <- function(angle, radius, ref) {
  # ref$angles is sorted, so use findInterval (O(log n)) instead of a
  # full linear scan to locate the nearest reference angle.
  fi <- findInterval(angle, ref$angles)
  if (fi < 1) fi <- 1
  if (fi >= length(ref$angles)) fi <- length(ref$angles) - 1
  idx <- if (abs(ref$angles[fi] - angle) <= abs(ref$angles[fi + 1] - angle)) fi else fi + 1
  base <- ref$path[idx, ]
  if (idx < nrow(ref$path)) {
    dx <- ref$path$x[idx + 1] - base$x
    dy <- ref$path$y[idx + 1] - base$y
  } else {
    dx <- base$x - ref$path$x[idx - 1]
    dy <- base$y - ref$path$y[idx - 1]
  }
  norm <- c(-dy, dx)
  nl <- sqrt(sum(norm^2))
  if (nl > 0) norm <- norm / nl
  offset <- radius - ref$r0
  c(x = base$x + norm[1] * offset,
    y = base$y + norm[2] * offset)
}

