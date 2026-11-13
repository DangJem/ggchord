# Curve coordinate mapping function
map_to_curve_many <- function(angle, radius, ref) {
  n <- length(angle)
  nref <- length(ref$angles)
  periodic <- isTRUE(ref$periodic)
  if (periodic) {
    period <- diff(range(ref$angles))
    angle <- min(ref$angles) + (angle - min(ref$angles)) %% period
  }

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

  if (periodic) {
    # The final reference duplicates the first. Use neighbours on both sides
    # of the seam rather than clamping to a one-sided endpoint tangent.
    canonical <- ((idx - 1L) %% (nref - 1L)) + 1L
    previous <- ((canonical - 2L) %% (nref - 1L)) + 1L
    following <- (canonical %% (nref - 1L)) + 1L
    dx <- ref$path$x[following] - ref$path$x[previous]
    dy <- ref$path$y[following] - ref$path$y[previous]
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
  if (isTRUE(ref$periodic)) return(as.numeric(map_to_curve_many(angle, radius, ref)[1L, ]))
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


# Only genuinely closed drawing paths may use periodic feature geometry.
# A circular view with a visible gap still has two separate display ends.
ggchord_closed_reference <- function(ref) {
  path <- ref$path
  if (nrow(path) < 4L) return(FALSE)
  extent <- max(diff(range(path$x)), diff(range(path$y)), 1e-8)
  sqrt((path$x[1L]-tail(path$x,1L))^2 +
    (path$y[1L]-tail(path$y,1L))^2) <= extent*1e-8 &&
    abs(diff(range(ref$angles))-2*pi) < 1e-6
}
