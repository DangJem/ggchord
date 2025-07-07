#' Generate curved sequence paths
#'
#' Generates smooth sequence paths (supporting straight lines, arcs, and custom curvatures) based on start angle, end angle, radius, and curvature.
#'
#' @param start_angle Numeric, start angle (in radians)
#' @param end_angle Numeric, end angle (in radians)
#' @param radius Numeric, path radius
#' @param curvature Numeric, curvature (0 = straight line, 1 = standard arc, >1 = more curved)
#' @param n_points Integer, number of points in the path (controls smoothness), default 100
#' @return data.frame containing columns x, y (coordinates of points on the path)
#' @keywords internal
#'
#' @importFrom stats median setNames
#' @importFrom utils tail
#' @importFrom graphics text
#'
generate_curvature_path <- function(start_angle, end_angle, radius, curvature, n_points = 100) {
  # Restrict curvature range
  curvature <- max(-0.99, min(10, curvature))

  # Calculate center angle and total angle
  center_angle <- (start_angle + end_angle) / 2
  total_angle <- end_angle - start_angle

  # Generate uniformly distributed angles
  angles <- seq(start_angle, end_angle, length.out = n_points)

  # Ensure fixed positions for start and end points
  fixed_start_x <- radius * cos(start_angle)
  fixed_start_y <- radius * sin(start_angle)
  fixed_end_x <- radius * cos(end_angle)
  fixed_end_y <- radius * sin(end_angle)

  # Adjust path based on curvature
  if (curvature == 0) {
    # Straight line
    t <- seq(0, 1, length.out = n_points)
    x <- fixed_start_x + t * (fixed_end_x - fixed_start_x)
    y <- fixed_start_y + t * (fixed_end_y - fixed_start_y)
  } else if (curvature == 1) {
    # Original circular path
    x <- radius * cos(angles)
    y <- radius * sin(angles)
  } else {
    # Use Bézier curve for smooth curvature, ensuring fixed start and end points
    # Control point position adjusted based on curvature
    mid_angle <- center_angle
    mid_radius_factor <- 1 + (curvature - 1) * 0.3  # Reduce curvature's impact on control points
    mid_radius <- radius * mid_radius_factor

    # Control point coordinates
    control_x <- mid_radius * cos(mid_angle)
    control_y <- mid_radius * sin(mid_angle)

    # Calculate path points using Bézier curve
    x <- numeric(n_points)
    y <- numeric(n_points)

    for (i in 1:n_points) {
      t <- (i - 1) / (n_points - 1)

      # Quadratic Bézier curve formula
      x[i] <- (1 - t)^2 * fixed_start_x + 2 * (1 - t) * t * control_x + t^2 * fixed_end_x
      y[i] <- (1 - t)^2 * fixed_start_y + 2 * (1 - t) * t * control_y + t^2 * fixed_end_y
    }
  }

  # Ensure no NaN or Inf values
  x[is.na(x) | is.infinite(x)] <- 0
  y[is.na(y) | is.infinite(y)] <- 0

  data.frame(x = x, y = y)
}


#' Generate Bézier curve points
#'
#' Generates points on a Bézier curve based on start point, end point, and control points (for smooth ribbons)
#'
#' @param p0 Numeric vector (length 2), start point coordinates (x, y)
#' @param p3 Numeric vector (length 2), end point coordinates (x, y)
#' @param c1 Numeric vector (length 2), first control point coordinates (x, y)
#' @param c2 Numeric vector (length 2), second control point coordinates (x, y)
#' @param n Integer, number of curve points (controls smoothness), default 100
#' @return data.frame containing columns x, y (coordinates of points on the Bézier curve)
#' @keywords internal
bezier_pts <- function(p0, p3, c1, c2, n = 100) {
  t <- seq(0, 1, length.out = n)
  bx <- (1 - t)^3*p0[1] + 3*(1 - t)^2*t*c1[1] + 3*(1 - t)*t^2*c2[1] + t^3*p3[1]
  by <- (1 - t)^3*p0[2] + 3*(1 - t)^2*t*c1[2] + 3*(1 - t)*t^2*c2[2] + t^3*p3[2]
  data.frame(x = bx, y = by)
}


#' Generate major axis tick breakpoints
#'
#' Generates uniform and visually appealing major tick positions based on sequence length and target tick count (avoids excessively short end ticks)
#'
#' @param max_value Numeric, sequence length (maximum value)
#' @param n Integer, target number of ticks, default 5
#' @param tol Numeric (0-1), tolerance threshold for end tick length (proportion of the median length of other ticks), default 0.5
#' @return Numeric vector, major tick positions (including 0 and max_value)
#' @keywords internal
breakPointsFunc <- function(max_value, n = 5, tol = 0.5) {
  if (max_value <= 0) return(c(0, max_value))

  # 1. Generate approximately n ticks using pretty()
  ticks <- pretty(c(0, max_value), n = n)
  ticks <- ticks[ticks >= 0 & ticks <= max_value]
  ticks <- sort(unique(c(0, ticks, max_value)))

  # 2. If the last segment is too small (less than median of other segments * tol), remove the penultimate tick
  if (length(ticks) >= 3) {
    d <- diff(ticks)
    # Median of other segments (excluding the last one)
    med <- median(d[-length(d)])
    # Remove penultimate point if the last segment is too small
    if (d[length(d)] < med * tol) {
      ticks <- ticks[-(length(ticks) - 1)]
    }
  }

  return(ticks)
}


#' Custom gene arrow legend drawing function
#'
#' Generates gene arrow-shaped legend symbols (polygons) for ggplot2 legends
#'
#' @param data Legend data (contains aesthetic mapping parameters like fill, colour, size)
#' @param params Legend parameters (automatically passed by ggplot2)
#' @param size Legend symbol size
#' @return grid::polygonGrob object, gene arrow-shaped legend symbol
#' @keywords internal
draw_key_gene_arrow <- function(data, params, size) {
  # Five vertices: rectangle bottom-left, rectangle top-left, rectangle top-right, arrow tip, rectangle bottom-right
  x_pts <- unit(c(0.1, 0.1, 0.6, 0.9, 0.6), "npc")
  y_pts <- unit(c(0.2, 0.8, 0.8, 0.5, 0.2), "npc")

  polygonGrob(
    x = x_pts, y = y_pts,
    gp = gpar(
      fill = alpha(if_null_else(data$fill, "grey"), if_null_else(data$alpha, 1)),
      col = if_null_else(data$colour, "black"),
      lwd = if_null_else(data$size, 0.5) * .pt
    )
  )
}
