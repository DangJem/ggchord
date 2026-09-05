#' Generate curved sequence paths
#'
#' Generates smooth sequence paths (supporting straight lines, arcs, and custom curvatures) based on start angle, end angle, radius, and curvature.
#'
#' @param start_angle Numeric, start angle (in radians)
#' @param end_angle Numeric, end angle (in radians)
#' @param radius Numeric, path radius
#' @param curvature Numeric, curvature (0 = straight chord, 1 = standard arc, negative = opposite bow)
#' @param n_points Integer, number of points in the path (controls smoothness), default 100
#' @return data.frame containing columns x, y (coordinates of points on the path)
#' @keywords internal
#'
#' @importFrom stats median setNames
#' @importFrom utils tail
#' @importFrom graphics text
#'
generate_curvature_path <- function(start_angle, end_angle, radius, curvature, n_points = 100) {
  if (length(curvature) != 1L || !is.finite(curvature)) {
    ggchord_stop("curvature must be one finite number")
  }
  angles <- seq(start_angle, end_angle, length.out = n_points)
  center <- (start_angle + end_angle) / 2
  half_span <- (end_angle - start_angle) / 2
  # Scale the signed perpendicular distance from the endpoint chord. This
  # is continuous through 0 and 1; -c reflects +c across that chord.
  # No asymmetric clipping: large positive and negative values stay distinct.
  t <- seq(0, 1, length.out = n_points)
  chord_along <- radius * (2 * t - 1) * sin(half_span)
  circle_along <- radius * sin(angles - center)
  # Fade tangential reparameterisation continuously to the straight chord.
  # This also keeps curvature=0 inside its endpoints for arcs over 180 degrees.
  along <- chord_along + min(abs(curvature), 1) * (circle_along - chord_along)
  across <- radius * (cos(half_span) +
    curvature * (cos(angles - center) - cos(half_span)))
  x <- across * cos(center) - along * sin(center)
  y <- across * sin(center) + along * cos(center)
  if (any(!is.finite(c(x, y)))) ggchord_stop("Curvature produces non-finite coordinates")
  data.frame(x = x, y = y)
}



#' Generate Bezier curve points
#'
#' Generates points on a Bezier curve based on start point, end point, and control points (for smooth ribbons)
#'
#' @param p0 Numeric vector (length 2), start point coordinates (x, y)
#' @param p3 Numeric vector (length 2), end point coordinates (x, y)
#' @param c1 Numeric vector (length 2), first control point coordinates (x, y)
#' @param c2 Numeric vector (length 2), second control point coordinates (x, y)
#' @param n Integer, number of curve points (controls smoothness), default 100
#' @return data.frame containing columns x, y (coordinates of points on the Bezier curve)
#' @keywords internal
bezier_pts <- function(p0, p3, c1, c2, n = 100) {
  t <- seq(0, 1, length.out = n)
  bx <- (1 - t)^3*p0[1] + 3*(1 - t)^2*t*c1[1] + 3*(1 - t)*t^2*c2[1] + t^3*p3[1]
  by <- (1 - t)^3*p0[2] + 3*(1 - t)^2*t*c1[2] + 3*(1 - t)*t^2*c2[2] + t^3*p3[2]
  cbind(x = bx, y = by)
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
  # Match geom_gene(): a constant-width body followed by a linearly tapered
  # head, without the wider shoulder used by a conventional block arrow.
  # Strand overrides mirror the same polygon horizontally.
  x <- c(0.10, 0.62, 0.90, 0.62, 0.10)
  if (identical(as.character(data$strand %||% "+")[1], "-")) x <- 1 - x
  x_pts <- grid::unit(x, "npc")
  y_pts <- grid::unit(c(0.32, 0.32, 0.50, 0.68, 0.68), "npc")

  grid::polygonGrob(
    x = x_pts, y = y_pts,
    gp = grid::gpar(
      fill = ggplot2::alpha(if_null_else(data$fill, "grey"), if_null_else(data$alpha, 1)),
      col = if_null_else(data$colour, "#353A3E"),
      lwd = if_null_else(data$linewidth %||% data$size, 0.5) * ggplot2::.pt
    )
  )
}

#' Key glyph for sequence legends
#'
#' Draws the path symbol only when the key data contains colour; otherwise returns a blank (prevents ggplot2 4.x
#' from mixing unrelated layers into other legends with default grey/black symbols).
#' @keywords internal
key_glyph_seq <- function(data, params, size) {
  data$colour <- data$seq_colour %||% data$colour
  if (is.null(data$colour)) return(ggplot2::zeroGrob())
  col <- ggplot2::alpha(data$colour, data$alpha %||% 1)
  lwd <- (data$linewidth %||% data$size %||% 0.8) * ggplot2::.pt
  grid::grobTree(
    grid::segmentsGrob(
      x0 = grid::unit(0.12, "npc"), x1 = grid::unit(0.74, "npc"),
      y0 = grid::unit(0.5, "npc"), y1 = grid::unit(0.5, "npc"),
      gp = grid::gpar(col = col, lwd = lwd, lineend = "round")
    ),
    grid::polygonGrob(
      x = grid::unit(c(0.70, 0.90, 0.70), "npc"),
      y = grid::unit(c(0.34, 0.5, 0.66), "npc"),
      gp = grid::gpar(fill = col, col = col)
    )
  )
}

#' Key glyph for ribbon legends
#'
#' Draws the polygon symbol only when the key data contains fill; otherwise returns a blank.
#' @keywords internal
key_glyph_ribbon <- function(data, params, size) {
  data$fill <- data$ribbon_fill %||% data$fill
  data$colour <- data$ribbon_colour %||% data$colour
  data$alpha <- data$ribbon_alpha %||% data$alpha
  if (is.null(data$fill)) return(ggplot2::zeroGrob())
  ggplot2::draw_key_polygon(data, params, size)
}

#' Key glyph for gene arrow legends
#'
#' Draws the gene arrow only when the key data contains fill; otherwise returns a blank.
#' @keywords internal
key_glyph_gene <- function(data, params, size) {
  data$fill <- data$gene_fill %||% data$feature_fill %||% data$fill
  if (is.null(data$fill)) return(ggplot2::zeroGrob())
  draw_key_gene_arrow(data, params, size)
}

#' Key glyph for generic genomic features
#' @keywords internal
key_glyph_feature <- function(data, params, size) {
  data$fill <- data$feature_fill %||% data$fill
  if (is.null(data$fill)) return(ggplot2::zeroGrob())

  shape <- as.character(data$feature_shape %||% "arrow")[1]
  if (!shape %in% c("arrow", "block", "chevron", "lollipop")) {
    shape <- "arrow"
  }
  col <- ggplot2::alpha(data$colour %||% "#353A3E", data$alpha %||% 1)
  fill <- ggplot2::alpha(data$fill, data$alpha %||% 1)
  lwd <- (data$linewidth %||% data$size %||% 0.25) * ggplot2::.pt
  gp <- grid::gpar(col = col, fill = fill, lwd = lwd, linejoin = "round")

  if (identical(shape, "arrow")) {
    return(draw_key_gene_arrow(data, params, size))
  }
  if (identical(shape, "block")) {
    return(grid::rectGrob(
      x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
      width = grid::unit(0.76, "npc"), height = grid::unit(0.34, "npc"), gp = gp
    ))
  }
  if (identical(shape, "chevron")) {
    return(grid::polygonGrob(
      x = grid::unit(c(0.10, 0.62, 0.90, 0.62, 0.10, 0.34), "npc"),
      y = grid::unit(c(0.32, 0.32, 0.50, 0.68, 0.68, 0.50), "npc"),
      gp = gp
    ))
  }

  grid::grobTree(
    grid::segmentsGrob(
      x0 = grid::unit(0.18, "npc"), x1 = grid::unit(0.62, "npc"),
      y0 = grid::unit(0.5, "npc"), y1 = grid::unit(0.5, "npc"),
      gp = grid::gpar(col = col, lwd = max(lwd, 0.7), lineend = "round")
    ),
    grid::circleGrob(
      x = grid::unit(0.70, "npc"), y = grid::unit(0.5, "npc"),
      r = grid::unit(0.18, "npc"), gp = gp
    )
  )
}
