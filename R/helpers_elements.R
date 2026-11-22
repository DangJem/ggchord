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
#'
#' @noRd
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
#' @noRd
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
#' @param n,tol Retained internal compatibility arguments. The plasmid cadence
#'   is selected from 1/2/2.5/5 times a power of ten with at most ten ticks.
#' @return Numeric vector of major tick positions. The maximum/seam is omitted.
#' @keywords internal
#' @noRd
breakPointsFunc <- function(max_value, n = 5, tol = 0.5, available_mm = NULL) {
  if (max_value <= 1) return(max_value)
  # Endpoints are scientific positions, not a relabelled zero.
  context <- ggchord_render_context()
  panel <- context$panel_mm
  if (is.null(panel) || any(!is.finite(panel))) panel <- ggchord_render_dimensions() * .7
  budget <- available_mm %||% if (all(is.finite(panel))) min(panel) else 100
  label_mm <- tryCatch(convertWidth(grobWidth(textGrob(
    as.character(max_value), gp=gpar(fontsize=8))), "mm", valueOnly=TRUE),
    error=function(e) 10)
  count <- max(2L, min(as.integer(n) + 2L, floor(budget / max(label_mm + 5, 10))))
  middle <- pretty(c(1,max_value), n=max(1,count-2))
  middle <- middle[middle > 1 & middle < max_value]
  if (length(middle) > count-2) middle <- middle[unique(round(seq(1,length(middle),length.out=count-2)))]
  sort(unique(c(1,middle,max_value)))
}


#' Key glyph for sequence legends
#'
#' Draws the path symbol only when the key data contains colour; otherwise returns a blank (prevents ggplot2 4.x
#' from mixing unrelated layers into other legends with default grey/black symbols).
#' @keywords internal
#' @noRd
key_glyph_seq <- function(data, params, size) {
  data$colour <- data$colour %||% data$seq_colour
  if (is.null(data$colour)) return(ggplot2::zeroGrob())
  col <- ggplot2::alpha(data$colour, data$alpha %||% 1)
  lwd <- (data$linewidth %||% data$size %||% 0.8) * ggplot2::.pt
  grobTree(
    segmentsGrob(
      x0 = unit(0.12, "npc"), x1 = unit(0.74, "npc"),
      y0 = unit(0.5, "npc"), y1 = unit(0.5, "npc"),
      gp = gpar(col = col, lwd = lwd, lineend = "round")
    ),
    polygonGrob(
      x = unit(c(0.70, 0.90, 0.70), "npc"),
      y = unit(c(0.34, 0.5, 0.66), "npc"),
      gp = gpar(fill = col, col = col)
    )
  )
}

#' Key glyph for ribbon legends
#'
#' Draws the polygon symbol only when the key data contains fill; otherwise returns a blank.
#' @keywords internal
#' @noRd
key_glyph_ribbon <- function(data, params, size) {
  data$fill <- data$fill %||% data$ribbon_fill
  data$colour <- data$colour %||% data$ribbon_colour
  data$alpha <- data$alpha %||% data$ribbon_alpha
  data$linetype <- data$linetype %||% data$ribbon_linetype
  if (is.null(data$fill)) return(ggplot2::zeroGrob())
  ggplot2::draw_key_polygon(data, params, size)
}

#' Key glyph for generic genomic features
#' @keywords internal
#' @noRd
key_glyph_feature <- function(data, params, size) {
  key <- params$arrow_key %||% list(style = "flush")
  mapped <- as.character(data$shape)[1L]
  fallback <- data$feature_shape %||% "arrow"
  if (identical(fallback, "arrow")) fallback <- key$shape %||% fallback
  shape <- if (mapped %in% ggchord_feature_shapes) mapped else as.character(fallback)[1L]
  if (!shape %in% ggchord_feature_shapes) shape <- key$shape %||% "arrow"
  strand <- as.character(data$strand %||% "+")[1L]
  opacity <- data$alpha %||% 1
  colour <- scales::alpha(data$colour %||% "#353A3E", opacity)
  fill <- scales::alpha(data$fill %||% data$feature_fill %||% "grey", opacity)
  gp <- gpar(col = colour, fill = fill, lwd = (data$linewidth %||% .25) * ggplot2::.pt,
    lty = data$linetype %||% 1, linejoin = "round")
  if (identical(shape, "lollipop")) {
    return(grobTree(segmentsGrob(.50, .19, .50, .52, gp = gp),
      grid::circleGrob(.50, .61, r = unit(.17, "snpc"), gp = gp)))
  }
  fixed <- isTRUE(key$fixed_body_width) && shape %in% ggchord_block_arrow_shapes
  head_ratio <- key$head_ratio %||%
    if (fixed && identical(key$style, "shouldered")) ggchord_circular_arrow_head_ratio else 1
  parts <- ggchord_feature_geometry(shape, .1, .9, 1, .52, 1, NULL,
    arrow_head_length = .28, arrow_head_width = head_ratio,
    arrow_head_style = key$style, short_feature = "shrink", draw_head = strand != ".",
    draw_start_head = identical(strand, "+/-"), bidirectional = identical(strand, "+/-"),
    fixed_body_width = fixed)
  half_height <- max(vapply(parts, function(part) max(abs(part$radius - 1)), numeric(1)))
  if (identical(shape, "capped_line")) gp$fill <- gp$col
  glyphs <- lapply(parts, function(part) {
    x <- if (identical(strand, "-")) 1 - part$angle else part$angle
    y <- .5 + (part$radius - 1) * min(1, .34 / half_height)
    polygonGrob(unit(x, "npc"), unit(y, "npc"), gp = gp)
  })
  if (length(glyphs) == 1L) glyphs[[1L]] else do.call(grobTree, glyphs)
}
