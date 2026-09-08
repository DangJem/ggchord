# geom-seq.R - sequence arc layer
# Fetches pre-computed sequence arc data from the package environment and renders it with geom_path
# Sequence layout parameters are specified in this layer and stored for use at print time

GeomChordSeq <- ggplot2::ggproto(
  "GeomChordSeq", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(
    seq_colour = "#3D4348", alpha = 1, linewidth = 0.9,
    linetype = 1, .component = "path"
  ),
  draw_key = ggplot2::draw_key_path,
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        arrow = NULL) {
    grobs <- list()
    paths <- data[data$.component == "path", , drop = FALSE]
    bands <- data[data$.component == "band", , drop = FALSE]
    if (nrow(paths)) {
      paths$colour <- paths$seq_colour
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        paths, panel_params, coord, lineend = "round", arrow = arrow,
        na.rm = na.rm
      )
    }
    if (nrow(bands)) {
      bands$colour <- bands$seq_colour
      bands$fill <- bands$seq_colour
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPolygon$draw_panel(
        bands, panel_params, coord, na.rm = na.rm
      )
    }
    do.call(grid::grobTree, grobs)
  }
)

ggchord_seq_style_geometry <- function(seq_arcs, params) {
  style <- params$seq_style %||% "auto"
  if (identical(style, "auto")) {
    style <- if (isTRUE(params$circular)) "double" else "single"
  }
  gap <- params$seq_backbone_gap %||% 0.025
  width <- params$seq_backbone_width %||% 0.035
  pieces <- lapply(names(seq_arcs), function(id) {
    arc <- seq_arcs[[id]]
    n <- nrow(arc)
    prev <- pmax(1L, seq_len(n) - 1L)
    next_ <- pmin(n, seq_len(n) + 1L)
    dx <- arc$x[next_] - arc$x[prev]
    dy <- arc$y[next_] - arc$y[prev]
    norm <- sqrt(dx^2 + dy^2)
    norm[norm <= 1e-12] <- 1
    nx <- -dy / norm; ny <- dx / norm
    flip <- nx * arc$x + ny * arc$y < 0
    nx[flip] <- -nx[flip]; ny[flip] <- -ny[flip]
    shifted <- function(amount, group, component) {
      out <- arc
      out$x <- out$x + nx * amount
      out$y <- out$y + ny * amount
      out$group <- paste(id, group, sep = "\r")
      out$.component <- component
      out
    }
    switch(
      style,
      single = shifted(0, 1L, "path"),
      double = rbind(
        shifted(gap / 2, 1L, "path"),
        shifted(-gap / 2, 2L, "path")
      ),
      band = {
        outer <- shifted(width / 2, 1L, "band")
        inner <- shifted(-width / 2, 1L, "band")
        out <- rbind(outer, inner[nrow(inner):1L, , drop = FALSE])
        out$group <- paste(id, "band", sep = "\r")
        out
      }
    )
  })
  do.call(rbind, pieces)
}

#' Add a sequence arc layer
#'
#' Draws arcs (or straight lines, depending on the curvature setting) representing sequences in the chord diagram.
#' Sequence layout parameters (order, orientation, radius, curvature, colors, etc.) are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param seq_order Optional character vector. Specifies the drawing order of sequences
#' @param seq_orientation Optional numeric (1 or -1). Sequence orientation, default 1
#' @param seq_gap Optional numeric. Gap proportion between sequences, default 0.03
#' @param seq_radius Optional numeric (> 0). Sequence arc radius, default 1.0
#' @param seq_curvature Optional numeric. Signed arc bow (0=straight, 1=standard arc, negative=opposite bow), default 1.0. Finite positive and negative values are accepted without clipping; magnitudes above 1 amplify the bow.
#' @param seq_style Backbone geometry. `"auto"` uses a double backbone in
#'   [coord_circular()] and the historical single backbone in [coord_chord()].
#' @param seq_backbone_gap Gap between double backbone lines.
#' @param seq_backbone_width Width of a band backbone.
#' @param linewidth Arc line width, default 0.9
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_path()}
#'
#' @return A ggplot2 layer.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq()
#' p
geom_seq <- function(mapping = NULL, data = NULL,
                     seq_order = NULL,
                     seq_orientation = NULL,
                     seq_gap = NULL,
                     seq_radius = NULL,
                     seq_curvature = NULL,
                     seq_style = c("auto", "single", "double", "band"),
                     seq_backbone_gap = 0.025,
                     seq_backbone_width = 0.035,
                     linewidth = 0.9,
                     position = "identity",
                     show.legend = TRUE,
                     inherit.aes = FALSE,
                     ...) {
  show_legend_supplied <- !missing(show.legend)
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  arrow_supplied <- "arrow" %in% names(dots)
  seq_style <- match.arg(seq_style)
  if (!is.numeric(seq_backbone_gap) || length(seq_backbone_gap) != 1L ||
      !is.finite(seq_backbone_gap) || seq_backbone_gap < 0 ||
      !is.numeric(seq_backbone_width) || length(seq_backbone_width) != 1L ||
      !is.finite(seq_backbone_width) || seq_backbone_width <= 0) {
    ggchord_stop("geom_seq(): backbone gap/width must be finite non-negative/positive numbers")
  }
  dots <- ggchord_colour_dots(dots, "geom_seq()")
  ggchord_reject_retired(dots, "geom_seq()", c(
    seq_labels = "geom_seq_label(labels = ...) or a label scale",
    seq_colors = "scale_seq_colour_manual(values = ...)",
    show_legend = "show.legend",
    legend_position = "guides(seq_colour = guide_ggchord_legend(position = ...))"
  ))
  removed_group_args <- intersect(
    names(dots),
    c("seq_group", "seq_group_gap", "seq_group_labels",
      "seq_group_label_radius", "seq_group_colors")
  )
  if (length(removed_group_args) > 0L) {
    ggchord_stop(
      "Sequence grouping was removed in v0.10.0; remove argument(s): ",
      paste(removed_group_args, collapse = ", ")
    )
  }
  removed_group_aes <- intersect(
    names(mapping), c("seq_group", "group_colour")
  )
  if (length(removed_group_aes) > 0L) {
    ggchord_stop(
      "Sequence grouping was removed in v0.10.0; remove aesthetic(s): ",
      paste(removed_group_aes, collapse = ", ")
    )
  }

  # The layout is computed at build time (ggplot_build.ggchord). The
  # parameters are attached to the layer itself so that the plot object is
  # fully self-contained.
  lyr <- ggplot2::layer(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             accver = character(0)),
    mapping     = ggplot2::aes(
      x = x, y = y, group = group, seq_colour = accver,
      .component = I(.component)
    ),
    stat        = "identity",
    geom        = GeomChordSeq,
    position    = position,
    show.legend = if (identical(show.legend, TRUE)) {
                    c(seq_colour = TRUE, fill = FALSE)
                  } else show.legend,
    inherit.aes = inherit.aes,
    check.param = FALSE,
    key_glyph   = key_glyph_seq,
    params      = c(
      list(linewidth = linewidth),
      if (!arrow_supplied) list(
        arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm"))
      ) else list(),
      dots
    )
  )
  lyr$ggchord_type <- "seq"
  lyr$ggchord_params <- list(
    type                  = "seq",
    seq_order             = seq_order,
    seq_labels            = NULL,
    seq_orientation       = seq_orientation,
    seq_gap               = seq_gap,
    seq_radius            = seq_radius,
    seq_curvature         = seq_curvature,
    seq_style             = seq_style,
    seq_backbone_gap      = as.numeric(seq_backbone_gap),
    seq_backbone_width    = as.numeric(seq_backbone_width),
    seq_arrow_supplied    = arrow_supplied,
    seq_legend_supplied   = show_legend_supplied,
    seq_colors            = NULL,
    legend_position       = NULL
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping, c("accver", "length", "seq_ring")
  )
  lyr
}
