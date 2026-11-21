# geom-seq.R - sequence arc layer
# Fetches pre-computed sequence arc data from the package environment and renders it with geom_path
# Backbone appearance is owned by each sequence layer.

GeomChordSeq <- ggproto(
  "GeomChordSeq", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = aes(
    colour = "#3D4348", fill = "grey60", alpha = 1, linewidth = 0.9,
    linetype = 1, .component = "path"
  ),
  draw_key = ggplot2::draw_key_path,
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        arrow = NULL) {
    grobs <- list()
    paths <- data[data$.component == "path", , drop = FALSE]
    bands <- data[data$.component == "band", , drop = FALSE]
    breaks <- data[data$.component == "break", , drop = FALSE]
    if (nrow(paths)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        paths, panel_params, coord, lineend = "round", arrow = arrow,
        na.rm = na.rm
      )
    }
    if (nrow(bands)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPolygon$draw_panel(
        bands, panel_params, coord
      )
    }
    if (nrow(breaks)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        breaks, panel_params, coord, lineend = "butt", arrow = NULL,
        na.rm = na.rm)
    }
    do.call(grobTree, grobs)
  }
)

ggchord_seq_style_geometry <- function(seq_arcs, params) {
  style <- params$seq_style
  if (!length(seq_arcs)) return(data.frame(x = numeric(), y = numeric(),
    accver = character(), group = character(), .component = character()))
  ids <- names(seq_arcs)
  physical <- params$seq_thickness_mm
  if (!is.null(names(physical)) && (anyDuplicated(names(physical)) ||
      any(!names(physical) %in% (params$geometry_reference_ids %||% ids))))
    ggchord_stop("geom_seq(): invalid or unknown thickness keys")
  thickness <- process_sequence_param(physical, ids, "thickness", 1)
  offset <- ggchord_resolve_position_offset(params$feature_position$offset %||% 0,
    params$geometry_reference_ids %||% ids, "offset", "geom_seq()")
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
    if (all(abs(dy) < 1e-12)) {
      nx[] <- 0; ny[] <- 1
    } else {
      flip <- nx * arc$x + ny * arc$y < 0
      nx[flip] <- -nx[flip]; ny[flip] <- -ny[flip]
    }
    physical_scale <- ggchord_units_per_mm(params$normal_axis)
    width <- thickness[[id]] * physical_scale
    shift <- offset[[id]][[1L]] * physical_scale
    shifted <- function(amount, group, component) {
      out <- arc
      out$x <- out$x + nx * (amount + shift)
      out$y <- out$y + ny * (amount + shift)
      out$group <- paste(id, group, sep = "\r")
      out$.component <- component
      out
    }
    switch(
      style,
      single = shifted(0, 1L, "path"),
      double = rbind(
        shifted(width / 2, 1L, "path"),
        shifted(-width / 2, 2L, "path")
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

#' Draw sequence backbones
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional reference table or function selecting from the registered reference.
#' @param geometry Backbone style, either `"line"`, `"double"` or `"band"`.
#' @param thickness Physical thickness in mm or absolute grid units. Named values select references.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details The source is the complete reference table with `accver` and `length`. A backbone can select an authoritative subset. An explicit complete table in geom_seq() can register a reference when the plot has no global data.
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' ggchord() +
#'   geom_seq(data = seqs, colour = "#187B77") +
#'   geom_gene(data = genes) +
#'   coord_collinear()
#'
#' ggchord(seqs) +
#'   geom_seq(linewidth = 2, colour = "white") +
#'   geom_seq(linewidth = 0.5, colour = "black") +
#'   coord_collinear()
geom_seq <- function(mapping = NULL, data = NULL, geometry = "line",
    thickness = 1, stat = "identity", position = "identity", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geometry <- match.arg(geometry, c("line", "double", "band"))
  if (!ggchord_identity_stat(stat)) ggchord_stop("geom_seq(): stat must be identity")
  position <- ggchord_as_feature_position(position, "geom_seq()")
  if (ggchord_local_position(position) && !identical(position$placement, "offset"))
    ggchord_stop("geom_seq(): position must be identity or position_offset()")
  if (geometry == "line" && !missing(thickness))
    ggchord_stop("geom_seq(): line uses linewidth, thickness applies to double or band")
  thickness <- ggchord_physical_mm(thickness, "geom_seq(thickness)", nonnegative = TRUE)
  if (any(unlist(thickness, use.names = FALSE) <= 0))
    ggchord_stop("geom_seq(): thickness must be positive")
  dots <- ggchord_colour_dots(list(...), "geom_seq()")
  if (geometry != "band" && "fill" %in% c(names(mapping), names(dots)))
    ggchord_stop("geom_seq(): fill applies only to band")
  if (any(c("accver", "length", "sequence_unit", "topology", "seq_ring") %in% names(mapping)))
    ggchord_stop("geom_seq(): reference identity is defined by the table")
  ggchord_check_geom_params(dots, GeomChordSeq, "geom_seq()")
  layer <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), group = character()),
    mapping = aes(x = x, y = y, group = group, .component = I(.component)),
    stat = "identity", geom = GeomChordSeq, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, check.param = FALSE,
    params = c(list(na.rm = na.rm), dots)
  )
  layer$ggchord_type <- "seq"
  layer$ggchord_params <- list(id = id, type = "seq", physical_layout = geometry != "line" ||
    ggchord_local_position(position),
    seq_style = if (geometry == "line") "single" else geometry,
    seq_backbone_gap = 0, seq_backbone_width = 0,
    feature_position = position, seq_thickness_mm = thickness)
  ggchord_capture_layer_input(layer, data, mapping, character())
}
