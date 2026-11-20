#' Draw true aligned character cells
#' @export
#' @param mapping Aesthetic mapping from `ggplot2::aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param representation Draw text, filled tiles, or both.
#' @param overflow Hide text, use tile fallback, or report a failure if the character does not fit.
#' @param line_width Maximum number of displayed alignment columns per wrapped row.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Provide an aligned character table with source or alignment identity. Display wrapping does not change alignment columns, and a gap has no source position. Pair text and tiles explicitly with a residue fill scale when colour conveys alphabet symbols.
#' @examples
#' library(ggplot2)
#'
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), alphabet = "DNA")
#' ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(aes(fill = symbol)) +
#'   scale_fill_residue(alphabet = "DNA") +
#'   ggplot2::coord_fixed(ratio = 1)
geom_residue <- function(mapping = NULL, data = NULL, representation = "tile_text",
    overflow = "tile", line_width = NULL, stat = "identity", position = "identity",
    ..., na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  representation <- match.arg(representation, c("text", "tile", "tile_text"))
  overflow_supplied <- !missing(overflow)
  overflow <- match.arg(overflow, c("tile", "hide", "error"))
  if (!ggchord_identity_stat(stat) || !ggchord_is_data_identity(position))
    ggchord_stop("geom_residue(): use identity Stat and Position")
  if (representation == "tile" && overflow_supplied)
    ggchord_stop("geom_residue(): overflow applies only when text is drawn")
  if (!is.null(line_width)) ggchord_label_integer(line_width, "line_width")
  dots <- ggchord_colour_dots(list(...), "geom_residue()")
  ggchord_residue_validate_styles(c(names(mapping), names(dots)), representation, overflow)
  ggchord_check_geom_params(dots, GeomResidue, "geom_residue()")
  layer <- ggplot2::layer(mapping = mapping, data = data, stat = "identity", geom = GeomResidue,
    position = "identity", show.legend = show.legend, inherit.aes = inherit.aes,
    params = c(list(representation = representation, overflow = overflow,
      line_width = line_width, na.rm = na.rm), dots))
  if (!is.null(id)) layer <- layer_id(layer, id)
  ggplot2::ggproto("LayerResidue", layer,
    compute_aesthetics = function(self, data, plot) {
      mapping <- self$computed_mapping
      fallback <- c(x = "alignment_column", y = "accver", label = "symbol")
      for (field in setdiff(names(fallback), names(mapping)))
        mapping[[field]] <- rlang::new_quosure(rlang::sym(fallback[[field]]), emptyenv())
      fixed <- intersect(names(self$aes_params), names(mapping))
      mapping[fixed] <- NULL
      ggchord_residue_validate_styles(c(names(mapping), names(self$aes_params)),
        representation, overflow)
      self$computed_mapping <- mapping
      value <- ggplot2::ggproto_parent(layer, self)$compute_aesthetics(data, plot)
      if (nrow(value)) {
        if (anyDuplicated(data.frame(x = value$x, y = value$y)))
          ggchord_stop("geom_residue(): one character is required per alignment cell")
        label <- as.character(value$label)
        bad <- !is.na(label) & (nchar(label) != 1L | grepl("[[:space:]]", label))
        if (any(bad))
          ggchord_stop("geom_residue(): each label must contain one residue or gap character")
        if (!is.null(data$is_gap) && !identical(as.logical(data$is_gap), label %in% c("-", ".")))
          ggchord_stop("geom_residue(): symbols conflict with gap metadata")
        if (!is.null(data$is_gap) && !is.null(data$source_position) &&
            any(data$is_gap & !is.na(data$source_position)))
          ggchord_stop("geom_residue(): gaps cannot have source positions")
      }
      if (!is.null(self$geom_params$line_width) && nrow(value)) {
        display <- ggchord_residue_wrap(value$x, value$y, self$geom_params$line_width)
        value$x <- display$column
        value$y <- display$row
        self$ggchord_residue_wrap <- list(origin = display$origin, width = self$geom_params$line_width)
      }
      value
    })
}

ggchord_residue_validate_styles <- function(styles, representation, overflow) {
  if (representation == "tile" &&
      any(c("colour", "size", "family", "fontface", "lineheight") %in% styles))
    ggchord_stop("geom_residue(): tile geometry does not use text styles")
  if (representation == "text" && overflow != "tile" &&
      any(c("fill", "cell_colour", "cell_linewidth") %in% styles))
    ggchord_stop("geom_residue(): text uses cell styles only with overflow = 'tile'")
}

#' Draw true alignment characters
#' @export
#' @param mapping Aesthetic mapping from `ggplot2::aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param overflow Hide text, use tile fallback, or report a failure if the character does not fit.
#' @param line_width Maximum number of displayed alignment columns per wrapped row.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Provide an aligned character table with source or alignment identity. Display wrapping does not change alignment columns, and a gap has no source position. Pair text and tiles explicitly with a residue fill scale when colour conveys alphabet symbols.
#' @examples
#' library(ggplot2)
#'
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), alphabet = "DNA")
#' ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue_text() +
#'   ggplot2::coord_fixed(ratio = 1)
geom_residue_text <- function(mapping = NULL, data = NULL, overflow = "hide",
    line_width = NULL, stat = "identity", position = "identity", ...,
    na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_residue(mapping, data, representation = "text", overflow = overflow,
    line_width = line_width, stat = stat, position = position, ...,
    na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes, id = id)
}

#' Draw true alignment cells as tiles
#' @export
#' @param mapping Aesthetic mapping from `ggplot2::aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param line_width Maximum number of displayed alignment columns per wrapped row.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Provide an aligned character table with source or alignment identity. Display wrapping does not change alignment columns, and a gap has no source position. Pair text and tiles explicitly with a residue fill scale when colour conveys alphabet symbols.
#' @examples
#' library(ggplot2)
#'
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), alphabet = "DNA")
#' ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue_tile(aes(fill = symbol)) +
#'   scale_fill_residue(alphabet = "DNA") +
#'   ggplot2::coord_fixed(ratio = 1)
geom_residue_tile <- function(mapping = NULL, data = NULL, line_width = NULL,
    stat = "identity", position = "identity", ..., na.rm = FALSE,
    show.legend = FALSE, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  if ("overflow" %in% names(list(...))) ggchord_stop("geom_residue_tile(): overflow applies to text")
  geom_residue(mapping, data, representation = "tile", line_width = line_width,
    stat = stat, position = position, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.LayerResidue <- function(object, plot, ...) {
  # Clone before supplying defaults: reusing a layer must not change its mapping.
  original <- object
  object <- ggplot2::ggproto(NULL, original)
  plot <- NextMethod()
  if (!is.null(object$geom_params$line_width)) {
    y_scale <- plot$scales$get_scales("y")
    if (is.null(y_scale)) plot <- plot + ggplot2::scale_y_discrete() else
      if (!y_scale$is_discrete()) ggchord_stop("geom_residue(): line_width requires a discrete y scale")
    if (is.null(plot$scales$get_scales("x"))) {
      width <- object$geom_params$line_width
      plot <- plot + ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .5),
        breaks = function(limits) {
          breaks <- unique(c(1, pretty(c(1, min(width, limits[2L])))))
          breaks[breaks >= 1 & breaks <= width & breaks <= limits[2L]]
        })
    }
    if (is.null(plot$labels$x)) plot <- plot + ggplot2::labs(x = "Column within line")
    if (is.null(plot$labels$y)) plot <- plot + ggplot2::labs(y = NULL)
  }
  plot
}

# Display rows are separate from source/alignment coordinates.
ggchord_residue_wrap <- function(x, y, width, origin = NULL) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
      any(x < 1 | x != floor(x) | x > 2^53 - 1) || length(y) != length(x) || anyNA(y))
    ggchord_stop("geom_residue(): wrapping requires positive integer x positions and nonmissing sequence rows")
  origin <- origin %||% min(x)
  block <- floor((x - origin) / width)
  start <- origin + block * width
  end <- pmin(start + width - 1, max(x))
  label <- paste0(as.character(y), " [", start, "-", end, "]")
  order <- order(block, if (is.factor(y)) as.integer(y) else match(y, unique(y)))
  list(column = x - start + 1, row = factor(label, levels = rev(unique(label[order]))), origin = origin)
}

GeomResidue <- ggplot2::ggproto("GeomResidue", ggplot2::Geom,
  required_aes = c("x", "y", "label"),
  default_aes = ggplot2::aes(colour = "black", fill = "grey95", alpha = 1,
    size = 3, family = "mono", fontface = 1, lineheight = 1.2,
    cell_colour = NA_character_, cell_linewidth = .2),
  extra_params = c("na.rm", "representation", "overflow", "line_width"),
  setup_data = function(data, params) {
    data$xmin <- data$x - .5; data$xmax <- data$x + .5
    data$ymin <- as.numeric(data$y) - .5; data$ymax <- as.numeric(data$y) + .5
    data
  },
  draw_key = ggplot2::draw_key_rect,
  draw_panel = function(data, panel_params, coord, representation = "tile_text",
      overflow = "tile", line_width = NULL, na.rm = FALSE) {
    if (!inherits(coord, "CoordCartesian") || isTRUE(coord$ggchord_coord))
      ggchord_stop("geom_residue(): Cartesian coordinates are required")
    for (axis in c("x", "y")) {
      scale <- panel_params[[axis]]$scale
      if (!is.null(scale) && !scale$is_discrete() &&
          !scale$get_transformation()$name %in% c("identity", "reverse"))
        ggchord_stop("geom_residue(): position scales must be identity or reverse")
    }
    if (!nrow(data)) return(grid::nullGrob())
    grid::gTree(data = data, panel_params = panel_params, coord = coord,
      overflow = overflow, representation = representation, cl = "ggchord_residue_grob")
  })

#' @export
#' @importFrom grid makeContent
makeContent.ggchord_residue_grob <- function(x) {
  d <- x$data
  center <- x$coord$transform(d, x$panel_params)
  corners <- function(sign) {
    value <- d; value$x <- as.numeric(d$x) + sign * .5
    value$y <- as.numeric(d$y) + sign * .5
    x$coord$transform(value, x$panel_params)
  }
  a <- corners(-1); b <- corners(1)
  widths <- abs(b$x - a$x); heights <- abs(b$y - a$y)
  panel_w <- grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
  panel_h <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
  fit <- rep(FALSE, nrow(d))
  if (x$representation != "tile") {
    key <- paste(d$label, d$size, d$family, d$fontface, d$lineheight, sep = "\r")
    metrics <- lapply(which(!duplicated(key)), function(i) {
      glyph <- grid::textGrob(as.character(d$label[i]), gp = grid::gpar(
        fontsize = d$size[i] * 72.27 / 25.4, fontfamily = d$family[i],
        fontface = d$fontface[i], lineheight = d$lineheight[i]))
      c(width = grid::convertWidth(grid::grobWidth(glyph), "mm", valueOnly = TRUE),
        height = grid::convertHeight(grid::grobHeight(glyph), "mm", valueOnly = TRUE))
    })
    metrics <- do.call(rbind, metrics)[match(key, unique(key)), , drop = FALSE]
    fit <- widths * panel_w >= metrics[, "width"] * 1.15 &
      heights * panel_h >= metrics[, "height"] * 1.15
    if (any(!fit) && x$overflow == "error")
      ggchord_stop("geom_residue(): characters do not fit on this device; increase the canvas or reduce text size")
  }
  tile_rows <- if (x$representation %in% c("tile", "tile_text")) rep(TRUE, nrow(d)) else
    if (x$overflow == "tile") !fit else rep(FALSE, nrow(d))
  tiles <- if (!any(tile_rows)) grid::nullGrob() else grid::rectGrob(
    x = center$x[tile_rows], y = center$y[tile_rows],
    width = widths[tile_rows], height = heights[tile_rows], default.units = "native",
    gp = grid::gpar(col = scales::alpha(d$cell_colour[tile_rows], d$alpha[tile_rows]),
      lwd = d$cell_linewidth[tile_rows] * ggplot2::.pt,
      fill = scales::alpha(d$fill[tile_rows], d$alpha[tile_rows])))
  letters <- if (!any(fit)) grid::nullGrob() else grid::textGrob(as.character(d$label[fit]),
    x = center$x[fit], y = center$y[fit], default.units = "native",
    gp = grid::gpar(col = scales::alpha(d$colour[fit], d$alpha[fit]),
      fontsize = d$size[fit] * 72.27 / 25.4, fontfamily = d$family[fit],
      fontface = d$fontface[fit], lineheight = d$lineheight[fit]))
  grid::setChildren(x, grid::gList(tiles, letters))
}
