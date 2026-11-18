# A Cartesian character marker; alignment/source mapping belongs to adapters.

#' Draw physically readable residue cells
#'
#' A standard Cartesian ggplot layer for local nucleotide/protein excerpts.
#' Map `x` to alignment columns (or independently known source positions), `y`
#' to sequence rows and `label` to actual characters. It does not compute an
#' alignment or project through ggchord's biological backbone coordinates.
#' Use [register_ggchord_alignment()] to associate exact source positions with
#' a Fusion/magnify detail panel. Half-cell extents train position scales so
#' first/last cells remain fully visible by default.
#' Native scales must be identity or reverse; nonlinear scales and non-Cartesian
#' coordinates are rejected rather than misrepresenting single-position cells.
#'
#' At drawing time, each cell's final physical width/height is compared with
#' its actual glyph metrics. `overflow = "tile"` omits letters that cannot fit
#' and retains their mapped fill cells; `"error"` requires readable letters.
#' Explicit font sizes are never silently reduced. The decision is remade at every target
#' size, including standard ggsave; DPI alone does not supply physical space.
#' Cells meet at their edges to form a continuous alignment matrix. When the
#' plot still uses its default Cartesian coordinate system, the layer supplies
#' a fixed 1:1 aspect ratio so rows remain as compact as alignment columns.
#' Explicit coordinates take priority. Without a size mapping or fixed size,
#' characters grow uniformly to fit their cells, with a 3 mm readability floor;
#' smaller cells use the same overflow policy rather than shrinking the text.
#' When fill is not supplied, label is mapped to fill. A chemistry scale is
#' inferred from the data's explicit alphabet column (DNA/RNA/AA; otherwise AA).
#' Existing user fill mappings, fixed fills and scales take priority.
#' See [scale_fill_residue()] for explicit alphabet and palette control.
#'
#' @param mapping,data Standard ggplot2 layer mapping and data. Required
#'   aesthetics are x, y and label. Ordinary colour, fill, alpha, size, family
#'   and fontface mappings/scales are supported. Explicit `size` is in mm;
#'   omitted size fits cells uniformly with a 3 mm readability floor.
#' @param overflow `"tile"` (default) or `"error"` when letters do not fit.
#' @param line_width Optional positive integer number of columns per displayed
#'   line. NULL (default) keeps one line. Lines start at the smallest displayed
#'   x value; gaps occupy columns. Multiple sequences share the same breaks.
#'   Row labels show the original column range; the x axis counts within a line.
#'   Wrapping changes display only, preserving original input/source coordinates.
#'   Registered Fusion selections must fit one displayed line; split a selection
#'   crossing lines into explicit source intervals.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Fixed text/cell aesthetics and na.rm, passed to the layer.
#' @return A standard ggplot2 layer. Cells span one x unit and one y unit.
#'   Use discrete y rows or unit-spaced rows.
#' @seealso [read_aligned_fasta()], [map_alignment_position()]
#' @export
#' @examples
#' # Explicit synthetic format example; real-data workflow is in the tutorial.
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), "DNA")
#' ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(aes(fill = symbol))
#' long <- as_ggchord_alignment(c(A = strrep("ACGT", 25)), "DNA")
#' ggplot2::ggplot(long, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(line_width = 50)
geom_residue <- function(mapping = NULL, data = NULL, overflow = c("tile", "error"),
    position = "identity", show.legend = FALSE, inherit.aes = TRUE, line_width = NULL, ...) {
  if (ggchord_local_position(position)) ggchord_stop("geom_residue(): position adjusts alignment cells; use a data position")
  if (!is.null(line_width) && (!is.numeric(line_width) || length(line_width) != 1L ||
      !is.finite(line_width) || line_width < 1 || line_width != floor(line_width) || line_width > 2^53 - 1))
    ggchord_stop("geom_residue(): line_width must be a positive integer or NULL")
  layer <- ggplot2::layer(mapping = mapping, data = data, stat = "identity", geom = GeomResidue,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = c(list(overflow = match.arg(overflow), auto_size = TRUE, line_width = line_width), list(...)))
  ggplot2::ggproto("LayerResidue", layer,
    compute_aesthetics = function(self, data, plot) {
      mapping <- self$computed_mapping
      fixed_fill <- "fill" %in% names(self$aes_params)
      if (!fixed_fill && is.null(mapping$fill) && !is.null(mapping$label)) {
        mapping$fill <- mapping$label
        self$computed_mapping <- mapping
      }
      symbol_fill <- !is.null(mapping$fill) && !is.null(mapping$label) &&
        identical(rlang::quo_get_expr(mapping$fill), rlang::quo_get_expr(mapping$label))
      if (!fixed_fill && symbol_fill && is.null(plot$scales$get_scales("fill"))) {
        alphabet <- unique(as.character(data$alphabet))
        if (length(alphabet) != 1L || !alphabet %in% c("DNA", "RNA", "AA")) alphabet <- "AA"
        plot$scales$add(scale_fill_residue(alphabet))
      }
      self$geom_params$auto_size <- is.null(mapping$size) && !"size" %in% names(self$aes_params)
      # Resolve ordinary mappings/recycling before wrapping and scale training.
      value <- ggplot2::ggproto_parent(layer, self)$compute_aesthetics(data, plot)
      if (!is.null(self$geom_params$line_width) && nrow(value)) {
        if (!all(c("x", "y") %in% names(value)))
          ggchord_stop("geom_residue(): line_width requires input x and y mappings")
        display <- ggchord_residue_wrap(value$x, value$y, self$geom_params$line_width)
        value$x <- display$column
        value$y <- display$row
        self$ggchord_residue_wrap <- list(origin = display$origin, width = self$geom_params$line_width)
      }
      value
    })
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.LayerResidue <- function(object, plot, ...) {
  # Clone before supplying defaults: reusing a layer must not change its mapping.
  original <- object
  object <- ggplot2::ggproto(NULL, original)
  default_coord <- isTRUE(plot$coordinates$default)
  plot <- NextMethod()
  if (default_coord) plot <- plot + ggplot2::coord_fixed(ratio = 1)
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
  default_aes = ggplot2::aes(colour = "#111111", fill = "#E5EDF2", alpha = NA,
    size = 3, family = "sans", fontface = 2, lineheight = 1.2),
  extra_params = c("na.rm", "overflow", "auto_size", "line_width"),
  setup_data = function(data, params) {
    data$xmin <- data$x - .5; data$xmax <- data$x + .5
    data$ymin <- as.numeric(data$y) - .5; data$ymax <- as.numeric(data$y) + .5
    data
  },
  draw_key = ggplot2::draw_key_rect,
  draw_panel = function(data, panel_params, coord, overflow = "tile", auto_size = TRUE, line_width = NULL, na.rm = FALSE) {
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
      overflow = overflow, auto_size = auto_size, cl = "ggchord_residue_grob")
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
  key <- paste(d$label, d$size, d$family, d$fontface, d$lineheight, sep = "\r")
  metrics <- lapply(which(!duplicated(key)), function(i) {
    glyph <- grid::textGrob(as.character(d$label[i]), gp = grid::gpar(
      fontsize = d$size[i] * 72.27 / 25.4, fontfamily = d$family[i],
      fontface = d$fontface[i], lineheight = d$lineheight[i]))
    c(width = grid::convertWidth(grid::grobWidth(glyph), "mm", valueOnly = TRUE),
      height = grid::convertHeight(grid::grobHeight(glyph), "mm", valueOnly = TRUE))
  })
  metrics <- do.call(rbind, metrics)[match(key, unique(key)), , drop = FALSE]
  if (isTRUE(x$auto_size)) {
    # One size for the whole matrix, including gaps, rather than per-letter zoom.
    # Both dimensions constrain the same uniform size.
    growth <- max(1, min(min(widths * panel_w) * .80 / max(metrics[, "width"], .Machine$double.eps),
      min(heights * panel_h) * .80 / max(metrics[, "height"], .Machine$double.eps)))
    d$size <- d$size * growth
    metrics <- metrics * growth
  }
  fit <- widths * panel_w >= metrics[, "width"] * 1.15 &
    heights * panel_h >= metrics[, "height"] * 1.15
  if (any(!fit) && x$overflow == "error") ggchord_stop(
    "geom_residue(): target panel cannot fit ", sum(!fit),
    " character(s) at the requested font size; focus fewer columns, enlarge the canvas or use overflow = 'tile'")
  tiles <- grid::rectGrob(x = center$x, y = center$y, width = widths, height = heights,
    default.units = "native", gp = grid::gpar(col = NA, fill = scales::alpha(d$fill, d$alpha)))
  letters <- if (!any(fit)) grid::nullGrob() else grid::textGrob(as.character(d$label[fit]),
    x = center$x[fit], y = center$y[fit], default.units = "native",
    gp = grid::gpar(col = scales::alpha(d$colour[fit], d$alpha[fit]),
      fontsize = d$size[fit] * 72.27 / 25.4, fontfamily = d$family[fit],
      fontface = d$fontface[fit], lineheight = d$lineheight[fit]))
  grid::setChildren(x, grid::gList(tiles, letters))
}
