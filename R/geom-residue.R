# A Cartesian character marker; alignment/source mapping belongs to adapters.

#' Draw physically readable residue cells
#'
#' A standard Cartesian ggplot layer for local nucleotide/protein excerpts.
#' Map `x` to alignment columns (or independently known source positions), `y`
#' to sequence rows and `label` to actual characters. It does not compute an
#' alignment or project through ggchord's biological backbone coordinates.
#' Native scales must be identity or reverse; nonlinear scales and non-Cartesian
#' coordinates are rejected rather than misrepresenting single-position cells.
#'
#' At drawing time, each cell's final physical width/height is compared with
#' its actual glyph metrics. `overflow = "tile"` omits letters that cannot fit
#' and retains their mapped fill cells; `"error"` requires readable letters.
#' Font size is never silently reduced. The decision is remade at every target
#' size, including standard ggsave; DPI alone does not supply physical space.
#'
#' @param mapping,data Standard ggplot2 layer mapping and data. Required
#'   aesthetics are x, y and label. Ordinary colour, fill, alpha, size, family
#'   and fontface mappings/scales are supported. `size` is in mm.
#' @param overflow `"tile"` (default) or `"error"` when letters do not fit.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Fixed text/cell aesthetics and na.rm, passed to the layer.
#' @return A standard ggplot2 layer. Cell width is one x unit and height one y
#'   unit with small internal spacing. Use discrete y rows or unit-spaced rows.
#' @seealso [read_aligned_fasta()], [map_alignment_position()]
#' @export
#' @examples
#' # Explicit synthetic format example; real-data workflow is in the tutorial.
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), "DNA")
#' ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(aes(fill = symbol))
geom_residue <- function(mapping = NULL, data = NULL, overflow = c("tile", "error"),
    position = "identity", show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(mapping = mapping, data = data, stat = "identity", geom = GeomResidue,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = c(list(overflow = match.arg(overflow)), list(...)))
}

GeomResidue <- ggplot2::ggproto("GeomResidue", ggplot2::Geom,
  required_aes = c("x", "y", "label"),
  default_aes = ggplot2::aes(colour = "#202020", fill = "#E5EDF2", alpha = NA,
    size = 3, family = "mono", fontface = 1, lineheight = 1.2),
  extra_params = c("na.rm", "overflow"),
  draw_key = ggplot2::draw_key_label,
  draw_panel = function(data, panel_params, coord, overflow = "tile", na.rm = FALSE) {
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
      overflow = overflow, cl = "ggchord_residue_grob")
  })

#' @export
#' @importFrom grid makeContent
makeContent.ggchord_residue_grob <- function(x) {
  d <- x$data
  center <- x$coord$transform(d, x$panel_params)
  corners <- function(sign) {
    value <- d; value$x <- as.numeric(d$x) + sign * .46
    value$y <- as.numeric(d$y) + sign * .40
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
