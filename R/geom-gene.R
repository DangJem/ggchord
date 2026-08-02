# geom-gene.R - gene arrow layer
# Fetches pre-computed gene arrow polygons and label data from the package environment
# Gene parameters are specified in this layer and stored for use at print time

#' Add a gene arrow layer
#'
#' Draws gene annotation arrows on the chord diagram. Gene layout parameters (offset, width, color scheme, etc.) are specified here.
#' The gene fill scale is kept independent from the ribbon's fill scale via a
#' separate internal aesthetic used by the ribbon layer.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param gene_offset Optional numeric/vector/list. Radial offset of gene arrows, default 0.1
#' @param gene_width Optional numeric/vector/list. Width of gene arrows, default 0.05
#' @param gene_color_scheme Character. "strand" or "manual", default "strand"
#' @param gene_colors Optional color vector. Fill color of gene arrows
#' @param gene_order Optional character vector. Display order of genes in the legend
#' @param gene_label_show Logical. Whether to show gene labels, default FALSE
#' @param gene_label_size Numeric. Label font size, default 2.5
#' @param gene_label_rotation Optional numeric/vector/list. Label rotation angle, default 0
#' @param gene_label_radial_offset Optional numeric/vector/list. Radial offset of labels, default 0
#' @param gene_label_circum_offset Optional numeric/vector/list. Circumferential offset of labels, default 0
#' @param gene_label_circum_limit Optional logical/vector/list. Whether to limit circumferential offset, default TRUE
#' @param show_legend Whether to show the legend, default TRUE
#' @param show_label Whether to show gene labels (overrides gene_label_show), default NULL
#' @param label_size Label font size (overrides gene_label_size), default NULL
#' @param legend_position Position of this layer's legend (the Strand or Gene
#'   Annotation legend): one of "left", "right", "top", "bottom" or "inside",
#'   default "right". Pass NULL to let the legend follow
#'   \code{theme(legend.position = ...)} together with the other legends.
#' @param gene_label_repel Logical, default FALSE. When TRUE, overlapping gene
#'   labels are automatically pushed apart (collision detection + automatic
#'   avoidance) so they do not cover each other.
#' @param gene_label_wrap Numeric or NULL, default NULL. When set, long gene
#'   annotations are wrapped at this many characters (e.g. 15), which makes the
#'   labels narrower and less prone to overlap.
#' @param gene_label_max_overlaps Numeric, default Inf. With
#'   \code{gene_label_repel = TRUE}, labels that still overlap more than this
#'   many other labels after de-overlapping are hidden (ggrepel-style). Use a
#'   finite value to declutter crowded plots.
#' @param gene_label_seed Numeric, default 123. Seed used by the de-overlap
#'   algorithm for reproducible results.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A list of ggplot2 layers
#' @export
geom_gene <- function(mapping = NULL, data = NULL,
                      gene_offset = NULL,
                      gene_width = NULL,
                      gene_color_scheme = NULL,
                      gene_colors = NULL,
                      gene_order = NULL,
                      gene_label_show = NULL,
                      gene_label_size = NULL,
                      gene_label_rotation = NULL,
                      gene_label_radial_offset = NULL,
                      gene_label_circum_offset = NULL,
                      gene_label_circum_limit = NULL,
                      show_legend = TRUE,
                      show_label = NULL,
                      label_size = NULL,
                      legend_position = "right",
                      gene_label_repel = FALSE,
                      gene_label_wrap = NULL,
                      gene_label_max_overlaps = Inf,
                      gene_label_seed = 123,
                      ...) {
  layers <- list()

  # Manual colors map by annotation; default colors map by strand.
  gene_scheme <- gene_color_scheme %||% "strand"
  fill_mapping <- if (identical(gene_scheme, "manual")) {
    aes(x = x, y = y, group = group, fill = anno)
  } else {
    aes(x = x, y = y, group = group, fill = strand)
  }

  # The polygon layer carries the gene parameters so that the plot object is
  # self-contained (scales are added at build time by ggplot_build.ggchord).
  poly_layer <- ggplot2::layer(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             group = integer(0),
                             strand = factor(character(0), levels = c("+", "-")),
                             anno = character(0), ord = integer(0)),
    mapping     = fill_mapping,
    stat        = "identity",
    geom        = GeomPolygon,
    position    = "identity",
    show.legend = if (identical(show_legend, TRUE)) {
                    c(fill = TRUE, colour = FALSE)
                  } else show_legend,
    inherit.aes = FALSE,
    check.param = FALSE,
    key_glyph   = key_glyph_gene,
    params      = list(color = "black", ...)
  )
  poly_layer$ggchord_type <- "gene_poly"
  poly_layer$ggchord_params <- list(
    type                      = "gene",
    gene_offset               = gene_offset,
    gene_width                = gene_width,
    gene_color_scheme         = gene_color_scheme,
    gene_colors               = gene_colors,
    gene_order                = gene_order,
    gene_label_show           = gene_label_show,
    gene_label_size           = gene_label_size,
    gene_label_rotation       = gene_label_rotation,
    gene_label_radial_offset  = gene_label_radial_offset,
    gene_label_circum_offset  = gene_label_circum_offset,
    gene_label_circum_limit   = gene_label_circum_limit,
    show_label_override       = show_label,
    label_size_override       = label_size,
    legend_position           = legend_position,
    gene_label_repel          = gene_label_repel,
    gene_label_wrap           = gene_label_wrap,
    gene_label_max_overlaps   = gene_label_max_overlaps,
    gene_label_seed           = gene_label_seed
  )
  layers[[length(layers) + 1]] <- poly_layer

  # Placeholder Text layer (for gene label injection; real data is injected at print time)
  text_layer <- geom_text(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             text_x = numeric(0), text_y = numeric(0),
                             text = character(0), text_angle = numeric(0),
                             hjust = numeric(0), vjust = numeric(0),
                             size = numeric(0)),
    mapping     = aes(x = text_x, y = text_y, label = text,
                      angle = text_angle, hjust = hjust, vjust = vjust,
                      size = size),
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  text_layer$ggchord_type <- "gene_text"
  layers[[length(layers) + 1]] <- text_layer

  layers
}
