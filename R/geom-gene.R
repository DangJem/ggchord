# geom-gene.R - gene arrow layer and gene label layer
# Fetches pre-computed gene arrow polygons / label data from the package
# environment. Gene parameters are specified in this layer and stored for use
# at print time.

# ---------------------------------------------------------------------------
# geom_gene(): gene arrow polygons
# ---------------------------------------------------------------------------

#' Add a gene arrow layer
#'
#' Draws gene annotation arrows on the chord diagram. Gene layout parameters
#' (offset, width, color scheme, etc.) are specified here.
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
#' @param show_legend Whether to show the legend, default TRUE
#' @param legend_position Position of this layer's legend (the Strand or Gene
#'   Annotation legend): one of "left", "right", "top", "bottom" or "inside",
#'   default "right". Pass NULL to let the legend follow
#'   \code{theme(legend.position = ...)} together with the other legends.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A list of ggplot2 layers. To annotate the genes with their labels,
#'   add a \code{\link{geom_gene_label}()} layer.
#' @export
geom_gene <- function(mapping = NULL, data = NULL,
                      gene_offset = NULL,
                      gene_width = NULL,
                      gene_color_scheme = NULL,
                      gene_colors = NULL,
                      gene_order = NULL,
                      show_legend = TRUE,
                      legend_position = "right",
                      ...) {
  # Backward compatibility: gene label parameters used to live here. Point the
  # user to the dedicated layer instead of silently ignoring them.
  legacy_label_args <- intersect(
    names(list(...)),
    c("gene_label_show", "gene_label_size", "gene_label_rotation",
      "gene_label_radial_offset", "gene_label_circum_offset",
      "gene_label_circum_limit", "gene_label_repel", "gene_label_wrap",
      "gene_label_max_overlaps", "gene_label_seed", "show_label", "label_size")
  )
  if (length(legacy_label_args) > 0) {
    warning(
      "Gene label arguments (",
      paste(legacy_label_args, collapse = ", "),
      ") have moved to the dedicated `geom_gene_label()` layer and are ignored ",
      "here. Add `geom_gene_label()` after `geom_gene()` to show gene labels.",
      call. = FALSE
    )
  }

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
    type              = "gene",
    gene_offset       = gene_offset,
    gene_width        = gene_width,
    gene_color_scheme = gene_color_scheme,
    gene_colors       = gene_colors,
    gene_order        = gene_order,
    legend_position   = legend_position
  )
  layers[[length(layers) + 1]] <- poly_layer

  layers
}

# ---------------------------------------------------------------------------
# geom_gene_label(): gene annotation labels
# ---------------------------------------------------------------------------

#' Add a gene label layer
#'
#' Draws the gene annotation labels on a chord diagram. This layer is
#' independent from \code{\link{geom_gene}()}: add it after \code{geom_gene()}
#' to annotate the gene arrows with their texts.
#'
#' Long annotations can be wrapped with \code{gene_label_wrap}. For automatic
#' de-overlapping (with leader lines), use
#' \code{\link{geom_gene_label_repel}()} instead.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param gene_label_size Numeric. Label font size, default 2.5
#' @param gene_label_rotation Optional numeric/vector/list. Label rotation angle, default 0
#' @param gene_label_radial_offset Optional numeric/vector/list. Radial offset of labels, default 0
#' @param gene_label_circum_offset Optional numeric/vector/list. Circumferential offset of labels, default 0
#' @param gene_label_circum_limit Optional logical/vector/list. Whether to limit circumferential offset, default TRUE
#' @param gene_label_wrap Numeric or NULL, default NULL. When set, long gene
#'   annotations are wrapped at this many characters (e.g. 15), which makes the
#'   labels narrower and less prone to overlap.
#' @param show_legend Whether to show the legend, default FALSE
#' @param ... Additional arguments passed to \code{geom_text()}
#'
#' @return A list of ggplot2 layers. To let the labels avoid each other and the
#'   genes (with leader lines), use \code{\link{geom_gene_label_repel}()}
#'   instead.
#' @export
geom_gene_label <- function(mapping = NULL, data = NULL,
                            gene_label_size = NULL,
                            gene_label_rotation = NULL,
                            gene_label_radial_offset = NULL,
                            gene_label_circum_offset = NULL,
                            gene_label_circum_limit = NULL,
                            gene_label_wrap = NULL,
                            show_legend = FALSE,
                            ...) {
  # Placeholder text layer (real data is injected at print time)
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
    show.legend = show_legend,
    ...
  )
  text_layer$ggchord_type <- "gene_text"
  text_layer$ggchord_params <- list(
    type                     = "gene_label",
    gene_label_size          = gene_label_size,
    gene_label_rotation      = gene_label_rotation,
    gene_label_radial_offset = gene_label_radial_offset,
    gene_label_circum_offset = gene_label_circum_offset,
    gene_label_circum_limit  = gene_label_circum_limit,
    gene_label_wrap          = gene_label_wrap
  )
  list(text_layer)
}

# ---------------------------------------------------------------------------
# geom_gene_label_repel(): ggrepel-style gene labels
# ---------------------------------------------------------------------------

#' Add a repelled gene label layer (ggrepel-style)
#'
#' Like \code{\link{geom_gene_label}()}, but the labels are placed with a
#' force-based simulation that pushes them away from the genes and from each
#' other (similar to \code{ggrepel::geom_text_repel()}). Labels that move far
#' enough from their anchor are connected to it with a leader line, and labels
#' that still overlap too many others can be hidden.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param gene_label_size Numeric. Label font size, default 2.5
#' @param gene_label_rotation Optional numeric/vector/list. Label rotation angle, default 0
#' @param gene_label_radial_offset Optional numeric/vector/list. Radial offset of labels, default 0
#' @param gene_label_circum_offset Optional numeric/vector/list. Circumferential offset of labels, default 0
#' @param gene_label_circum_limit Optional logical/vector/list. Whether to limit circumferential offset, default TRUE
#' @param gene_label_wrap Numeric or NULL, default NULL. When set, long gene
#'   annotations are wrapped at this many characters (e.g. 15).
#' @param max_overlaps Numeric, default Inf. Hide labels that still overlap
#'   more than this many other labels after repulsion (ggrepel-style
#'   decluttering). Use a finite value to clean up crowded plots.
#' @param box_padding Numeric, default 0.25. Extra padding around each label
#'   box (data units).
#' @param point_padding Numeric, default 0.1. Extra padding around the anchor
#'   points (data units).
#' @param min_segment_length Numeric, default 0.05. Labels that moved less than
#'   this distance (data units) from their anchor do not draw a leader line.
#'   Keep it small so that every label is connected to its gene.
#' @param force Numeric, default 1. Strength of the repulsive forces.
#' @param seed Numeric, default 123. Random seed for reproducibility.
#' @param gene_label_orientation Character, default "arc". One of
#'   \code{"arc"} (text rotated along the sequence arc) or \code{"horizontal"}
#'   (all labels are drawn horizontally).
#' @param gene_label_segment Character, default "line". Leader line style: a
#'   straight \code{"line"} from the gene to the label, or an L-shaped
#'   \code{"elbow"} (a short segment outward, then a horizontal segment to the
#'   label). Elbow segment lengths adapt to each label's position and text
#'   width, so labels can be placed freely.
#' @param gene_label_side Character, default "auto". Which side of the arc the
#'   labels sit on. \code{"auto"} keeps the strand-based placement
#'   (same as before); \code{"outside"} moves labels that would be inside the
#'   chord (where they can overlap the ribbons) to the outside of their arc;
#'   \code{"inside"} does the opposite. Labels moved to the other side are
#'   connected with a dashed leader line (see \code{gene_label_segment_linetype}).
#' @param gene_label_segment_linetype Character or numeric, default "auto".
#'   Leader-line linetype. \code{"auto"} draws solid lines, except for labels
#'   that were moved to the other side of their arc, which are drawn dashed.
#'   Any other valid ggplot2 linetype (e.g. \code{"solid"}, \code{"dashed"},
#'   \code{"dotted"}, or a numeric dash pattern) is used for all leader lines.
#' @param show_legend Whether to show the legend, default FALSE
#' @param ... Additional arguments passed to \code{geom_text()}
#'
#' @return A list of ggplot2 layers (a leader-line layer and a text layer).
#' @export
geom_gene_label_repel <- function(mapping = NULL, data = NULL,
                                  gene_label_size = NULL,
                                  gene_label_rotation = NULL,
                                  gene_label_radial_offset = NULL,
                                  gene_label_circum_offset = NULL,
                                  gene_label_circum_limit = NULL,
                                  gene_label_wrap = NULL,
                                  max_overlaps = Inf,
                                  box_padding = 0.25,
                                  point_padding = 0.1,
                                  min_segment_length = 0.05,
                                  force = 1,
                                  seed = 123,
                                  gene_label_orientation = "arc",
                                  gene_label_segment = "line",
                                  gene_label_side = "auto",
                                  gene_label_segment_linetype = "auto",
                                  show_legend = FALSE,
                                  ...) {
  gene_label_orientation <- match.arg(gene_label_orientation,
                                      c("arc", "horizontal"))
  gene_label_segment <- match.arg(gene_label_segment, c("line", "elbow"))
  gene_label_side <- match.arg(gene_label_side, c("auto", "inside", "outside"))
  gene_label_segment_linetype <- validate_gene_segment_linetype(
    gene_label_segment_linetype
  )
  layers <- list()

  # Leader line layer (from the anchor to the repelled label position)
  seg_layer <- geom_segment(
    data        = data.frame(x0 = numeric(0), y0 = numeric(0),
                             x1 = numeric(0), y1 = numeric(0),
                             group = integer(0),
                             linetype = character(0)),
    mapping     = aes(x = x0, y = y0, xend = x1, yend = y1, group = group,
                      linetype = I(linetype)),
    inherit.aes = FALSE,
    show.legend = FALSE,
    colour      = "grey50",
    linewidth   = 0.3
  )
  seg_layer$ggchord_type <- "gene_label_segment"
  seg_layer$ggchord_params <- list(type = "gene_label_segment")
  layers[[length(layers) + 1]] <- seg_layer

  # Text layer (drawn at the repelled positions)
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
    show.legend = show_legend,
    ...
  )
  text_layer$ggchord_type <- "gene_text_repel"
  text_layer$ggchord_params <- list(
    type                     = "gene_label_repel",
    gene_label_size          = gene_label_size,
    gene_label_rotation      = gene_label_rotation,
    gene_label_radial_offset = gene_label_radial_offset,
    gene_label_circum_offset = gene_label_circum_offset,
    gene_label_circum_limit  = gene_label_circum_limit,
    gene_label_wrap          = gene_label_wrap,
    max_overlaps             = max_overlaps,
    box_padding              = box_padding,
    point_padding            = point_padding,
    min_segment_length       = min_segment_length,
    force                    = force,
    seed                     = seed,
    gene_label_orientation   = gene_label_orientation,
    gene_label_segment       = gene_label_segment,
    gene_label_side          = gene_label_side,
    gene_label_segment_linetype = gene_label_segment_linetype
  )
  layers[[length(layers) + 1]] <- text_layer

  layers
}
