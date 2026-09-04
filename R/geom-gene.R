# geom-gene.R - gene arrow layer and gene label layer
# Fetches pre-computed gene arrow polygons / label data from the package
# environment. Gene parameters are specified in this layer and stored for use
# at print time.

# ---------------------------------------------------------------------------
# geom_gene(): gene arrow polygons
# ---------------------------------------------------------------------------

gene_geom <- rename_geom_aes(ggplot2::GeomPolygon, renames = c(fill = "gene_fill"))

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
#' @param position Position adjustment. Use [position_feature_stack()] to place
#'   overlapping genes or features on radial lanes.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A ggplot2 layer. To annotate the genes with their labels,
#'   add a \code{\link{geom_gene_label}()} layer.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(gene_data_example)
#' p <- ggchord(seq_data_example, gene_data = gene_data_example) +
#'   geom_seq() + geom_gene()
#' p
geom_gene <- function(mapping = NULL, data = NULL,
                      gene_offset = NULL,
                      gene_width = NULL,
                      position = "identity",
                      show.legend = TRUE,
                      inherit.aes = FALSE,
                      ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  ggchord_reject_retired(dots, "geom_gene()", c(
    gene_color_scheme = "aes(gene_fill = ...) and scale_gene_fill_manual()",
    gene_colors = "scale_gene_fill_manual(values = ...)",
    gene_order = "scale_gene_fill_manual(limits = ...)",
    show_legend = "show.legend",
    legend_position = "guides(gene_fill = guide_ggchord_legend(position = ...))"
  ))
  dots <- ggchord_colour_dots(dots, "geom_gene()")
  # Backward compatibility: gene label parameters used to live here. Point the
  # user to the dedicated layer instead of silently ignoring them.
  legacy_label_args <- intersect(
    names(dots),
    c("gene_label_show", "gene_label_size", "gene_label_rotation",
      "gene_label_radial_offset", "gene_label_circum_offset",
      "gene_label_circum_limit", "gene_label_repel", "gene_label_wrap",
      "gene_label_max_overlaps", "gene_label_seed", "show_label", "label_size")
  )
  if (length(legacy_label_args) > 0) {
    ggchord_stop(
      "Removed geom_gene() label argument(s): ",
      paste(legacy_label_args, collapse = ", "),
      ". Add geom_gene_label() or geom_gene_label_repel() as a separate layer."
    )
  }

  # Manual colors map by annotation; default colors map by strand.
  gene_scheme <- if (!is.null(mapping) && "gene_fill" %in% names(mapping)) {
    "manual"
  } else {
    "strand"
  }
  fill_mapping <- if (identical(gene_scheme, "manual")) {
    ggplot2::aes(x = x, y = y, group = group, gene_fill = anno)
  } else {
    ggplot2::aes(x = x, y = y, group = group, gene_fill = strand)
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
    geom        = gene_geom,
    position    = position,
    show.legend = if (identical(show.legend, TRUE)) {
                    c(gene_fill = TRUE, colour = FALSE)
                  } else show.legend,
    inherit.aes = inherit.aes,
    check.param = FALSE,
    key_glyph   = key_glyph_gene,
    params      = c(
      if (!("colour" %in% names(dots))) list(colour = "#353A3E") else list(),
      if (!("linewidth" %in% names(dots))) list(linewidth = 0.25) else list(),
      dots
    )
  )
  poly_layer$ggchord_type <- "gene_poly"
  poly_layer$ggchord_params <- list(
    type              = "gene",
    gene_offset       = gene_offset,
    gene_width        = gene_width,
    gene_color_scheme = gene_scheme,
    gene_colors       = NULL,
    gene_order        = NULL,
    legend_position   = NULL
  )
  poly_layer <- ggchord_capture_layer_input(
    poly_layer, data, mapping,
    c("seq_id", "start", "end", "strand", "anno")
  )
  poly_layer
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
#' This is the deterministic, fixed-position label layer. By default labels
#' sit outside their sequence and labels that would collide are omitted in
#' input-row order. Use \code{gene_label_overlap = "nudge"} to retain the
#' legacy gentle adjustment, or \code{"allow"} to draw every label at its
#' requested position. For automatic arrangement with leader lines, use
#' \code{\link{geom_gene_label_repel}()} instead.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param gene_label_orientation Character, default \code{"horizontal"}. Text
#'   orientation relative to each sequence path: \code{"radial"},
#'   \code{"tangent"}, or \code{"horizontal"}.
#' @param gene_label_side Character, default \code{"outside"}. Place labels
#'   outside the chord, use the strand-based \code{"auto"} side, or force
#'   \code{"inside"} placement.
#' @param gene_label_overlap Character, default \code{"hide"}. Fixed-label
#'   collision policy: omit later conflicting labels (\code{"hide"}), apply
#'   the legacy gentle adjustment (\code{"nudge"}), or keep all requested
#'   positions (\code{"allow"}). Input row order therefore provides a simple
#'   way to prioritize labels in \code{"hide"} mode.
#' @param gene_label_rotation Optional numeric/vector/list. Label rotation angle, default 0
#' @param gene_label_radial_offset Optional numeric/vector/list. Radial offset of labels, default 0.04
#' @param gene_label_circum_offset Optional numeric/vector/list. Circumferential offset of labels, default 0
#' @param gene_label_circum_limit Optional logical/vector/list. Whether to limit circumferential offset, default TRUE
#' @param gene_label_wrap Numeric or NULL, default NULL. When set, long gene
#'   annotations are wrapped at this many characters (e.g. 15), which makes the
#'   labels narrower and less prone to overlap.
#' @param position Position adjustment. Supply the same
#'   [position_feature_stack()] used by a gene/feature layer when labels are
#'   drawn from separate data or without a geometry layer.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_text()}
#'
#' @return A ggplot2 layer. To let the labels avoid each other and the
#'   genes (with leader lines), use \code{\link{geom_gene_label_repel}()}
#'   instead.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(gene_data_example)
#' p <- ggchord(seq_data_example, gene_data = gene_data_example) +
#'   geom_seq() + geom_gene() + geom_gene_label()
#' p
geom_gene_label <- function(mapping = NULL, data = NULL,
                            gene_label_orientation = "horizontal",
                            gene_label_side = "outside",
                            gene_label_overlap = "hide",
                            gene_label_rotation = NULL,
                            gene_label_radial_offset = 0.04,
                            gene_label_circum_offset = NULL,
                            gene_label_circum_limit = NULL,
                            gene_label_wrap = NULL,
                            position = "identity",
                            show.legend = FALSE,
                            inherit.aes = FALSE,
                            ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  dots <- ggchord_colour_dots(dots, "geom_gene_label()")
  ggchord_reject_retired(dots, "geom_gene_label()", c(
    gene_label_size = "size",
    show_legend = "show.legend"
  ))

  gene_label_orientation <- match.arg(
    gene_label_orientation, c("radial", "tangent", "horizontal")
  )
  gene_label_side <- match.arg(
    gene_label_side, c("outside", "auto", "inside")
  )
  gene_label_overlap <- match.arg(
    gene_label_overlap, c("hide", "nudge", "allow")
  )

  # Placeholder text layer (real data is injected at print time)
  text_layer <- do.call(ggplot2::geom_text, c(list(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             text_x = numeric(0), text_y = numeric(0),
                             text = character(0), text_angle = numeric(0),
                             hjust = numeric(0), vjust = numeric(0),
                             size = numeric(0)),
    mapping     = ggplot2::aes(x = text_x, y = text_y, label = text,
                      angle = text_angle, hjust = hjust, vjust = vjust,
                      size = I(size)),
    inherit.aes = inherit.aes,
    position = position,
    show.legend = show.legend
  ), dots))
  text_layer$ggchord_type <- "gene_text"
  text_layer$ggchord_theme_element <- "ggchord.gene.label"
  text_layer$ggchord_params <- list(
    type                     = "gene_label",
    gene_label_size          = dots$size %||% NULL,
    gene_label_orientation   = gene_label_orientation,
    gene_label_side          = gene_label_side,
    gene_label_overlap       = gene_label_overlap,
    gene_label_rotation      = gene_label_rotation,
    gene_label_radial_offset = gene_label_radial_offset,
    gene_label_circum_offset = gene_label_circum_offset,
    gene_label_circum_limit  = gene_label_circum_limit,
    gene_label_wrap          = gene_label_wrap
  )
  text_layer <- ggchord_capture_layer_input(
    text_layer, data, mapping,
    c("seq_id", "start", "end", "strand", "anno")
  )
  text_layer
}

# ---------------------------------------------------------------------------
# geom_gene_label_repel(): deterministic automatic gene-label layouts
# ---------------------------------------------------------------------------

#' Add an automatically arranged gene label layer
#'
#' Like \code{\link{geom_gene_label}()}, but labels are placed by one of three
#' deterministic collision-avoiding layouts. The default \code{"aligned"}
#' layout uses orderly cardinal rails, \code{"radial"} uses compact local
#' offset tracks, and \code{"arc"} keeps text close to and rotated with the
#' sequence curve.
#'
#' The local outside/inside label concepts are informed by SnapGene and
#' Geneious, but ggchord uses generic mode names and an independent geometry
#' implementation. See
#' \href{https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map}{SnapGene feature labels}
#' and
#' \href{https://manual.geneious.com/en/latest/Sequences.html}{Geneious label options}.
#'
#' Low-level force, padding, orientation and segment arguments used by earlier
#' releases have been removed and now produce an error. Use
#' \code{gene_label_layout} for automatic placement, or
#' \code{\link{geom_gene_label}()} for manual rotation and offsets.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param gene_label_layout Character, default \code{"aligned"}. Label layout:
#'   \code{"aligned"} uses horizontal labels on orderly top, bottom, left and
#'   right rails; \code{"radial"} uses horizontal labels on the nearest
#'   collision-free local offset track; \code{"arc"} rotates labels along the
#'   sequence tangent and keeps them close to their genes.
#' @param gene_label_wrap Numeric or NULL, default NULL. When set, long gene
#'   annotations are wrapped at this many characters (e.g. 15).
#' @param max_overlaps Numeric, default Inf. Hide labels that still overlap
#'   more than this many other labels after repulsion (ggrepel-style
#'   decluttering). Use a finite value to clean up crowded plots.
#' @param gene_label_side Character, default "outside". Which side of the arc
#'   the labels sit on. \code{"auto"} keeps the strand-based placement
#'   (same as before); \code{"outside"} moves labels that would be inside the
#'   chord (where they can overlap the ribbons) to the outside of their arc;
#'   \code{"inside"} does the opposite. Labels moved to the other side are
#'   connected with a dashed leader line (see \code{gene_label_segment_linetype}).
#' @param gene_label_segment_linetype Character or numeric, default "auto".
#'   Leader-line linetype. \code{"auto"} draws solid lines, except for labels
#'   that were moved to the other side of their arc, which are drawn dashed.
#'   Any other valid ggplot2 linetype (e.g. \code{"solid"}, \code{"dashed"},
#'   \code{"dotted"}, or a numeric dash pattern) is used for all leader lines.
#' @param position Position adjustment. Supply the same
#'   [position_feature_stack()] used by a gene/feature layer when labels are
#'   drawn from separate data or without a geometry layer.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_text()}
#'
#' @return A ggplot2 layer that draws labels and leader lines together.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(gene_data_example)
#' p <- ggchord(seq_data_example, gene_data = gene_data_example) +
#'   geom_seq() + geom_gene() + geom_gene_label_repel()
#' p
geom_gene_label_repel <- function(mapping = NULL, data = NULL,
                                  gene_label_layout = "aligned",
                                  gene_label_wrap = NULL,
                                  gene_label_side = "outside",
                                  max_overlaps = Inf,
                                  gene_label_segment_linetype = "auto",
                                  position = "identity",
                                  show.legend = FALSE,
                                  inherit.aes = FALSE,
                                  ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  # Inspect the unevaluated call before R's partial argument matching can turn
  # the removed `gene_label_segment` into `gene_label_segment_linetype`.
  raw_argument_names <- names(as.list(sys.call())[-1L])
  dots <- list(...)
  dots <- ggchord_colour_dots(dots, "geom_gene_label_repel()")
  removed <- c(
    "gene_label_rotation", "gene_label_radial_offset",
    "gene_label_circum_offset", "gene_label_circum_limit",
    "box_padding", "point_padding", "min_segment_length", "force", "seed",
    "gene_label_orientation", "gene_label_segment"
  )
  supplied_removed <- intersect(c(names(dots), raw_argument_names), removed)
  if (length(supplied_removed) > 0) {
    positioning <- intersect(
      supplied_removed,
      c("gene_label_rotation", "gene_label_radial_offset",
        "gene_label_circum_offset", "gene_label_circum_limit")
    )
    replacement <- if (length(positioning) > 0) {
      " Use geom_gene_label() when manual rotation or offsets are required."
    } else {
      " Repulsion and leader geometry are now selected by gene_label_layout."
    }
    ggchord_stop(
      "Removed geom_gene_label_repel() argument(s): ",
      paste(supplied_removed, collapse = ", "), ".", replacement
    )
  }
  ggchord_reject_retired(dots, "geom_gene_label_repel()", c(
    gene_label_size = "size",
    show_legend = "show.legend"
  ))

  gene_label_layout <- match.arg(
    gene_label_layout, c("aligned", "radial", "arc")
  )
  gene_label_side <- match.arg(gene_label_side, c("auto", "inside", "outside"))
  gene_label_segment_linetype <- validate_gene_segment_linetype(
    gene_label_segment_linetype
  )
  empty <- data.frame(
    x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
    label = character(), .component = character(), group = integer(),
    linetype = character(), size = numeric(), angle = numeric(),
    hjust = numeric(), vjust = numeric()
  )
  lyr <- ggplot2::layer(
    data = empty,
    mapping = ggplot2::aes(
      x = x, y = y, xend = xend, yend = yend, label = label,
      group = group, linetype = I(linetype), size = I(size),
      angle = angle, hjust = hjust, vjust = vjust,
      .component = I(.component)
    ),
    stat = "identity",
    geom = GeomChordGeneLabelRepel,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    check.aes = FALSE,
    check.param = FALSE,
    params = list(
      segment_params = list(),
      text_params = dots,
      na.rm = dots$na.rm %||% FALSE
    )
  )
  lyr$ggchord_type <- "gene_label_repel"
  lyr$ggchord_theme_components <- c(
    segment_params = "ggchord.gene.label.segment",
    text_params = "ggchord.gene.label"
  )
  lyr$ggchord_params <- list(
    type                     = "gene_label_repel",
    gene_label_layout        = gene_label_layout,
    gene_label_size          = dots$size %||% NULL,
    gene_label_wrap          = gene_label_wrap,
    max_overlaps             = max_overlaps,
    gene_label_side          = gene_label_side,
    gene_label_segment_linetype = gene_label_segment_linetype
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping,
    c("seq_id", "start", "end", "strand", "anno")
  )
  lyr
}
