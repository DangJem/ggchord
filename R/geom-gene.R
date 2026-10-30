# geom-gene.R - gene arrow layer and gene label layer
# The gene shortcut shares feature construction; layer-owned inputs and
# parameters are resolved by the plot-owned layout calculation at build time.

# ---------------------------------------------------------------------------
# geom_gene(): gene arrow polygons
# ---------------------------------------------------------------------------

gene_geom <- rename_geom_aes(GeomChordFeatureBase, renames = c(fill = "gene_fill"))

#' Add a gene arrow layer
#'
#' A directed-gene shortcut using the generic \code{\link{geom_feature}()}
#' constructor with arrow geometry and the independent \code{gene_fill} role.
#' Gene input still requires a known \code{+} or \code{-} strand; use
#' \code{\link{geom_feature}()} for unknown strands.
#' The gene fill scale is kept independent from the ribbon's fill scale via a
#' separate internal aesthetic used by the ribbon layer.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param width Optional numeric/vector/list. Width of gene arrows, default 0.05
#' @param arrow_head_length,arrow_head_width Shared arrow-head dimensions in
#'   the sequence-local frame.
#' @param arrow_head_style Arrowhead construction: shouldered, flush, or a
#'   full triangular wedge.
#' @param short_feature Fallback for arrows too short to hold their requested
#'   head: automatic wedge/block selection, a wedge, or a block.
#' @param position Feature placement. Use `"identity"`, `"strand"`,
#'   `"plasmid"`, [position_strand()], [position_plasmid()], or
#'   [position_feature_stack()].
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
#'   geom_seq() + geom_gene(width = .05)
#' p
geom_gene <- function(mapping = NULL, data = NULL,
                      width = NULL,
                      arrow_head_length = 0.04, arrow_head_width = 1,
                      arrow_head_style = c("shouldered", "flush", "triangle"),
                      short_feature = c("auto", "wedge", "block"),
                      position = "identity", show.legend = TRUE,
                      inherit.aes = FALSE, ...) {
  dots <- ggchord_colour_dots(list(...), "geom_gene()")
  ggchord_reject_retired(dots, "geom_gene()", c(
    gene_width = "width",
    gene_offset = "position = position_strand(offset = ...)",
    gene_color_scheme = "aes(gene_fill = ...) and scale_gene_fill_manual()",
    gene_colors = "scale_gene_fill_manual(values = ...)",
    gene_order = "scale_gene_fill_manual(limits = ...)",
    show_legend = "show.legend",
    legend_position = "guides(gene_fill = guide_ggchord_legend(position = ...))"
  ))
  labels <- names(dots)[grepl("^gene_label_|^label_", names(dots)) |
    names(dots) == "show_label"]
  if (length(labels)) ggchord_stop(
    "Removed geom_gene() label arguments: ", paste(labels, collapse = ", "),
    ". Add geom_gene_label() or geom_gene_label_repel().")
  mapping <- ggchord_normalize_mapping(mapping)
  feature_mapping <- mapping
  if ("gene_fill" %in% names(feature_mapping)) {
    names(feature_mapping)[names(feature_mapping) == "gene_fill"] <- "feature_fill"
  }
  if ("gene_fill" %in% names(dots)) {
    names(dots)[names(dots) == "gene_fill"] <- "feature_fill"
  }
  layer <- do.call(geom_feature, c(list(
    mapping = feature_mapping, data = data, shape = "arrow", width = width,
    arrow_head_length = arrow_head_length,
    arrow_head_width = arrow_head_width, arrow_head_style = match.arg(arrow_head_style),
    short_feature = match.arg(short_feature), position = position,
    show.legend = show.legend, inherit.aes = inherit.aes), dots))
  # A directed-gene shortcut retains its independent role scale and strict
  # biological input contract while sharing the generic feature constructor.
  layer$geom <- gene_geom
  layer$geom_params[c("segment_boundaries", "boundary_colour", "boundary_linewidth",
    "boundary_linetype", "adaptive_outline")] <- NULL
  layer$aes_params$feature_shape <- NULL
  if ("feature_fill" %in% names(layer$aes_params)) {
    names(layer$aes_params)[names(layer$aes_params) == "feature_fill"] <- "gene_fill"
  }
  mapped <- "gene_fill" %in% names(mapping)
  layer$mapping <- if (mapped) {
    ggplot2::aes(x = x, y = y, group = group, gene_fill = anno)
  } else ggplot2::aes(x = x, y = y, group = group, gene_fill = strand)
  layer$show.legend <- if (identical(show.legend, TRUE)) {
    c(gene_fill = TRUE, colour = FALSE)
  } else show.legend
  layer$ggchord_params$is_feature <- FALSE
  layer$ggchord_params$feature_role <- "gene"
  layer$ggchord_params$gene_color_scheme <- if (mapped) "manual" else "strand"
  layer$ggchord_input_transform <- NULL
  ggchord_capture_layer_input(layer, data, mapping,
    c("accver", "start", "end", "strand", "anno"))
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
#' input-row order. Use \code{label_overlap = "nudge"} to retain the
#' legacy gentle adjustment, or \code{"allow"} to draw every label at its
#' requested position. For automatic arrangement with leader lines, use
#' \code{\link{geom_gene_label_repel}()} instead.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param label_orientation Character, default \code{"horizontal"}. Text
#'   orientation relative to each sequence path: \code{"radial"},
#'   \code{"tangent"}, or \code{"horizontal"}.
#' @param label_side Character, default \code{"outside"}. Place labels
#'   outside the chord, use the strand-based \code{"auto"} side, or force
#'   \code{"inside"} placement.
#' @param label_overlap Character, default \code{"hide"}. Fixed-label
#'   collision policy: omit later conflicting labels (\code{"hide"}), apply
#'   the legacy gentle adjustment (\code{"nudge"}), or keep all requested
#'   positions (\code{"allow"}). Input row order therefore provides a simple
#'   way to prioritize labels in \code{"hide"} mode.
#' @param label_rotation Optional numeric/vector/list. Label rotation angle, default 0
#' @param label_radial_offset Optional numeric/vector/list. Radial offset of labels, default 0.04
#' @param label_circum_offset Optional numeric/vector/list. Circumferential offset of labels, default 0
#' @param label_circum_limit Optional logical/vector/list. Whether to limit circumferential offset, default TRUE
#' @param label_wrap Numeric or NULL, default NULL. When set, long gene
#'   annotations are wrapped at this many characters (e.g. 15), which makes the
#'   labels narrower and less prone to overlap.
#' @param position Position adjustment. Supply the same
#'   [position_feature_stack()] used by a gene/feature layer when labels are
#'   drawn from separate data or without a geometry layer.
#'   In \code{coord_collinear()}, labels without an explicit Position share
#'   a matching geometry layer's lanes; explicit label Positions do not move
#'   the geometry layer.
#'   Unique `feature_id` columns allow label subsets and reordered labels when
#'   their accession, interval and strand agree with the geometry input.
#'   Multiple matching geometry layers require an explicit label Position or
#'   unambiguous feature IDs.
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
#'   geom_seq() + geom_gene() + geom_gene_label(label_side = "auto")
#' p
geom_gene_label <- function(mapping = NULL, data = NULL,
    label_orientation = "horizontal", label_side = "outside",
    label_overlap = "hide", label_rotation = NULL,
    label_radial_offset = 0.04, label_circum_offset = NULL,
    label_circum_limit = NULL, label_wrap = NULL,
    position = "identity", show.legend = FALSE, inherit.aes = FALSE, ...) {
  dots <- list(...)
  ggchord_reject_label_prefix(dots, "geom_gene_label()")
  args <- c(list(mapping = mapping, data = data,
    gene_label_orientation = label_orientation, gene_label_side = label_side,
    gene_label_overlap = label_overlap, gene_label_rotation = label_rotation,
    gene_label_radial_offset = label_radial_offset,
    gene_label_circum_offset = label_circum_offset,
    gene_label_circum_limit = label_circum_limit, gene_label_wrap = label_wrap,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes), dots)
  if (missing(position)) args$position <- NULL
  do.call(ggchord_label_layer, args)
}

#' @noRd
ggchord_label_layer <- function(mapping = NULL, data = NULL,
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
  position_missing <- missing(position)
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  position_object <- ggchord_as_feature_position(position, "geom_gene_label()")
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
    position = position_object,
    show.legend = show.legend
  ), dots))
  text_layer$ggchord_type <- "gene_text"
  text_layer$ggchord_theme_element <- "ggchord.gene.label"
  text_layer$ggchord_params <- list(
    type                     = "gene_label",
    gene_label_size          = dots$size %||% NULL,
    gene_label_family        = dots$family %||% "",
    gene_label_fontface      = dots$fontface %||% 1,
    gene_label_lineheight    = dots$lineheight %||% 1.2,
    gene_label_orientation   = gene_label_orientation,
    gene_label_side          = gene_label_side,
    gene_label_overlap       = gene_label_overlap,
    gene_label_rotation      = gene_label_rotation,
    gene_label_radial_offset = gene_label_radial_offset,
    gene_label_circum_offset = gene_label_circum_offset,
    gene_label_circum_limit  = gene_label_circum_limit,
    gene_label_wrap          = gene_label_wrap,
    feature_position         = position_object,
    position_supplied        = !position_missing
  )
  text_layer <- ggchord_capture_layer_input(
    text_layer, data, mapping,
    c("accver", "start", "end", "strand", "anno")
  )
  text_layer
}

# ---------------------------------------------------------------------------
# geom_gene_label_repel(): deterministic automatic gene-label layouts
# ---------------------------------------------------------------------------

#' Add an automatically arranged gene label layer
#'
#' Like \code{\link{geom_gene_label}()}, with deterministic external label
#' placement. The default \code{"radial"} keeps horizontal text on a contour
#' following each sequence; crowded names fan out along that contour.
#' \code{"auto"} aligns left/right names in columns and uses radial placement
#' elsewhere. \code{"arc"} rotates text along the sequence tangent.
#'
#' The local outside/inside label organisation adapts concepts used by
#' SnapGene and Geneious. ggchord keeps generic mode names and integrates those
#' ideas into its own geometry system. See
#' \href{https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map}{SnapGene feature labels}
#' and
#' \href{https://manual.geneious.com/en/latest/Sequences.html}{Geneious label options}.
#'
#' Low-level force, padding, orientation and segment arguments used by earlier
#' releases have been removed and now produce an error. Use
#' \code{label_layout} for automatic placement, or
#' \code{\link{geom_gene_label}()} for manual rotation and offsets.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param label_layout Character, default \code{"radial"}. Label layout:
#'   \code{"radial"} arranges horizontal external labels along sequence-relative
#'   contours, with straight normal leaders or a short normal departure followed
#'   by a long connector when displaced; \code{"auto"} uses shared vertical
#'   columns on the left/right and radial placement elsewhere; \code{"arc"}
#'   rotates text along the sequence tangent. The former \code{"aligned"}
#'   mode has been removed; use \code{"auto"} for side columns.
#' @param label_wrap Numeric or NULL, default NULL. When set, long gene
#'   annotations are wrapped at this many characters (e.g. 15), overriding
#'   \code{label_fit}.
#' @param label_fit Character, default \code{"wrap"}. How labels whose
#'   initial boxes are crowded or unusually wide are fitted before deterministic layout:
#'   \code{"wrap"} preserves the complete annotation over multiple lines,
#'   \code{"ellipsis"} uses a single shortened line, \code{"auto"} wraps and
#'   adds an ellipsis only if text exceeds \code{label_max_lines}, and
#'   \code{"none"} disables adaptive fitting. Explicit \code{"shrink"} uses
#'   internal-first feature placement and reduces the complete text only as
#'   far as \code{label_min_size}. Labels that still cannot fit use the
#'   existing adjacent/external placement or hiding policy, at their requested
#'   size. This choice uses feature placement rather than the external layout
#'   selected by \code{label_layout}; other fitting modes are unchanged.
#' @param label_min_size Minimum text size in millimetres for explicit
#'   \code{label_fit = "shrink"}, default 1.5. It must not exceed the
#'   requested text size. Source text, requested/final size and fitting reason
#'   are retained in exported label geometry.
#' @param label_max_lines Positive integer, default 2. Maximum number of
#'   lines used by adaptive \code{"wrap"} or \code{"auto"} fitting.
#' @param max_overlaps Numeric, default Inf. Hide labels that still overlap
#'   more than this many other labels after repulsion (ggrepel-style
#'   decluttering). Use a finite value to clean up crowded plots.
#' @param label_side Character, default "outside". Which side of the arc
#'   the labels sit on. \code{"auto"} keeps the strand-based placement
#'   (same as before); \code{"outside"} moves labels that would be inside the
#'   chord (where they can overlap the ribbons) to the outside of their arc;
#'   \code{"inside"} does the opposite. Labels moved to the other side are
#'   connected with a dashed leader line (see \code{label_segment_linetype}).
#' @param label_segment_linetype Character or numeric, default "auto".
#'   Leader-line linetype. \code{"auto"} draws solid lines, except for labels
#'   that were moved to the other side of their arc, which are drawn dashed.
#'   Any other valid ggplot2 linetype (e.g. \code{"solid"}, \code{"dashed"},
#'   \code{"dotted"}, or a numeric dash pattern) is used for all leader lines.
#' @param label_segment_overlap Character, default \code{"fade"}.
#'   How portions of a leader line covered by another label are drawn:
#'   \code{"fade"} keeps them at reduced opacity, \code{"clip"} removes them,
#'   and \code{"show"} leaves the complete leader visible. This does not affect
#'   text wrapping, ellipsis or label hiding.
#' @param label_segment_overlap_alpha Numeric in \code{[0, 1]}, default
#'   0.18. Relative opacity applied only to covered leader-line portions when
#'   \code{label_segment_overlap = "fade"}.
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
#' ggchord(seq_data_example, gene_data = gene_data_example) +
#'   geom_seq() + geom_gene() +
#'   geom_gene_label_repel(label_fit = "shrink", label_min_size = 1.5,
#'     size = 4) + coord_collinear()
geom_gene_label_repel <- function(mapping = NULL, data = NULL,
    label_layout = "radial", label_wrap = NULL, label_fit = "wrap",
    label_max_lines = 2L, label_side = "outside", max_overlaps = Inf,
    label_segment_overlap = "fade", label_segment_overlap_alpha = 0.18,
    label_segment_linetype = "auto", position = "identity",
    show.legend = FALSE, inherit.aes = FALSE, label_min_size = 1.5, ...) {
  dots <- list(...)
  ggchord_reject_label_prefix(dots, "geom_gene_label_repel()")
  args <- c(list(mapping = mapping, data = data,
    gene_label_layout = label_layout, gene_label_wrap = label_wrap,
    gene_label_fit = label_fit, gene_label_max_lines = label_max_lines,
    gene_label_side = label_side, max_overlaps = max_overlaps,
    gene_label_segment_overlap = label_segment_overlap,
    gene_label_segment_overlap_alpha = label_segment_overlap_alpha,
    gene_label_segment_linetype = label_segment_linetype,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    gene_label_min_size = label_min_size), dots)
  if (missing(position)) args$position <- NULL
  do.call(ggchord_label_repel_layer, args)
}

#' @noRd
ggchord_label_repel_layer <- function(mapping = NULL, data = NULL,
                                  gene_label_layout = "radial",
                                  gene_label_wrap = NULL,
                                  gene_label_fit = "wrap",
                                  gene_label_max_lines = 2L,
                                  gene_label_side = "outside",
                                  max_overlaps = Inf,
                                  gene_label_segment_overlap = "fade",
                                  gene_label_segment_overlap_alpha = 0.18,
                                  gene_label_segment_linetype = "auto",
                                  position = "identity",
                                  show.legend = FALSE,
                                  inherit.aes = FALSE,
                                  gene_label_min_size = 1.5,
                                  ...) {
  position_missing <- missing(position)
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  # Inspect the unevaluated call before R's partial argument matching can turn
  # the removed `gene_label_segment` into `gene_label_segment_linetype`.
  raw_argument_names <- names(as.list(sys.call())[-1L])
  dots <- list(...)
  position_object <- ggchord_as_feature_position(position, "geom_gene_label_repel()")
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
      " Repulsion and leader geometry are now selected by label_layout."
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

  if (identical(gene_label_layout, "aligned")) {
    ggchord_stop("geom_gene_label_repel(): removed argument value 'aligned'; ",
      "use label_layout = 'auto' for side columns")
  }
  gene_label_layout <- match.arg(
    gene_label_layout, c("radial", "auto", "arc")
  )
  gene_label_fit <- match.arg(
    gene_label_fit, c("wrap", "none", "ellipsis", "auto", "shrink")
  )
  if (!is.numeric(gene_label_min_size) || length(gene_label_min_size) != 1L ||
      !is.finite(gene_label_min_size) || gene_label_min_size <= 0) ggchord_stop(
    "label_min_size must be one finite positive text size in millimetres")
  if (identical(gene_label_fit, "shrink")) gene_label_layout <- "feature"
  if (!is.numeric(gene_label_max_lines) ||
      length(gene_label_max_lines) != 1L ||
      !is.finite(gene_label_max_lines) || gene_label_max_lines < 1 ||
      gene_label_max_lines != as.integer(gene_label_max_lines)) {
    ggchord_stop("label_max_lines must be one positive integer")
  }
  gene_label_max_lines <- as.integer(gene_label_max_lines)
  gene_label_side <- match.arg(gene_label_side, c("auto", "inside", "outside"))
  gene_label_segment_overlap <- match.arg(
    gene_label_segment_overlap, c("fade", "clip", "show")
  )
  if (!is.numeric(gene_label_segment_overlap_alpha) ||
      length(gene_label_segment_overlap_alpha) != 1L ||
      !is.finite(gene_label_segment_overlap_alpha) ||
      gene_label_segment_overlap_alpha < 0 ||
      gene_label_segment_overlap_alpha > 1) {
    ggchord_stop(
      "label_segment_overlap_alpha must be one finite number in [0, 1]"
    )
  }
  gene_label_segment_linetype <- validate_gene_segment_linetype(
    gene_label_segment_linetype
  )
  empty <- data.frame(
    x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
    label = character(), .component = character(), group = integer(),
    .draw_as_arc = logical(),
    linetype = character(), alpha = numeric(), size = numeric(), angle = numeric(),
    hjust = numeric(), vjust = numeric()
  )
  lyr <- ggplot2::layer(
    data = empty,
    mapping = ggplot2::aes(
      x = x, y = y, xend = xend, yend = yend, label = label,
      group = group, linetype = I(linetype), size = I(size),
      alpha = I(alpha),
      angle = angle, hjust = hjust, vjust = vjust,
      .component = I(.component), .draw_as_arc = I(.draw_as_arc)
    ),
    stat = "identity",
    geom = GeomChordGeneLabelRepel,
    position = position_object,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    check.aes = FALSE,
    check.param = FALSE,
    params = list(
      segment_params = list(),
      text_params = utils::modifyList(dots, list(
        .ggchord_resolved_size = identical(gene_label_fit, "shrink"))),
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
    gene_label_family        = dots$family %||% "",
    gene_label_fontface      = dots$fontface %||% 1,
    gene_label_lineheight    = dots$lineheight %||% 1.2,
    gene_label_wrap          = gene_label_wrap,
    gene_label_fit           = gene_label_fit,
    gene_label_min_size      = gene_label_min_size,
    gene_label_max_lines     = gene_label_max_lines,
    max_overlaps             = max_overlaps,
    gene_label_side          = gene_label_side,
    gene_label_segment_overlap = gene_label_segment_overlap,
    gene_label_segment_overlap_alpha = gene_label_segment_overlap_alpha,
    gene_label_segment_linetype = gene_label_segment_linetype,
    feature_position         = position_object,
    position_supplied        = !position_missing
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping,
    c("accver", "start", "end", "strand", "anno")
  )
  lyr
}

#' Reject renamed public label parameters with an executable migration path
#' @noRd
ggchord_reject_label_prefix <- function(dots, caller) {
  old <- names(dots)[grepl("^gene_label_", names(dots))]
  if (!length(old)) return(invisible(NULL))
  replacements <- sub("^gene_", "", old)
  replacements[old == "gene_label_size"] <- "size"
  replacements[old == "gene_label_max_overlaps"] <- "max_overlaps"
  replacements[old == "gene_label_segment"] <- "label_layout"
  manual <- old %in% c("gene_label_rotation", "gene_label_radial_offset",
    "gene_label_circum_offset", "gene_label_circum_limit", "gene_label_orientation")
  if (grepl("repel", caller)) replacements[manual] <- paste0(
    "geom_gene_label(", replacements[manual], " = ...) for manual placement")
  ggchord_stop(caller, ": parameters removed in v0.16.1: ",
    paste(paste0(old, " -> ", replacements), collapse = "; "))
}
