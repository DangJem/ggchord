# Feature-label layers are public role-specific wrappers around the shared gene
# label measurement, collision and leader-routing engine.

#' Label generic genomic features
#'
#' @param mapping,data Standard layer inputs. Map feature text with `label`.
#' @param label_orientation One of `"feature"`, `"horizontal"`, `"radial"`,
#'   or `"tangent"`. `"feature"` keeps both internal and adjacent labels
#'   tangent to the feature direction; compact labels move beside the feature.
#' @param label_side One of `"inside"`, `"outside"`, or `"auto"`.
#' @param label_overlap Fixed-label collision policy.
#' @param label_wrap Optional character wrapping width.
#' @inheritParams geom_gene_label
#' @param position Feature Position used for the label anchor.
#'   In \code{coord_collinear()}, an omitted Position shares the matching
#'   feature layer's lanes; an explicit Position changes only label anchors,
#'   leaving feature polygons in their own Position.
#'   If both inputs contain unique `feature_id` values, label inputs may be
#'   subsets or reordered; their accession, interval and strand must still
#'   agree. Without IDs, inheritance requires the same interval rows in order.
#'   Multiple matching geometry layers require an explicit label Position or
#'   unambiguous feature IDs.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional text parameters.
#' @return A ggplot2 layer.
#' @export
geom_feature_label <- function(
    mapping = NULL, data = NULL,
    label_orientation = c("feature", "horizontal", "radial", "tangent"),
    label_side = c("outside", "inside", "auto"),
    label_overlap = c("hide", "nudge", "allow"),
    label_wrap = NULL,
    label_rotation = NULL, label_radial_offset = .04,
    label_circum_offset = NULL, label_circum_limit = NULL,
    position = "identity",
    show.legend = FALSE, inherit.aes = FALSE, ...) {
  position_explicit <- !missing(position)
  dots <- list(...)
  ggchord_reject_label_prefix(dots, "geom_feature_label()")
  label_orientation <- match.arg(label_orientation)
  label_side <- match.arg(label_side)
  label_overlap <- match.arg(label_overlap)
  lyr <- ggchord_label_layer(
    mapping = mapping, data = data,
    gene_label_orientation = if (label_orientation == "feature")
      "tangent" else label_orientation,
    gene_label_side = label_side,
    gene_label_overlap = label_overlap,
    gene_label_wrap = label_wrap,
    gene_label_rotation = label_rotation,
    gene_label_radial_offset = label_radial_offset,
    gene_label_circum_offset = label_circum_offset,
    gene_label_circum_limit = label_circum_limit,
    position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, ...
  )
  lyr$ggchord_input_transform <- ggchord_feature_label_data
  lyr$ggchord_role_aes <- unique(c(
    lyr$ggchord_role_aes, "label", "feature_label", "feature_type"
  ))
  lyr$ggchord_theme_element <- "ggchord.feature.label"
  lyr$ggchord_params$is_feature_label <- TRUE
  lyr$ggchord_params$position_explicit <- position_explicit
  lyr$ggchord_params$gene_label_orientation <- label_orientation
  if (!any(c("colour", "color") %in% names(mapping)) &&
      !any(c("colour", "color") %in% names(dots))) {
    lyr$mapping[["colour"]] <- ggplot2::aes(
      colour = I(feature_label_colour)
    )$colour
  }
  lyr
}

#' Automatically arrange generic feature labels
#'
#' @param label_layout One of `"feature"`, `"radial"`, `"auto"`, or
#'   `"callout"`. The default `"feature"` follows plasmid-map convention:
#'   text is tried inside first, then assigned to a fixed-radius internal
#'   track. All non-external text is laid out character by character along its
#'   selected track radius rather than rotated as one straight string. A light
#'   leader appears when one or more intervening tracks separate the feature
#'   and its label.
#'   Choose another mode for unrestricted external callouts.
#' @param label_side One of `"inside"`, `"outside"`, or `"auto"`.
#' @param label_wrap,label_fit,label_max_lines Text fitting controls.
#'   `label_fit = "shrink"` explicitly uses internal-first placement without
#'   truncating the text. Text size is bounded by `label_min_size`; unresolved
#'   labels keep their requested font and use adjacent/external placement or
#'   hiding. Other fitting modes retain their existing behavior.
#' @param label_min_size Minimum text size in millimetres for explicit shrink,
#'   default 1.5. Must not exceed the requested text size. Original text,
#'   requested/final sizes and fitting reason are exported with the label.
#' @param max_overlaps Maximum unresolved overlaps.
#' @param external Whether the staged feature layout may place a label outside
#'   the circular backbone after inside and adjacent placement fail. When
#'   `FALSE`, unresolved labels are hidden instead of pushing inward forever.
#'   External labels are rendered as rounded callouts using a lightened form
#'   of the feature's resolved fill. When restriction-site labels coexist on a
#'   circular map, both retain their own layout semantics but register final
#'   boxes in one coordinate-owned exterior track region so they cannot occupy
#'   the same space.
#' @inheritParams geom_feature_label
#' @inheritParams geom_gene_label_repel
#' @return A composite text and leader-line layer.
#' @examples
#' seq <- data.frame(accver = "circle", length = 1000)
#' features <- data.frame(
#'   accver = "circle", start = c(100, 130), end = c(360, 155),
#'   directionality = c("forward", "nondirectional"),
#'   anno = c("long feature", "short feature")
#' )
#' tracks <- position_feature_stack(base_position = position_plasmid())
#' ggchord(seq, validate = "none") +
#'   geom_seq() +
#'   geom_feature_plasmid(data = features, position = tracks) +
#'   geom_feature_label_repel(data = features, position = tracks,
#'     external = TRUE) +
#'   coord_circular()
#'
#' # Explicitly fit complete text inside a feature with a readable lower bound.
#' features$label <- c("A longer feature annotation", "short feature")
#' ggchord(seq, validate = "none") + geom_seq() +
#'   geom_feature(data = features, width = .12) +
#'   geom_feature_label_repel(data = features, label_fit = "shrink",
#'     label_min_size = 1.5, size = 4) + coord_collinear()
#' @export
geom_feature_label_repel <- function(
    mapping = NULL, data = NULL,
    label_layout = c("feature", "radial", "auto", "callout"),
    label_wrap = NULL,
    label_fit = c("wrap", "none", "ellipsis", "auto", "shrink"),
    label_max_lines = 2L,
    label_side = c("outside", "inside", "auto"),
    max_overlaps = Inf,
    external = TRUE,
    position = "identity",
    show.legend = FALSE, inherit.aes = FALSE,
    label_min_size = 1.5, label_segment_overlap = "fade",
    label_segment_overlap_alpha = .18, label_segment_linetype = "auto", ...) {
  position_explicit <- !missing(position)
  dots <- list(...)
  ggchord_reject_label_prefix(dots, "geom_feature_label_repel()")
  label_layout <- match.arg(label_layout)
  label_fit <- match.arg(label_fit)
  label_side <- match.arg(label_side)
  if (!is.logical(external) || length(external) != 1L || is.na(external)) {
    ggchord_stop("geom_feature_label_repel(): external must be TRUE or FALSE")
  }
  if (identical(label_layout, "feature")) {
    lyr <- ggchord_label_repel_layer(
      mapping = mapping, data = data,
      gene_label_layout = "radial", gene_label_side = "inside",
      gene_label_wrap = label_wrap, gene_label_fit = label_fit,
      gene_label_min_size = label_min_size,
      gene_label_max_lines = label_max_lines,
      max_overlaps = max_overlaps,
      gene_label_segment_overlap = label_segment_overlap,
      gene_label_segment_overlap_alpha = label_segment_overlap_alpha,
      gene_label_segment_linetype = label_segment_linetype,
      position = position, show.legend = show.legend,
      inherit.aes = inherit.aes, ...
    )
    # `feature` is an internal staged mode of the shared composite renderer:
    # fit inside first, nudge only nearby labels, and emit leaders only when a
    # real radial gap separates the curved glyph envelope from its feature.
    lyr$ggchord_params$gene_label_layout <- "feature"
    lyr$ggchord_params$is_feature_label <- TRUE
    lyr$ggchord_params$position_explicit <- position_explicit
    lyr$ggchord_params$label_side_requested <- label_side
    lyr$ggchord_params$feature_label_external <- external
    lyr$ggchord_input_transform <- ggchord_feature_label_data
    lyr$ggchord_role_aes <- unique(c(
      lyr$ggchord_role_aes, "label", "feature_label", "feature_type"
    ))
    lyr$ggchord_theme_components <- c(
      segment_params = "ggchord.feature.label.segment",
      text_params = "ggchord.feature.label"
    )
    if (!any(c("colour", "color") %in% names(mapping)) &&
        !any(c("colour", "color") %in% names(dots))) {
      lyr$mapping[["colour"]] <- ggplot2::aes(
        colour = I(feature_label_colour)
      )$colour
    }
    lyr$mapping[["feature_label_mode"]] <- ggplot2::aes(
      feature_label_mode = I(feature_label_mode)
    )$feature_label_mode
    lyr$mapping[["feature_label_fill"]] <- ggplot2::aes(
      feature_label_fill = I(feature_label_fill)
    )$feature_label_fill
    return(lyr)
  }
  core_layout <- if (identical(label_layout, "callout")) "auto" else label_layout
  lyr <- ggchord_label_repel_layer(
    mapping = mapping, data = data,
    gene_label_layout = core_layout,
    gene_label_wrap = label_wrap,
    gene_label_fit = label_fit,
    gene_label_min_size = label_min_size,
    gene_label_max_lines = label_max_lines,
    gene_label_side = label_side,
    max_overlaps = max_overlaps,
    gene_label_segment_overlap = label_segment_overlap,
    gene_label_segment_overlap_alpha = label_segment_overlap_alpha,
    gene_label_segment_linetype = label_segment_linetype,
    position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, ...
  )
  lyr$ggchord_params$gene_label_layout <- if (label_fit == "shrink") "feature" else label_layout
  lyr$ggchord_params$is_feature_label <- TRUE
  lyr$ggchord_params$position_explicit <- position_explicit
  lyr$ggchord_params$feature_label_external <- external
  lyr$ggchord_input_transform <- ggchord_feature_label_data
  lyr$ggchord_role_aes <- unique(c(
    lyr$ggchord_role_aes, "label", "feature_label", "feature_type"
  ))
  lyr$ggchord_theme_components <- c(
    segment_params = "ggchord.feature.label.segment",
    text_params = "ggchord.feature.label"
  )
  if (!any(c("colour", "color") %in% names(mapping)) &&
      !any(c("colour", "color") %in% names(dots))) {
    lyr$mapping[["colour"]] <- ggplot2::aes(
      colour = I(feature_label_colour)
    )$colour
  }
  lyr
}

# The shared geometry engine stores the feature category in `anno` because it
# drives feature-fill scales. The shared label engine also consumes `anno`, so
# label wrappers deliberately replace it with the resolved feature label after
# the common input normalization has finished.
ggchord_feature_label_data <- function(x) {
  label <- if ("label" %in% names(x)) {
    x$label
  } else if ("feature_label" %in% names(x)) {
    x$feature_label
  } else if ("anno" %in% names(x)) {
    x$anno
  } else if ("feature_type" %in% names(x)) {
    x$feature_type
  } else if ("type" %in% names(x)) {
    x$type
  } else if ("category" %in% names(x)) {
    x$category
  } else rep(NA_character_, nrow(x))
  # A segmented biological feature may need several polygons but still owns
  # one label. Keep the original feature row for label placement.
  out <- ggchord_feature_data(x, "arrow", expand_segments = FALSE)
  out$anno <- as.character(label)
  out$label <- as.character(label)
  out
}
