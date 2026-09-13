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
#' @param position Feature Position used for the label anchor.
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
    position = "identity",
    show.legend = FALSE, inherit.aes = FALSE, ...) {
  dots <- list(...)
  label_orientation <- match.arg(label_orientation)
  label_side <- match.arg(label_side)
  label_overlap <- match.arg(label_overlap)
  lyr <- geom_gene_label(
    mapping = mapping, data = data,
    gene_label_orientation = if (label_orientation == "feature")
      "tangent" else label_orientation,
    gene_label_side = label_side,
    gene_label_overlap = label_overlap,
    gene_label_wrap = label_wrap,
    position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, ...
  )
  lyr$ggchord_input_transform <- ggchord_feature_label_data
  lyr$ggchord_role_aes <- unique(c(
    lyr$ggchord_role_aes, "label", "feature_label", "feature_type"
  ))
  lyr$ggchord_theme_element <- "ggchord.feature.label"
  lyr$ggchord_params$is_feature_label <- TRUE
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
#'   tangent text is tried inside first, then placed adjacent and nudged along
#'   the local tangent. A light leader appears only after material movement.
#'   Choose another mode for unrestricted external callouts.
#' @param label_side One of `"inside"`, `"outside"`, or `"auto"`.
#' @param label_wrap,label_fit,label_max_lines Text fitting controls.
#' @param max_overlaps Maximum unresolved overlaps.
#' @param external Whether the staged feature layout may place a label outside
#'   the circular backbone after inside and adjacent placement fail. When
#'   `FALSE`, unresolved labels are hidden instead of pushing inward forever.
#'   External labels are rendered as rounded callouts using a lightened form
#'   of the feature's resolved fill. When restriction-site labels coexist on a
#'   circular map, both layers are packed in one ordered perimeter layout.
#' @inheritParams geom_feature_label
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
#' @export
geom_feature_label_repel <- function(
    mapping = NULL, data = NULL,
    label_layout = c("feature", "radial", "auto", "callout"),
    label_wrap = NULL,
    label_fit = c("wrap", "none", "ellipsis", "auto"),
    label_max_lines = 2L,
    label_side = c("outside", "inside", "auto"),
    max_overlaps = Inf,
    external = TRUE,
    position = "identity",
    show.legend = FALSE, inherit.aes = FALSE, ...) {
  dots <- list(...)
  label_layout <- match.arg(label_layout)
  label_fit <- match.arg(label_fit)
  label_side <- match.arg(label_side)
  if (!is.logical(external) || length(external) != 1L || is.na(external)) {
    ggchord_stop("geom_feature_label_repel(): external must be TRUE or FALSE")
  }
  if (identical(label_layout, "feature")) {
    lyr <- geom_gene_label_repel(
      mapping = mapping, data = data,
      gene_label_layout = "radial", gene_label_side = "inside",
      gene_label_wrap = label_wrap, gene_label_fit = label_fit,
      gene_label_max_lines = label_max_lines,
      max_overlaps = max_overlaps,
      position = position, show.legend = show.legend,
      inherit.aes = inherit.aes, ...
    )
    # `feature` is an internal staged mode of the shared composite renderer:
    # fit inside first, nudge only nearby labels, and emit leaders only for
    # labels whose final displacement is visually meaningful.
    lyr$ggchord_params$gene_label_layout <- "feature"
    lyr$ggchord_params$is_feature_label <- TRUE
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
  lyr <- geom_gene_label_repel(
    mapping = mapping, data = data,
    gene_label_layout = core_layout,
    gene_label_wrap = label_wrap,
    gene_label_fit = label_fit,
    gene_label_max_lines = label_max_lines,
    gene_label_side = label_side,
    max_overlaps = max_overlaps,
    position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, ...
  )
  lyr$ggchord_params$gene_label_layout <- label_layout
  lyr$ggchord_params$is_feature_label <- TRUE
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
