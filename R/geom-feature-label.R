# Feature-label layers are public role-specific wrappers around the shared gene
# label measurement, collision and leader-routing engine.

#' Label generic genomic features
#'
#' @param mapping,data Standard layer inputs. Map feature text with `label`.
#' @param label_orientation One of `"horizontal"`, `"radial"`, or `"tangent"`.
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
    label_orientation = c("horizontal", "radial", "tangent"),
    label_side = c("outside", "inside", "auto"),
    label_overlap = c("hide", "nudge", "allow"),
    label_wrap = NULL,
    position = "identity",
    show.legend = FALSE, inherit.aes = FALSE, ...) {
  label_orientation <- match.arg(label_orientation)
  label_side <- match.arg(label_side)
  label_overlap <- match.arg(label_overlap)
  lyr <- geom_gene_label(
    mapping = mapping, data = data,
    gene_label_orientation = label_orientation,
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
  lyr
}

#' Automatically arrange generic feature labels
#'
#' @param label_layout One of `"radial"`, `"auto"`, or `"callout"`.
#' @param label_side One of `"inside"`, `"outside"`, or `"auto"`.
#' @param label_wrap,label_fit,label_max_lines Text fitting controls.
#' @param max_overlaps Maximum unresolved overlaps.
#' @inheritParams geom_feature_label
#' @return A composite text and leader-line layer.
#' @export
geom_feature_label_repel <- function(
    mapping = NULL, data = NULL,
    label_layout = c("radial", "auto", "callout"),
    label_wrap = NULL,
    label_fit = c("wrap", "none", "ellipsis", "auto"),
    label_max_lines = 2L,
    label_side = c("outside", "inside", "auto"),
    max_overlaps = Inf,
    position = "identity",
    show.legend = FALSE, inherit.aes = FALSE, ...) {
  label_layout <- match.arg(label_layout)
  label_fit <- match.arg(label_fit)
  label_side <- match.arg(label_side)
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
  lyr$ggchord_input_transform <- ggchord_feature_label_data
  lyr$ggchord_role_aes <- unique(c(
    lyr$ggchord_role_aes, "label", "feature_label", "feature_type"
  ))
  lyr$ggchord_theme_components <- c(
    segment_params = "ggchord.feature.label.segment",
    text_params = "ggchord.feature.label"
  )
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
  out <- ggchord_feature_data(x, "arrow")
  out$anno <- as.character(label)
  out$label <- as.character(label)
  out
}
