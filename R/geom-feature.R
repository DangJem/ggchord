# geom-feature.R - generic feature layer (v0.9.0)

GeomChordFeatureBase <- ggplot2::ggproto(
  "GeomChordFeatureBase", ggplot2::GeomPolygon
)
feature_polygon_geom <- rename_geom_aes(
  GeomChordFeatureBase, renames = c(fill = "feature_fill")
)
feature_geom_defaults <- feature_polygon_geom$default_aes
feature_geom_defaults$feature_shape <- "arrow"
feature_geom <- ggplot2::ggproto(
  "GeomChordFeature", feature_polygon_geom,
  default_aes = feature_geom_defaults,
  draw_key = function(data, params, size) {
    key_glyph_feature(data, params, size)
  }
)

#' Normalize common and mapped feature roles once at build time
#' @noRd
ggchord_feature_data <- function(data, fixed_shape = "arrow") {
  required <- c("accver", "start", "end", "strand")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    ggchord_stop(
      "geom_feature(): data is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  feature_value <- if ("feature_type" %in% names(data)) {
    data$feature_type
  } else if ("type" %in% names(data)) {
    data$type
  } else if ("category" %in% names(data)) {
    data$category
  } else if ("anno" %in% names(data)) {
    data$anno
  } else {
    rep(NA_character_, nrow(data))
  }
  out <- as.data.frame(data, stringsAsFactors = FALSE)
  out$accver <- as.character(out$accver)
  out$start <- as.numeric(out$start)
  out$end <- as.numeric(out$end)
  out$strand <- as.character(out$strand)
  out$anno <- as.character(feature_value)
  out$label <- if ("feature_label" %in% names(out)) {
    as.character(out$feature_label)
  } else if ("label" %in% names(out)) {
    as.character(out$label)
  } else {
    out$anno
  }
  out$.feature_shape_raw <- if ("feature_shape" %in% names(out)) {
    as.character(out$feature_shape)
  } else {
    rep(fixed_shape, nrow(out))
  }
  allowed <- c("arrow", "block", "chevron", "lollipop")
  if (anyNA(out$.feature_shape_raw) ||
      (!is.null(fixed_shape) && any(!out$.feature_shape_raw %in% allowed))) {
    ggchord_stop(
      "geom_feature(): feature shapes must be arrow, block, chevron, or lollipop"
    )
  }
  out$type <- out$anno
  out
}

#' Draw generic genomic features
#'
#' A general layer for CDS, tRNA, rRNA, repeat, CRISPR, promoter or user-defined
#' features. It supports directional arrows, blocks, notched chevrons and
#' lollipops, all generated against the sequence's real local curve. Feature
#' fill and geometry use independent role-specific scales.
#'
#' @param mapping Optional aesthetic mapping. Role aesthetics such as
#'   \code{accver}, \code{start}, \code{end} and \code{strand} may rename
#'   input columns; ordinary visual mappings are evaluated after geometry is
#'   generated.
#' @param data data.frame with \code{accver}, \code{start}, \code{end} and
#'   \code{strand}; optional \code{type}, \code{category} and \code{label}.
#' @param feature_shape Fixed feature geometry used when \code{feature_shape}
#'   is not mapped in \code{aes()}: \code{"arrow"}, \code{"block"},
#'   \code{"chevron"}, or \code{"lollipop"}. The default is \code{"arrow"}.
#'   Use \code{aes(feature_shape = type)} together with
#'   \code{scale_feature_shape_manual()} to map categories to geometry.
#' @param shape Optional concise alias for fixed `feature_shape`. Supplying
#'   both is an error.
#' @param feature_width Optional numeric/vector/list controlling width in the
#'   shared feature geometry engine.
#' @param feature_offset Deprecated placement input. Explicit values are
#'   translated to the former strand-separated geometry and emit a warning.
#' @param arrow_head_length,arrow_head_width Arrow-head dimensions in the
#'   sequence-local frame, shared with [geom_gene()].
#' @param short_feature Fallback for arrows too short to hold their requested
#'   head: automatic wedge/block selection, a wedge, or a block.
#' @param position Feature placement. Use `"identity"`, `"strand"`,
#'   `"plasmid"`, [position_strand()], [position_plasmid()], or
#'   [position_feature_stack()].
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional fixed polygon aesthetics.
#'
#' @return A ggplot2 layer.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' features <- data.frame(accver = "MT108731.1",
#'                        start = 1000, end = 4000,
#'                        strand = "+", type = "CDS")
#' p <- ggchord(seq_data_example) + geom_seq() +
#'   geom_feature(data = features)
#' p
geom_feature <- function(mapping = NULL, data = NULL,
                         feature_shape = "arrow",
                         shape = NULL,
                         feature_width = NULL,
                         feature_offset = NULL,
                         arrow_head_length = 0.04,
                         arrow_head_width = 1,
                         short_feature = c("auto", "wedge", "block"),
                         position = "identity",
                         show.legend = TRUE,
                         inherit.aes = FALSE,
                         ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  feature_shape_missing <- missing(feature_shape)
  mapping <- ggchord_normalize_mapping(mapping)
  dots <- list(...)
  if (!is.null(shape)) {
    if (!feature_shape_missing) {
      ggchord_stop("geom_feature(): supply only one of `shape` and `feature_shape`")
    }
    feature_shape <- shape
  }
  short_feature <- match.arg(short_feature)
  position_object <- ggchord_as_feature_position(position, "geom_feature()")
  if (!is.null(feature_offset)) {
    if (!ggchord_is_identity_position(position_object)) {
      ggchord_stop("geom_feature(): `feature_offset` cannot be combined with a non-identity position")
    }
    warning(
      "`feature_offset` placement is deprecated. Use a feature Position instead.",
      call. = FALSE
    )
  }
  if (!is.numeric(arrow_head_length) || length(arrow_head_length) != 1L ||
      !is.finite(arrow_head_length) || arrow_head_length < 0 ||
      !is.numeric(arrow_head_width) || length(arrow_head_width) != 1L ||
      !is.finite(arrow_head_width) || arrow_head_width <= 0) {
    ggchord_stop("geom_feature(): arrow head dimensions must be finite non-negative/positive numbers")
  }
  ggchord_reject_retired(dots, "geom_feature()", c(
    type = "aes(feature_type = ...)",
    category = "aes(feature_fill = ...)",
    label = "aes(feature_label = ...)",
    feature_colors = "scale_feature_fill_manual(values = ...)",
    feature_order = "scale_feature_fill_manual(limits = ...)",
    show_legend = "show.legend",
    legend_position = "guides(feature_fill = guide_ggchord_legend(position = ...))"
  ))

  allowed_shapes <- c("arrow", "block", "chevron", "lollipop")
  shape_mapped <- !is.null(mapping) && "feature_shape" %in% names(mapping)
  if (!shape_mapped &&
      (!is.character(feature_shape) || length(feature_shape) != 1L ||
       is.na(feature_shape) || !feature_shape %in% allowed_shapes)) {
    ggchord_stop(
      "geom_feature(): feature_shape must be 'arrow', 'block', ",
      "'chevron', or 'lollipop'"
    )
  }
  roles <- c(
    "accver", "start", "end", "strand", "feature_type", "feature_label",
    "feature_shape"
  )
  placeholder <- data.frame(
    accver = character(), start = numeric(), end = numeric(),
    strand = character(), anno = character(), label = character(),
    type = character(), .feature_shape_raw = character()
  )

  visual_mapping <- mapping
  if (!is.null(visual_mapping)) {
    visual_mapping <- visual_mapping[
      setdiff(names(visual_mapping), roles)
    ]
  }
  if (shape_mapped) {
    visual_mapping[["feature_shape"]] <- as.name(".feature_shape_raw")
  }
  lyr <- ggplot2::layer(
    data = placeholder,
    mapping = ggplot2::aes(
      x = x, y = y, group = group, feature_fill = anno,
      feature_shape = .feature_shape_raw
    ),
    stat = "identity", geom = feature_geom,
    position = position_object,
    show.legend = if (identical(show.legend, TRUE)) {
      c(feature_fill = TRUE, feature_shape = shape_mapped, colour = FALSE)
    } else show.legend,
    inherit.aes = inherit.aes, check.param = FALSE,
    key_glyph = key_glyph_feature,
    params = c(
      if (!("colour" %in% names(dots))) list(colour = "#353A3E") else list(),
      if (!("linewidth" %in% names(dots))) list(linewidth = 0.25) else list(),
      dots
    )
  )
  lyr$ggchord_type <- "gene_poly"
  lyr$ggchord_params <- list(
    type = "gene", gene_offset = feature_offset,
    gene_width = feature_width,
    arrow_head_length = as.numeric(arrow_head_length),
    arrow_head_width = as.numeric(arrow_head_width),
    short_feature = short_feature,
    feature_position = position_object,
    legacy_offset = feature_offset,
    feature_role = "feature",
    gene_color_scheme = "manual",
    gene_colors = NULL, gene_order = NULL,
    is_feature = TRUE,
    feature_shape_mapped = shape_mapped,
    feature_shape = feature_shape
  )
  lyr <- ggchord_capture_layer_input(lyr, data, mapping, roles)
  fixed_shape <- if (shape_mapped) NULL else feature_shape
  lyr$ggchord_input_transform <- function(x) {
    ggchord_feature_data(x, fixed_shape = fixed_shape)
  }
  lyr$ggchord_params$gene_data_override <- NULL
  if (!"feature_fill" %in% names(visual_mapping)) {
    lyr$mapping[["feature_fill"]] <- as.name("anno")
  }
  if (shape_mapped && !"feature_fill" %in% names(visual_mapping)) {
    # Both aesthetics describe the same feature type in the common case.
    # Using the same source column lets ggplot2 merge fill and shape into one
    # compact legend instead of printing two redundant "Feature" guides.
    lyr$mapping[["feature_fill"]] <- as.name(".feature_shape_raw")
  }
  if (shape_mapped) lyr$mapping[["feature_shape"]] <- as.name(".feature_shape_raw")
  if (!shape_mapped) lyr$aes_params$feature_shape <- feature_shape
  lyr$ggchord_obstacle_provider <- ggchord_entity_obstacles
  lyr
}
