# geom-feature.R - generic feature layer (v0.9.0)

feature_polygon_geom <- rename_geom_aes(
  GeomPolygon, renames = c(fill = "feature_fill")
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

#' Draw generic genomic features
#'
#' A general layer for CDS, tRNA, rRNA, repeat, CRISPR, promoter or user-defined
#' features. It supports directional arrows, blocks, notched chevrons and
#' lollipops, all generated against the sequence's real local curve. Feature
#' fill and geometry use independent role-specific scales.
#'
#' @param mapping Optional aesthetic mapping. Role aesthetics such as
#'   \code{seq_id}, \code{start}, \code{end} and \code{strand} may rename
#'   input columns; ordinary visual mappings are evaluated after geometry is
#'   generated.
#' @param data data.frame with \code{seq_id}, \code{start}, \code{end} and
#'   \code{strand}; optional \code{type}, \code{category} and \code{label}.
#' @param type Column name used as the feature type, default \code{"type"}.
#' @param category Optional column name used for colour grouping; defaults to
#'   \code{type}.
#' @param label Optional column name used for annotation text; defaults to
#'   \code{label} when present, otherwise \code{type}.
#' @param feature_shape Fixed feature geometry used when \code{feature_shape}
#'   is not mapped in \code{aes()}: \code{"arrow"}, \code{"block"},
#'   \code{"chevron"}, or \code{"lollipop"}. The default is \code{"arrow"}.
#'   Use \code{aes(feature_shape = type)} together with
#'   \code{scale_feature_shape_manual()} to map categories to geometry.
#' @param feature_colors Optional named color vector by feature value; unnamed
#'   vectors are recycled positionally.
#' @param feature_width Optional numeric or named vector controlling feature
#'   width; passed to \code{geom_gene(gene_width = ...)}.
#' @param feature_offset Optional numeric or named vector controlling feature
#'   offset; passed to \code{geom_gene(gene_offset = ...)}.
#' @param feature_order Optional feature order for the legend.
#' @param show_legend Logical. Show the feature legend, default TRUE.
#' @param legend_position Position of the feature legend: \code{"left"},
#'   \code{"right"}, \code{"top"}, \code{"bottom"} or \code{"inside"}.
#' @param ... Additional arguments passed to \code{geom_gene()}.
#'
#' @return A list of ggplot2 layers
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' features <- data.frame(seq_id = "MT108731.1",
#'                        start = 1000, end = 4000,
#'                        strand = "+", type = "CDS")
#' p <- ggchord(seq_data_example) + geom_seq() + geom_feature(features)
#' p
geom_feature <- function(mapping = NULL, data = NULL,
                         type = "type",
                         category = NULL,
                         label = "label",
                         feature_shape = "arrow",
                         feature_colors = NULL,
                         feature_width = NULL,
                         feature_offset = NULL,
                         feature_order = NULL,
                         show_legend = TRUE,
                         legend_position = "right",
                         ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  legend_position_supplied <- !missing(legend_position)
  if (legend_position_supplied) {
    ggchord_deprecate_once(
      "geom_feature(legend_position)",
      "guides(feature_fill = guide_ggchord_legend(position = ...))"
    )
  }

  # Preserve the pre-v0.9 positional `geom_feature(data)` call while exposing
  # the standard ggplot2 `mapping, data` signature.
  if (is.data.frame(mapping) && is.null(data)) {
    data <- mapping
    mapping <- NULL
  }

  holder <- list(
    ggchord_input_data = data,
    ggchord_input_mapping = mapping,
    ggchord_role_aes = c("seq_id", "start", "end", "strand", "type",
                         "category", "label")
  )
  data <- ggchord_resolve_layer_input(holder)
  if (!is.data.frame(data)) {
    ggchord_stop("geom_feature(): data must be a data.frame")
  }
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
  required <- c("seq_id", "start", "end", "strand")
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0) {
    ggchord_stop("geom_feature(): data must contain columns ",
                 paste(required, collapse = ", "))
  }

  value_col <- category %||% type
  if (!value_col %in% colnames(data)) {
    ggchord_stop("geom_feature(): value column '", value_col, "' not found in data")
  }

  gene_data <- as.data.frame(data, stringsAsFactors = FALSE)
  gene_data$seq_id <- as.character(data$seq_id)
  gene_data$start <- as.numeric(data$start)
  gene_data$end <- as.numeric(data$end)
  gene_data$strand <- as.character(data$strand)
  # `anno` is currently the value consumed by geom_gene() for its fill. Keep
  # the display label separate so a category column cannot be overwritten by
  # an unrelated label column.
  gene_data$anno <- as.character(data[[value_col]])
  gene_data$label <- if (label %in% colnames(data)) {
    as.character(data[[label]])
  } else {
    gene_data$anno
  }
  if (shape_mapped) {
    shape_value <- tryCatch(
      rlang::eval_tidy(mapping[["feature_shape"]], data = data),
      error = function(e) ggchord_stop(
        "Cannot evaluate `feature_shape` in layer mapping: ",
        conditionMessage(e)
      )
    )
    if (length(shape_value) == 1L && nrow(data) != 1L) {
      shape_value <- rep(shape_value, nrow(data))
    }
    if (length(shape_value) != nrow(data) || anyNA(shape_value)) {
      ggchord_stop(
        "Mapped `feature_shape` must return one non-missing value per input row"
      )
    }
    gene_data$.feature_shape_raw <- as.character(shape_value)
  } else {
    gene_data$.feature_shape_raw <- rep(feature_shape, nrow(gene_data))
  }
  # Preserve the original value columns for traceability.
  if (type %in% colnames(data)) gene_data$type <- as.character(data[[type]])
  if (!is.null(category) && category %in% colnames(data)) {
    gene_data$category <- as.character(data[[category]])
  }

  vals <- unique(gene_data$anno)
  if (is.null(feature_colors)) {
    pal <- chord_default_palette(length(vals))
    names(pal) <- vals
  } else if (is.null(names(feature_colors))) {
    if (length(feature_colors) == 1) {
      pal <- setNames(rep(feature_colors, length(vals)), vals)
    } else if (length(feature_colors) == length(vals)) {
      pal <- setNames(as.character(feature_colors), vals)
    } else {
      ggchord_stop("geom_feature(): feature_colors must be length 1 or match the number of unique features")
    }
  } else {
    unknown <- setdiff(names(feature_colors), vals)
    if (length(unknown) > 0) {
      ggchord_stop("geom_feature(): feature_colors contains unknown value(s): ",
                   paste(unknown, collapse = ", "))
    }
    pal <- chord_default_palette(length(vals))
    names(pal) <- vals
    pal[names(feature_colors)] <- as.character(feature_colors)
  }

  visual_mapping <- mapping
  if (!is.null(visual_mapping)) {
    visual_mapping <- visual_mapping[
      setdiff(names(visual_mapping), holder$ggchord_role_aes)
    ]
  }
  if (shape_mapped) {
    visual_mapping[["feature_shape"]] <- as.name(".feature_shape_raw")
  }
  gene_args <- list(
    mapping = visual_mapping,
    data = gene_data,
    gene_offset = feature_offset,
    gene_width = feature_width,
    gene_color_scheme = "manual",
    gene_colors = pal,
    gene_order = feature_order,
    show_legend = show_legend
  )
  layers <- do.call(geom_gene, c(gene_args, list(...)))
  for (lyr in layers) {
    if (legend_position_supplied) {
      lyr$ggchord_params$legend_position <- legend_position
    }
    lyr$ggchord_params$gene_data_override <- gene_data
    lyr$ggchord_params$is_feature <- TRUE
    lyr$ggchord_params$feature_shape_mapped <- shape_mapped
    lyr$ggchord_params$feature_shape <- feature_shape
    names(lyr$mapping)[names(lyr$mapping) == "gene_fill"] <- "feature_fill"
    if (shape_mapped && is.null(category)) {
      # Both aesthetics describe the same feature type in the common case.
      # Using the same source column lets ggplot2 merge fill and shape into one
      # compact legend instead of printing two redundant "Feature" guides.
      lyr$mapping[["feature_fill"]] <- as.name(".feature_shape_raw")
    }
    if (is.logical(lyr$show.legend) && !is.null(names(lyr$show.legend))) {
      names(lyr$show.legend)[names(lyr$show.legend) == "gene_fill"] <-
        "feature_fill"
      if (shape_mapped && isTRUE(show_legend)) {
        lyr$show.legend <- c(lyr$show.legend, feature_shape = TRUE)
      }
    }
    lyr$geom <- feature_geom
    if (!shape_mapped) lyr$aes_params$feature_shape <- feature_shape
    lyr$ggchord_legacy_scales <- NULL
    lyr <- ggchord_add_legacy_scale(
      lyr,
      (!missing(feature_colors) && !is.null(feature_colors)) ||
        (!missing(feature_order) && !is.null(feature_order)),
      "feature_colors/feature_order", "feature_fill",
      "aes(feature_fill = ...) + scale_feature_fill_manual()"
    )
  }
  layers
}
