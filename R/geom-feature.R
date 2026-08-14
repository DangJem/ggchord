# geom-feature.R - generic feature layer (v0.9.0)

#' Draw generic genomic features
#'
#' A thin, backwards-compatible convenience layer for CDS, tRNA, rRNA, repeat,
#' CRISPR, promoter or user-defined features. It prepares a gene-compatible
#' table from a \code{type} / \code{category} / \code{label} specification and
#' reuses the proven \code{geom_gene()} geometry and scales.
#'
#' @param data data.frame with \code{seq_id}, \code{start}, \code{end} and
#'   \code{strand}; optional \code{type}, \code{category} and \code{label}.
#' @param type Column name used as the feature type, default \code{"type"}.
#' @param category Optional column name used for colour grouping; defaults to
#'   \code{type}.
#' @param label Optional column name used for annotation text; defaults to
#'   \code{label} when present, otherwise \code{type}.
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
#' data(seq_data_example)
#' features <- data.frame(seq_id = "MT108731.1",
#'                        start = 1000, end = 4000,
#'                        strand = "+", type = "CDS")
#' p <- ggchord(seq_data_example) + geom_seq() + geom_feature(features)
#' ggplot2::ggplot_build(p)
geom_feature <- function(data,
                         type = "type",
                         category = NULL,
                         label = "label",
                         feature_colors = NULL,
                         feature_width = NULL,
                         feature_offset = NULL,
                         feature_order = NULL,
                         show_legend = TRUE,
                         legend_position = "right",
                         ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!is.data.frame(data)) {
    ggchord_stop("geom_feature(): data must be a data.frame")
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

  gene_data <- data.frame(
    seq_id = as.character(data$seq_id),
    start = as.numeric(data$start),
    end = as.numeric(data$end),
    strand = as.character(data$strand),
    stringsAsFactors = FALSE
  )
  if (label %in% colnames(data)) {
    gene_data$anno <- as.character(data[[label]])
  } else if (value_col %in% colnames(data)) {
    gene_data$anno <- as.character(data[[value_col]])
  } else {
    gene_data$anno <- "feature"
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

  layers <- geom_gene(
    data = gene_data,
    gene_offset = feature_offset,
    gene_width = feature_width,
    gene_color_scheme = "manual",
    gene_colors = pal,
    gene_order = feature_order,
    show_legend = show_legend,
    legend_position = legend_position,
    ...
  )
  for (lyr in layers) {
    lyr$ggchord_params$gene_data_override <- gene_data
  }
  layers
}
