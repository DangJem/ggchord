# geom-feature.R - generic feature layer (v0.9.0)

feature_polygon_geom <- rename_geom_aes(
  ggplot2::GeomPolygon, renames = c(fill = "feature_fill")
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
  required <- c("seq_id", "start", "end", "strand")
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
  } else if ("anno" %in% names(data)) {
    data$anno
  } else {
    ggchord_stop(
      "geom_feature(): data must contain `type` or map `feature_type` in aes()"
    )
  }
  out <- as.data.frame(data, stringsAsFactors = FALSE)
  out$seq_id <- as.character(out$seq_id)
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
#'   \code{seq_id}, \code{start}, \code{end} and \code{strand} may rename
#'   input columns; ordinary visual mappings are evaluated after geometry is
#'   generated.
#' @param data data.frame with \code{seq_id}, \code{start}, \code{end} and
#'   \code{strand}; optional \code{type}, \code{category} and \code{label}.
#' @param feature_shape Fixed feature geometry used when \code{feature_shape}
#'   is not mapped in \code{aes()}: \code{"arrow"}, \code{"block"},
#'   \code{"chevron"}, or \code{"lollipop"}. The default is \code{"arrow"}.
#'   Use \code{aes(feature_shape = type)} together with
#'   \code{scale_feature_shape_manual()} to map categories to geometry.
#' @param feature_width Optional numeric or named vector controlling feature
#'   width; passed to \code{geom_gene(gene_width = ...)}.
#' @param feature_offset Optional numeric or named vector controlling feature
#'   offset; passed to \code{geom_gene(gene_offset = ...)}.
#' @param position Position adjustment passed to [geom_gene()]. Use
#'   [position_feature_stack()] to stack overlaps on radial lanes.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_gene()}.
#'
#' @return A ggplot2 layer.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' features <- data.frame(seq_id = "MT108731.1",
#'                        start = 1000, end = 4000,
#'                        strand = "+", type = "CDS")
#' p <- ggchord(seq_data_example) + geom_seq() +
#'   geom_feature(data = features)
#' p
geom_feature <- function(mapping = NULL, data = NULL,
                         feature_shape = "arrow",
                         feature_width = NULL,
                         feature_offset = NULL,
                         position = "identity",
                         show.legend = TRUE,
                         inherit.aes = FALSE,
                         ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  mapping <- ggchord_normalize_mapping(mapping)
  dots <- list(...)
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
    "seq_id", "start", "end", "strand", "feature_type", "feature_label",
    "feature_shape"
  )
  placeholder <- data.frame(
    seq_id = character(), start = numeric(), end = numeric(),
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
  gene_args <- list(
    mapping = visual_mapping,
    data = placeholder,
    gene_offset = feature_offset,
    gene_width = feature_width,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes
  )
  lyr <- do.call(geom_gene, c(gene_args, dots))
  lyr$ggchord_input_data <- data
  lyr$ggchord_input_mapping <- mapping
  lyr$ggchord_role_aes <- roles
  fixed_shape <- if (shape_mapped) NULL else feature_shape
  lyr$ggchord_input_transform <- function(x) {
    ggchord_feature_data(x, fixed_shape = fixed_shape)
  }
  lyr$ggchord_params$gene_data_override <- NULL
  lyr$ggchord_params$is_feature <- TRUE
  lyr$ggchord_params$gene_color_scheme <- "manual"
  lyr$ggchord_params$feature_shape_mapped <- shape_mapped
  lyr$ggchord_params$feature_shape <- feature_shape
  names(lyr$mapping)[names(lyr$mapping) == "gene_fill"] <- "feature_fill"
  if (!"feature_fill" %in% names(visual_mapping)) {
    lyr$mapping[["feature_fill"]] <- as.name("anno")
  }
  if (shape_mapped && !"feature_fill" %in% names(visual_mapping)) {
    # Both aesthetics describe the same feature type in the common case.
    # Using the same source column lets ggplot2 merge fill and shape into one
    # compact legend instead of printing two redundant "Feature" guides.
    lyr$mapping[["feature_fill"]] <- as.name(".feature_shape_raw")
  }
  if (is.logical(lyr$show.legend) && !is.null(names(lyr$show.legend))) {
    names(lyr$show.legend)[names(lyr$show.legend) == "gene_fill"] <-
      "feature_fill"
      if (shape_mapped && isTRUE(show.legend)) {
      lyr$show.legend <- c(lyr$show.legend, feature_shape = TRUE)
    }
  }
  lyr$geom <- feature_geom
  if (!shape_mapped) lyr$aes_params$feature_shape <- feature_shape
  lyr
}
