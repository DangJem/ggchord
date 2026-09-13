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
  extra_params = c(
    "na.rm", "segment_boundaries", "boundary_colour",
    "boundary_linewidth", "boundary_linetype", "adaptive_outline"
  ),
  draw_key = function(data, params, size) {
    key_glyph_feature(data, params, size)
  },
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        segment_boundaries = TRUE,
                        boundary_colour = "#666666",
                        boundary_linewidth = .25,
                        boundary_linetype = "22",
                        adaptive_outline = TRUE) {
    component <- data$.component %||% rep("polygon", nrow(data))
    polygons <- data[is.na(component) | component == "polygon", , drop = FALSE]
    boundaries <- data[component == "boundary", , drop = FALSE]
    points <- data[component == "point", , drop = FALSE]
    grobs <- list()
    if (nrow(polygons)) {
      if (isTRUE(adaptive_outline) && "feature_shape" %in% names(polygons)) {
        small <- polygons$feature_shape %in% c(
          "compact_arrow", "promoter_arrow", "primer_arrow", "marker",
          "chevron", "lollipop"
        )
        polygons$linewidth[small] <- pmin(polygons$linewidth[small], .18)
      }
      if ("feature_fill_explicit" %in% names(polygons)) {
        explicit <- !is.na(polygons$feature_fill_explicit) &
          nzchar(polygons$feature_fill_explicit)
        polygons$feature_fill[explicit] <-
          polygons$feature_fill_explicit[explicit]
        polygons$fill <- polygons$feature_fill
      }
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPolygon$draw_panel(
        polygons, panel_params, coord
      )
    }
    if (isTRUE(segment_boundaries) && nrow(boundaries)) {
      boundaries$colour <- boundary_colour
      boundaries$linewidth <- boundary_linewidth
      boundaries$linetype <- boundary_linetype
      if ("boundary_linetype" %in% names(boundaries)) {
        explicit <- !is.na(boundaries$boundary_linetype) &
          nzchar(as.character(boundaries$boundary_linetype))
        boundaries$linetype[explicit] <-
          as.character(boundaries$boundary_linetype[explicit])
      }
      if ("boundary_draw_linetype" %in% names(boundaries)) {
        draw_style <- !is.na(boundaries$boundary_draw_linetype) &
          nzchar(as.character(boundaries$boundary_draw_linetype))
        boundaries$linetype[draw_style] <-
          as.character(boundaries$boundary_draw_linetype[draw_style])
      }
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        boundaries, panel_params, coord, lineend = "butt",
        linejoin = "round", na.rm = na.rm
      )
    }
    if (nrow(points)) {
      points$linewidth <- pmax(points$linewidth, .55)
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        points, panel_params, coord, lineend = "butt", linejoin = "round",
        na.rm = na.rm
      )
    }
    do.call(grid::grobTree, grobs)
  }
)

#' Normalize common and mapped feature roles once at build time
#' @noRd
ggchord_feature_data <- function(data, fixed_shape = "arrow",
                                 expand_segments = TRUE) {
  direction_column <- intersect(c("directionality", "direction"), names(data))[1L]
  if (!is.na(direction_column)) {
    direction <- tolower(as.character(data[[direction_column]]))
    direction_map <- c(
      forward = "+", reverse = "-", bidirectional = "+/-", none = ".",
      nondirectional = "."
    )
    if (anyNA(direction) || any(!direction %in% names(direction_map))) {
      ggchord_stop(
        "geom_feature(): direction must be forward, reverse, ",
        "bidirectional, none, or nondirectional"
      )
    }
    data$strand <- unname(direction_map[direction])
  }
  if (isTRUE(expand_segments)) data <- ggchord_expand_feature_segments(data)
  required <- c("accver", "start", "end", "strand")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    ggchord_stop(
      "geom_feature(): data is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  feature_value <- if (all(c("feature_color", "anno") %in% names(data))) {
    # Built-in common features carry their authoring colour and use their
    # biological name as the palette key; this preserves distinctions such as
    # AmpR/TcR versus rop even though all three are CDS records.
    data$anno
  } else if ("feature_type" %in% names(data)) {
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
  out$.feature_biological_strand <- out$strand
  out$strand[out$strand %in% c(".", "+/-")] <- "+"
  out$anno <- as.character(feature_value)
  out$label <- if ("feature_label" %in% names(out)) {
    as.character(out$feature_label)
  } else if ("label" %in% names(out)) {
    as.character(out$label)
  } else {
    out$anno
  }
  if (!"feature_class" %in% names(out)) {
    semantic_source <- if ("feature_type" %in% names(out)) {
      out$feature_type
    } else if ("type" %in% names(out)) {
      out$type
    } else if ("category" %in% names(out)) {
      out$category
    } else out$anno
    out$feature_class <- ggchord_feature_classes(semantic_source)
  } else {
    out$feature_class <- ggchord_feature_classes(out$feature_class)
  }
  out$.feature_shape_raw <- if ("feature_shape" %in% names(out)) {
    as.character(out$feature_shape)
  } else {
    rep(fixed_shape, nrow(out))
  }
  nondirectional <- out$.feature_biological_strand == "." &
    out$.feature_shape_raw == "arrow"
  out$.feature_shape_raw[nondirectional] <- "block"
  allowed <- c(
    "arrow", "compact_arrow", "promoter_arrow", "primer_arrow", "marker",
    "block", "chevron", "lollipop"
  )
  if (anyNA(out$.feature_shape_raw) ||
      (!is.null(fixed_shape) && any(!out$.feature_shape_raw %in% allowed))) {
    ggchord_stop(
      "geom_feature(): unknown feature shape"
    )
  }
  out$type <- out$anno
  if ("display_priority" %in% names(out)) {
    if (!is.numeric(out$display_priority) || anyNA(out$display_priority) ||
        any(!is.finite(out$display_priority))) {
      ggchord_stop("geom_feature(): display_priority must be finite numeric values")
    }
  }
  if ("prioritized_display" %in% names(out) &&
      (!is.logical(out$prioritized_display) || anyNA(out$prioritized_display))) {
    ggchord_stop("geom_feature(): prioritized_display must contain TRUE or FALSE")
  }
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
#'   `display_priority` or logical `prioritized_display`. As an alternative to
#'   \code{strand}, \code{directionality} or \code{direction} may contain
#'   \code{"forward"}, \code{"reverse"}, \code{"bidirectional"}, or
#'   \code{"none"}.
#' @param feature_shape Fixed feature geometry used when \code{feature_shape}
#'   is not mapped in \code{aes()}: \code{"arrow"},
#'   \code{"compact_arrow"}, \code{"promoter_arrow"},
#'   \code{"primer_arrow"}, \code{"marker"}, \code{"block"},
#'   \code{"chevron"}, or \code{"lollipop"}. The default is \code{"arrow"}.
#'   Use \code{aes(feature_shape = type)} together with
#'   \code{scale_feature_shape_manual()} to map categories to geometry.
#' @param shape Optional concise alias for fixed `feature_shape`. Supplying
#'   both is an error.
#' @param feature_width Optional numeric/vector/list controlling width in the
#'   shared feature geometry engine. In a circular feature-label composition,
#'   the base width responds to final label size within bounded limits.
#' @param feature_offset Deprecated placement input. Explicit values are
#'   translated to the former strand-separated geometry and emit a warning.
#' @param arrow_head_length,arrow_head_width Arrow-head dimensions in the
#'   sequence-local frame, shared with [geom_gene()]. Values above one for
#'   `arrow_head_width` create a protruding block-arrow head.
#' @param arrow_head_style Arrowhead construction: `"shouldered"` keeps the
#'   requested protruding head, `"flush"` keeps it level with the body, and
#'   `"triangle"` draws a full triangular wedge.
#' @param short_feature Fallback for arrows too short to hold their requested
#'   head: automatic wedge/block selection, a wedge, or a block.
#' @param segment_boundaries Draw explicit internal marks carried by a
#'   segmented feature. Contiguous source-segment joins use dotted dividers;
#'   a non-solid segment-level \code{line_style} may override them, while a
#'   `cleavage_arrows` coordinate adds or merges a biological cut mark. These
#'   marks never create additional arrowheads.
#' @param boundary_colour,boundary_linewidth,boundary_linetype Appearance of
#'   internal segment joins.
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
                         arrow_head_style = c("shouldered", "flush", "triangle"),
                         short_feature = c("auto", "wedge", "block"),
                         segment_boundaries = TRUE,
                         boundary_colour = "#666666",
                         boundary_linewidth = .25,
                         boundary_linetype = "22",
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
  arrow_head_style <- match.arg(arrow_head_style)
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
  if (!is.logical(segment_boundaries) || length(segment_boundaries) != 1L ||
      is.na(segment_boundaries) || !is.character(boundary_colour) ||
      length(boundary_colour) != 1L || !is.numeric(boundary_linewidth) ||
      length(boundary_linewidth) != 1L || !is.finite(boundary_linewidth) ||
      boundary_linewidth < 0) {
    ggchord_stop("geom_feature(): invalid segment-boundary appearance")
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

  allowed_shapes <- c(
    "arrow", "compact_arrow", "promoter_arrow", "primer_arrow", "marker",
    "block", "chevron", "lollipop"
  )
  shape_mapped <- !is.null(mapping) && "feature_shape" %in% names(mapping)
  if (!shape_mapped &&
      (!is.character(feature_shape) || length(feature_shape) != 1L ||
       is.na(feature_shape) || !feature_shape %in% allowed_shapes)) {
    ggchord_stop(
      "geom_feature(): unknown feature_shape"
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
      list(
        segment_boundaries = segment_boundaries,
        boundary_colour = boundary_colour,
        boundary_linewidth = boundary_linewidth,
        boundary_linetype = boundary_linetype,
        adaptive_outline = !("linewidth" %in% names(dots) ||
          (!is.null(mapping) && "linewidth" %in% names(mapping)))
      ),
      dots
    )
  )
  lyr$ggchord_type <- "gene_poly"
  lyr$ggchord_params <- list(
    type = "gene", gene_offset = feature_offset,
    gene_width = feature_width,
    arrow_head_length = as.numeric(arrow_head_length),
    arrow_head_width = as.numeric(arrow_head_width),
    arrow_head_style = arrow_head_style,
    short_feature = short_feature,
    feature_position = position_object,
    legacy_offset = feature_offset,
    feature_role = "feature",
    gene_color_scheme = "manual",
    gene_colors = NULL, gene_order = NULL,
    is_feature = TRUE,
    feature_shape_mapped = shape_mapped,
    feature_shape = feature_shape,
    feature_fill_mapping = visual_mapping[["feature_fill"]],
    feature_fill_fixed = dots$feature_fill %||% NULL
  )
  if (!"feature_fill" %in% names(visual_mapping) && is.data.frame(data) &&
      all(c("anno", "feature_color") %in% names(data))) {
    usable <- !is.na(data$anno) & nzchar(as.character(data$anno)) &
      !is.na(data$feature_color) & nzchar(as.character(data$feature_color))
    source_colours <- stats::setNames(
      as.character(data$feature_color[usable]), as.character(data$anno[usable])
    )
    lyr$ggchord_params$gene_colors <- source_colours[!duplicated(
      names(source_colours), fromLast = TRUE
    )]
  }
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

#' Draw SnapGene-style features on a circular plasmid
#'
#' A concise plasmid-map preset built on [geom_feature()]. It uses the same
#' generic feature engine, while selecting a compact circular track,
#' protruding shouldered block arrows, and dashed internal segment joins.
#'
#' @inheritParams geom_feature
#' @param feature_width Track width; the default matches the compact feature
#'   bands used by desktop plasmid editors.
#' @param position Feature position. `NULL` (the default) automatically stacks
#'   overlapping features on compact lanes inside the plasmid backbone. Reuse
#'   the same explicit [position_feature_stack()] object for a separate label
#'   layer so polygons and labels share their lanes. Optional numeric
#'   `display_priority` (larger first) or logical `prioritized_display` can
#'   override the default longest-first ordering for overlapping features.
#' @return A ggplot2 layer.
#' @export
geom_feature_plasmid <- function(mapping = NULL, data = NULL,
                                 feature_shape = "arrow", shape = NULL,
                                 feature_width = .058,
                                 arrow_head_length = .055,
                                 arrow_head_width = 1.32,
                                 arrow_head_style = "shouldered",
                                 short_feature = "auto",
                                 segment_boundaries = TRUE,
                                 boundary_colour = "#666666",
                                 boundary_linewidth = .25,
                                 boundary_linetype = "22",
                                 position = NULL,
                                 show.legend = FALSE,
                                 inherit.aes = FALSE, ...) {
  if (is.null(position)) {
    position <- position_feature_stack(
      spacing = .10, base_position = position_plasmid()
    )
  }
  geom_feature(
    mapping = mapping, data = data, feature_shape = feature_shape,
    shape = shape, feature_width = feature_width,
    arrow_head_length = arrow_head_length,
    arrow_head_width = arrow_head_width,
    arrow_head_style = arrow_head_style,
    short_feature = short_feature,
    segment_boundaries = segment_boundaries,
    boundary_colour = boundary_colour,
    boundary_linewidth = boundary_linewidth,
    boundary_linetype = boundary_linetype,
    position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, ...
  )
}
