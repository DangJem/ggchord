# layout-export.R - stable public layout export

#' Export computed ggchord geometry
#'
#' Builds a ggchord plot and returns selected computed geometry in a stable,
#' layer-aware structure. Every exported row records its source layer and,
#' where applicable, its source input row. Original mapped columns already
#' attached to the geometry are preserved.
#' This is the recommended interface for consuming layout tables;
#' [get_chord_layout()] remains available for the complete layout list.
#' Both retain their documented fields. See [ggchord_layout_contract] for
#' units, provenance and stability boundaries.
#'
#' @param plot A ggchord plot.
#' @param include Character vector selecting any of \code{"seq"},
#'   \code{"ribbon"}, \code{"gene"}, \code{"feature"}, \code{"labels"},
#'   \code{"alignment"}, \code{"link"}, \code{"restriction"}, \code{"track"}, and \code{"axis"}. Point links and
#'   restriction sites require explicit selection.
#' Dotplots default to `alignment` and `axis`; `axis_layout` contains reversible
#' axis offsets, directions and source origins. Alignment exports retain all four
#' endpoints, layer identity and source row; omissions appear in `connection_report`.
#' @param original_data Logical. Include the resolved per-layer input tables
#'   under \code{original_data}, default \code{FALSE}.
#'
#' @return An object of class \code{ggchord_layout_export}. Selected component
#'   names contain data frames. Circular layouts additionally contain an
#'   \code{annotation_registry} data frame. Its inner rows expose the physical
#'   track, radial width and boundaries; outer rows expose region, sector,
#'   annotation class, dominant and actual band, spill reason, leader-crossing
#'   count, slot, bounding box and leader corridor.
#'   Explicitly inside primer and restriction callouts
#'   register measured inner boxes rather than exterior bands or slots.
#'   The \code{metadata} member
#'   records rotation, aspect ratio, limits, native `sequence_unit` (bp/aa), and
#'   transformation state.
#'   Collinear layouts instead include `sequence_layout`, containing each
#'   sequence's horizontal range, backbone y, orientation, `units_per_position`,
#'   and `sequence_unit`. The compatibility field `units_per_bp` is NA for aa. It also records
#'   reversible x offset, anchor and source origin. `orientation_diagnostics`
#'   records directional support, confidence and unresolved cases.
#'   Position lanes remain local offsets and are not circular physical tracks.
#'   Their `label_diagnostics` records placed or hidden labels after shared
#'   linear collision handling, including source layer, row, reason and leader
#'   collisions. `connection_report` records retained and omitted ribbon input
#'   rows with a selection reason. When converted inputs are present,
#'   `conversion_reports` retains their conversion attributes, including source
#'   class, column renaming, sequence information and object metadata.
#'   Sequence break marks are included under `seq`.
#'   Segmented linear displays retain one `sequence_layout` row per window,
#'   including `row_id`, segment order, break flags and display `gap_bp`.
#'   Source accession, origin/end/full length, local orientation, x offset and
#'   units per bp describe each reversible piece; gaps have no source mapping.
#'   Automatic segment break paths have `.decoration = "segment_break"` in
#'   the sequence geometry. `gap_bp` is display spacing, not omitted length.
#'   Quantitative geometry is selected with `track`; `track_registry` records
#'   each ID/sequence band, source layers, numerical limits, unit and physical
#'   measurement status. `track_diagnostics` records retained/missing/clipped
#'   source rows. Numerical track axes and explicit bp reference bars are
#'   included under `axis`; `axis_selection` records selected sequence axes.
#'   Explicit label shrinking adds `label_fit_diagnostics`, preserving original
#'   text, requested/minimum/final sizes and fallback or hiding reasons.
#'   `render_context` distinguishes layout-only estimates from actual panel
#'   measurements. For a target-size measured report use [diagnose_ggchord()].
#'   `source_registry` and `layer_sources` retain constructor and per-layer
#'   annotation identities without overwriting user columns. Explicit layer
#'   tables/functions have layer-local source scope. `annotation_duplicates`
#'   reports shared constructor rows drawn by multiple geometry layers.
#' @md
#' @export
#'
#' @section New coordinate fields:
#' Radial synteny layouts include `radial_frames` (tangent/normal, angle,
#' local origin and inner radius; theta is in radians) and `radial_endpoints` (full HSP fronts).
#' Metadata records cyclic adjacency, display rotation and closed one-based
#' interval boundaries; these never imply a biological sequence closure.
#' Dotplots include independent `axis_layout` frames and original endpoint
#' rows under `alignment`; axis selection and unknown directions appear in
#' `connection_report`.
#'
#' @section Ribbon drawing:
#' `ribbon_rendering` records per-layer fill/tip mode. Exported polygon fronts
#' retain the full source intervals. An optional subject tip is decorative;
#' its shoulder translates the target endpoint span to 88 percent of the
#' channel. Continuous fills use a device-resolved raster mesh and vector
#' outlines. A bounded fill budget warns when dense layers reduce resolution.
#'
#' @examples
#' library(dplyr)
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq() +
#'   geom_link_ribbon()
#' exported <- export_ggchord_layout(p, include = c("seq", "ribbon"))
#' names(exported)
export_ggchord_layout <- function(
    plot,
    include = c("seq", "ribbon", "gene", "feature", "labels", "axis", "track"),
    original_data = FALSE) {
  default_include <- missing(include)
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!inherits(plot, "ggchord") || is.null(plot$ggchord)) {
    ggchord_stop("export_ggchord_layout(): plot must be a ggchord object")
  }
  allowed <- c("seq", "link", "ribbon", "gene", "feature", "labels",
               "restriction", "axis", "track", "alignment")
  if (!is.character(include) || anyNA(include) ||
      any(!include %in% allowed)) {
    ggchord_stop(
      "export_ggchord_layout(): include must contain only: ",
      paste(allowed, collapse = ", ")
    )
  }
  include <- unique(include)
  if (!is.logical(original_data) || length(original_data) != 1L ||
      is.na(original_data)) {
    ggchord_stop("export_ggchord_layout(): original_data must be TRUE or FALSE")
  }

  if (isTRUE(plot$coordinates$ggchord_dotplot)) {
    if (default_include) include <- c("alignment", "axis")
    return(ggchord_export_dotplot(plot, include, original_data))
  }
  plot <- ggchord_register_layers(plot)
  layout <- get_chord_layout(plot, build = TRUE)
  exported <- stats::setNames(vector("list", length(include)), include)
  for (nm in include) exported[[nm]] <- data.frame()
  kept_inputs <- list()

  classify_component <- function(component, layer) {
    switch(component,
      seq = "seq",
      seq_break = "seq",
      seq_scale = "axis",
      track = "track",
      ribbon = "ribbon",
      link = "link",
      gene_poly = if (isTRUE(layer$ggchord_params$is_feature)) {
        "feature"
      } else {
        "gene"
      },
      seq_region = "feature",
      restriction_site = "restriction",
      gene_text = "labels",
      gene_text_repel = "labels",
      gene_label_segment = "labels",
      gene_label_repel = "labels",
      seq_label = "labels",
      seq_center_label = "labels",
      axis_line = "axis",
      axis_seg = "axis",
      axis_text = "axis",
      axis = "axis",
      NA_character_
    )
  }

  seen_layer_ids <- character(0)
  for (layer in plot$layers) {
    layer_id <- layer$ggchord_layer_id %||% ""
    if (!nzchar(layer_id)) next
    if (layer_id %in% seen_layer_ids) next
    seen_layer_ids <- c(seen_layer_ids, layer_id)
    matching <- Filter(function(x) {
      identical(x$ggchord_layer_id %||% "", layer_id)
    }, plot$layers)
    if (any(vapply(matching, function(x) {
      isTRUE(x$ggchord_params$is_feature)
    }, logical(1)))) {
      layer$ggchord_params$is_feature <- TRUE
    }
    registry <- layout$layer_geometry[[layer_id]]
    if (is.null(registry)) next

    for (component in names(registry)) {
      category <- classify_component(component, layer)
      if (is.na(category) || !category %in% include) next
      geometry <- registry[[component]]
      if (!is.data.frame(geometry) || nrow(geometry) == 0L) next
      geometry <- as.data.frame(geometry, stringsAsFactors = FALSE)
      input <- layout$layer_inputs[[layer_id]][[component]]
      if (!"source_row" %in% names(geometry)) {
        geometry$source_row <- if (!is.null(input) &&
            all(c("accver") %in% names(geometry)) &&
            "accver" %in% names(input)) {
          as.integer(match(
            as.character(geometry$accver), as.character(input$accver)
          ))
        } else {
          rep(NA_integer_, nrow(geometry))
        }
      }
      geometry$layer_id <- rep(layer_id, nrow(geometry))
      geometry$component <- rep(component, nrow(geometry))
      exported[[category]] <- if (nrow(exported[[category]]) == 0L &&
          ncol(exported[[category]]) == 0L) {
        geometry
      } else {
        ggchord_rbind_fill(list(exported[[category]], geometry))
      }

      if (isTRUE(original_data)) {
        if (!is.null(input)) {
          if (is.null(kept_inputs[[layer_id]])) kept_inputs[[layer_id]] <- list()
          kept_inputs[[layer_id]][[component]] <- input
        }
      }
    }
  }

  # The sequence axis is now an automatic decoration rather than a public
  # layer, so export it directly from the computed layout under a stable
  # synthetic layer id.
  if ("axis" %in% include && !any(exported$axis$component %in%
      c("axis", "axis_line", "axis_seg", "axis_text"))) {
    automatic_axis <- ggchord_axis_geometry(layout)
    if (nrow(automatic_axis) > 0L) {
      automatic_axis$source_row <- NA_integer_
      automatic_axis$layer_id <- ".automatic-axis"
      automatic_axis$component <- automatic_axis$.component
      exported$axis <- ggchord_rbind_fill(list(exported$axis, automatic_axis))
    }
  }
  if ("axis" %in% include && is.data.frame(layout$track_axes) && nrow(layout$track_axes)) {
    d <- layout$track_axes
    d$layer_id <- ".track-axes"; d$component <- "track_axis"
    exported$axis <- ggchord_rbind_fill(list(exported$axis, d))
  }
  exported$track_registry <- layout$track_registry %||% data.frame()
  exported$track_diagnostics <- layout$track_diagnostics %||% data.frame()
  exported$track_statistics <- layout$track_statistics %||% list()
  if (isTRUE(original_data)) exported$track_source_data <- layout$track_source_data %||% list()
  exported$source_registry <- ggchord_source_registry(plot$ggchord$sources %||% list())
  exported$layer_sources <- ggchord_annotation_source_maps(layout)
  exported$annotation_duplicates <- ggchord_annotation_duplicate_sources(exported$layer_sources)
  exported$axis_selection <- layout$axis_selection %||% data.frame()
  exported$render_context <- layout$render_context %||% ggchord_estimated_render_context(layout)
  exported$label_fit_diagnostics <- ggchord_label_fit_diagnostics(layout)
  exported$label_visibility <- ggchord_label_visibility(layout)

  coord <- plot$coordinates
  exported$connection_report <- layout$connection_report %||% data.frame()
  ribbon_layers <- Filter(function(layer) identical(layer$ggchord_type,"ribbon"),plot$layers)
  exported$ribbon_rendering <- ggchord_rbind_fill(lapply(ribbon_layers,function(layer) data.frame(
    layer_id=layer$ggchord_layer_id,tip=layer$ggchord_params$tip %||% "none",
    fill=if(ggchord_ribbon_uses_flow(layer,plot)) "endpoint_gradient" else "user_or_value_fill",
    geometry_basis=if(identical(layer$ggchord_params$link_branch %||% "none","none")) "full_source_interval_fronts" else "source_fronts_and_explicit_artificial_branch_junctions",
    shoulder_basis="translated_target_endpoint_span_at_channel_fraction_0.88",
    marker_basis="decorative_subject_front_midpoint_not_a_single_match",stringsAsFactors=FALSE))) %||% data.frame()

  fitted <- if (isTRUE(layout$collinear) && identical(coord$fit, "geometry")) {
    ggchord_geometry_limits(layout)
  } else ggchord_adaptive_limits(layout)
  metadata <- list(
    coordinate_space = "ggchord_layout",
    units = "data",
    sequence_unit = layout$sequence_unit %||% "bp",
    length_mode = layout$length_mode %||% "proportional",
    rotation = layout$rotation %||%
      if (isTRUE(coord$ggchord_coord)) coord$rotation else 0,
    rotation_applied = TRUE,
    ratio = coord$ratio %||% 1,
    coord_transform_applied = FALSE,
    fit = coord$fit %||% "labels",
    xlim = coord$user_xlim %||% fitted$xlim,
    ylim = coord$user_ylim %||% fitted$ylim,
    clip = coord$clip %||% "off"
  )
  if (isTRUE(coord$ggchord_circular)) {
    metadata$coordinate <- "circular"
    metadata$gap <- coord$circular_gap
    metadata$direction <- coord$circular_direction
    exported$annotation_registry <-
      layout$circular_annotation_registry %||%
      ggchord_empty_circular_annotation_registry()
  } else if (isTRUE(coord$ggchord_collinear)) {
    metadata$coordinate <- "collinear"
    if (isTRUE(coord$ggchord_sankey)) metadata$coordinate <- "sankey"
    metadata$genomic_scale <- "linear"
    metadata$normal_direction <- "positive_up"
    metadata$seq_axis <- layout$seq_axis %||% "horizontal"
    if (identical(metadata$seq_axis, "vertical")) metadata$normal_direction <- "positive_left"
    if (isTRUE(coord$ggchord_sankey)) metadata$normal_direction <- "positive_outward"
    if (isTRUE(coord$ggchord_radial)) {
      metadata$coordinate <- "radial_synteny"
      metadata$normal_direction <- "counterclockwise_from_outward_spoke"
      metadata$start_angle <- coord$radial_start_angle
      metadata$clockwise <- coord$radial_clockwise
      metadata$inner_radius <- coord$radial_inner_radius
      metadata$adjacency <- "cyclic"
      metadata$interval_fronts <- "closed_one_based_interval_boundaries"
      exported$radial_frames <- ggchord_restore_instance_data(
        transform(layout$radial_frames, accver = key), layout$display_instances)
      exported$radial_endpoints <- ggchord_restore_instance_data(layout$radial_endpoints, layout$display_instances)
    }
    exported$sequence_layout <- layout$sequence_layout
    if (identical(layout$sequence_unit, "aa")) {
      exported$sequence_layout$units_per_bp <- NA_real_
      if ("gap_bp" %in% names(exported$sequence_layout)) {
        exported$sequence_layout$gap_bp <- NA_real_
      }
    }
    exported$label_diagnostics <- layout$collinear_label_diagnostics %||% data.frame()
    exported$orientation_diagnostics <- layout$orientation_diagnostics %||% data.frame()
  } else {
    metadata$coordinate <- "chord"
  }

  exported$metadata <- metadata
  ref <- layout$sequence_reference
  display <- if (isTRUE(layout$collinear)) {
    if (identical(layout$seq_axis, "vertical")) layout$sequence_layout$axis_max - layout$sequence_layout$axis_min else
      layout$sequence_layout$xmax - layout$sequence_layout$xmin
  } else abs(ref$ends[layout$seqs] - ref$starts[layout$seqs])
  exported$sequence_display <- data.frame(accver = layout$seqs,
    source_length = unname(ref$lens[layout$seqs]), display_length = unname(display),
    display_scale_factor = unname(display / ref$lens[layout$seqs]),
    length_mode = metadata$length_mode, sequence_unit = metadata$sequence_unit,
    display_unit = if (isTRUE(layout$collinear)) "layout" else "radian")
  conversions <- lapply(plot$ggchord$data, attr, which = "ggchord_conversion")
  conversions <- Filter(Negate(is.null), conversions)
  for (id in names(layout$layer_inputs)) {
    values <- lapply(layout$layer_inputs[[id]], attr, which = "ggchord_conversion")
    values <- Filter(Negate(is.null), values)
    if (length(values)) conversions[[id]] <- values
  }
  if (length(conversions)) exported$conversion_reports <- conversions
  if (isTRUE(original_data)) exported$original_data <- kept_inputs
  if (!is.null(layout$display_instances)) {
    instances <- layout$display_instances
    exported$display_instances <- instances[, intersect(c("instance_id", "accver", "declared_order", "side"), names(instances))]
    exported$instance_routing <- layout$instance_routing %||% data.frame()
    exported$instance_omissions <- layout$instance_omissions %||% data.frame()
    for (field in intersect(names(exported), c(include, "sequence_display", "sequence_layout",
        "track_registry", "track_diagnostics", "axis_selection", "connection_report", "label_diagnostics")))
      exported[[field]] <- ggchord_restore_instance_data(exported[[field]], instances)
    for (field in intersect(names(exported), c("original_data", "track_statistics")))
      exported[[field]] <- lapply(exported[[field]], function(x) {
        if (is.data.frame(x)) ggchord_restore_instance_data(x, instances) else
          lapply(x, ggchord_restore_instance_data, instances = instances)
      })
  }
  class(exported) <- c("ggchord_layout_export", "list")
  exported
}
