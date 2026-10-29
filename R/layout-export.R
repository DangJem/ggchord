# layout-export.R - stable public layout export

#' Export computed ggchord geometry
#'
#' Builds a ggchord plot and returns selected computed geometry in a stable,
#' layer-aware structure. Every exported row records its source layer and,
#' where applicable, its source input row. Original mapped columns already
#' attached to the geometry are preserved.
#'
#' @param plot A ggchord plot.
#' @param include Character vector selecting any of \code{"seq"},
#'   \code{"ribbon"}, \code{"gene"}, \code{"feature"}, \code{"labels"},
#'   \code{"link"}, \code{"restriction"}, \code{"track"}, and \code{"axis"}. Point links and
#'   restriction sites require explicit selection.
#' @param original_data Logical. Include the resolved per-layer input tables
#'   under \code{original_data}, default \code{FALSE}.
#'
#' @return An object of class \code{ggchord_layout_export}. Selected component
#'   names contain data frames. Circular layouts additionally contain an
#'   \code{annotation_registry} data frame. Its inner rows expose the physical
#'   track, radial width and boundaries; outer rows expose region, sector,
#'   annotation class, dominant and actual band, spill reason, leader-crossing
#'   count, slot, bounding box and leader corridor. The \code{metadata} member
#'   records rotation, aspect ratio, limits, coordinate units and
#'   transformation state.
#'   Collinear layouts instead include `sequence_layout`, containing each
#'   sequence's horizontal range, backbone y, orientation, units per bp,
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
#' @md
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq() + geom_link_ribbon()
#' exported <- export_ggchord_layout(p, include = c("seq", "ribbon"))
#' names(exported)
export_ggchord_layout <- function(
    plot,
    include = c("seq", "ribbon", "gene", "feature", "labels", "axis", "track"),
    original_data = FALSE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!inherits(plot, "ggchord") || is.null(plot$ggchord)) {
    ggchord_stop("export_ggchord_layout(): plot must be a ggchord object")
  }
  allowed <- c("seq", "link", "ribbon", "gene", "feature", "labels",
               "restriction", "axis", "track")
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
  exported$axis_selection <- layout$axis_selection %||% data.frame()
  exported$render_context <- layout$render_context %||% ggchord_estimated_render_context(layout)
  fits <- lapply(names(layout$layer_layouts), function(id) {
    # Dependency sub-layouts may carry compatibility label diagnostics even
    # when their owner only draws a backbone or feature. Export only the
    # actual text owner's diagnostics, including owners whose text was hidden.
    if (!any(names(layout$layer_geometry[[id]]) %in%
        c("gene_text", "gene_text_repel", "gene_label_repel"))) return(NULL)
    d <- layout$layer_layouts[[id]]$label_fit_diagnostics
    if (!is.data.frame(d) || !nrow(d)) return(NULL)
    resolved <- lapply(layout$layer_geometry[[id]], function(g) {
      if (!is.data.frame(g) || !nrow(g) || !"label_size_requested" %in% names(g)) return(NULL)
      if (".component" %in% names(g)) g <- g[g$.component %in% c("text", "arc_text"), , drop = FALSE]
      if (!nrow(g)) return(NULL)
      unique(g[intersect(c("source_row", "label_size_final", "label_fit_reason"), names(g))])
    })
    resolved <- ggchord_rbind_fill(resolved)
    if (is.data.frame(resolved) && nrow(resolved)) {
      at <- match(d$source_row, resolved$source_row)
      for (field in intersect(c("label_size_final", "label_fit_reason"), names(resolved)))
        d[[field]][!is.na(at)] <- resolved[[field]][at[!is.na(at)]]
      d$label_fit_reason[is.na(at)] <- "hidden_after_layout"
    } else d$label_fit_reason <- "hidden_after_layout"
    d$layer_id <- id
    d
  })
  exported$label_fit_diagnostics <- ggchord_rbind_fill(fits) %||% data.frame()

  coord <- plot$coordinates
  exported$connection_report <- layout$connection_report %||% data.frame()
  fitted <- if (isTRUE(layout$collinear) && identical(coord$fit, "geometry")) {
    ggchord_geometry_limits(layout)
  } else ggchord_adaptive_limits(layout)
  metadata <- list(
    coordinate_space = "ggchord_layout",
    units = "data",
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
    metadata$genomic_scale <- "linear"
    metadata$normal_direction <- "positive_up"
    exported$sequence_layout <- layout$sequence_layout
    exported$label_diagnostics <- layout$collinear_label_diagnostics %||% data.frame()
    exported$orientation_diagnostics <- layout$orientation_diagnostics %||% data.frame()
  } else {
    metadata$coordinate <- "chord"
  }

  exported$metadata <- metadata
  conversions <- lapply(plot$ggchord$data, attr, which = "ggchord_conversion")
  conversions <- Filter(Negate(is.null), conversions)
  for (id in names(layout$layer_inputs)) {
    values <- lapply(layout$layer_inputs[[id]], attr, which = "ggchord_conversion")
    values <- Filter(Negate(is.null), values)
    if (length(values)) conversions[[id]] <- values
  }
  if (length(conversions)) exported$conversion_reports <- conversions
  if (isTRUE(original_data)) exported$original_data <- kept_inputs
  class(exported) <- c("ggchord_layout_export", "list")
  exported
}
