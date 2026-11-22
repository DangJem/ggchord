#' Export the source, mapped, statistical or rendered stages of a plot
#' @param plot A ggchord plot, Fusion or registered character plot.
#' @param layers Optional layer IDs in the root, or a table with list-column
#'   panel_path and character layer_id. An empty path selects the root.
#' @param stages Stage names: source, mapped, statistic and geometry.
#' @param width,height Optional target device dimensions, specified together.
#' @param units Device dimensions in inches, millimetres or centimetres.
#' @param device A png, pdf or svg device, a device function, or NULL.
#' @param dpi Raster device resolution.
#' @return Per-layer stage tables, geometry components, stage status, reference
#'   and provenance information, diagnostics and measurement conditions.
#' @details Fusion export records failed panels by path and failed ribbons by
#' path and layer ID. It retains the stages and geometry of other panels and
#' ribbons. Ribbon tips are returned as decoration geometry separately from
#' the source fronts.
#' `measurement$status` describes the device used for measurement. With the
#' geometry stage, `geometry$metrics` reports each panel's measurement status,
#' geometry metric status, convergence and number of refinement passes. A
#' render using a private estimate device can converge on that device while
#' its geometry metrics remain estimates for the requested output. Native
#' character panels report measured cell positions; glyphs are fitted at draw
#' time. `converged` is missing when physical refinement was not evaluated.
#' Consult `stage_status` separately for the availability of layer data.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' export_ggchord_layout(p_basic, stages = c("source", "statistic"))
#' export_ggchord_layout(p_basic, stages = "geometry",
#'   width = 8, height = 5, units = "in", device = "pdf")
export_ggchord_layout <- function(plot, layers = NULL, stages = "geometry",
    width = NULL, height = NULL, units = "in", device = NULL, dpi = 300) {
  caller <- "export_ggchord_layout()"
  ggchord_output_plot(plot, caller)
  choices <- c("source", "mapped", "statistic", "geometry")
  if (!is.character(stages) || !length(stages) || anyNA(stages) ||
      any(!stages %in% choices) || anyDuplicated(stages))
    ggchord_stop(caller, ": stages must be unique names from ", paste(choices, collapse = ", "))
  request <- ggchord_output_request(width, height, units, device, dpi, caller)
  selections <- ggchord_output_selection(layers, caller)
  catalogue <- new.env(parent = emptyenv())
  catalogue$layers <- list()
  catalogue$reference <- list()
  failures <- list()
  record_failure <- function(path, error, layer_id = NULL) {
    failures[[length(failures) + 1L]] <<- list(path = path, layer_id = layer_id,
      message = conditionMessage(error))
  }
  catalogue$on_error <- record_failure
  prospective <- ggchord_output_prepare(plot, character(), stages, catalogue)
  entries <- catalogue$layers
  if (!is.null(selections)) {
    at <- vapply(selections, function(selection) {
      hit <- which(vapply(entries, function(entry) identical(entry$layer_id, selection$layer_id) &&
        identical(entry$panel_path, selection$panel_path), logical(1)))
      if (length(hit) != 1L)
        ggchord_stop(caller, ": unknown layer '", selection$layer_id, "' in the selected path")
      hit
    }, integer(1))
    if (anyDuplicated(at)) ggchord_stop(caller, ": select each layer only once")
    entries <- entries[at]
  }

  need_geometry <- "geometry" %in% stages && !is.null(request$size)
  need_build <- need_geometry || any(stages %in% c("mapped", "statistic"))
  failure <- NULL
  measurement_device <- NULL
  diagnostic <- data.frame()
  results <- NULL
  geometry <- list(layers = list(), panels = data.frame(), fronts = data.frame(),
    paths = data.frame(), text = data.frame(), zoom_regions = data.frame(), zoom_paths = data.frame())
  if ("source" %in% stages)
    for (entry in entries) ggchord_output_capture_source(entry)
  if (need_build) {
    results <- tryCatch({
      if (need_geometry) {
        ggchord_with_output_device(request, function() {
          measurement_device <<- names(dev.cur())
          result <- ggchord_output_build_geometry(prospective, record_failure)
          ggchord_output_attach_geometry(entries, result)
          result
        })
      } else {
        ggchord_output_build_stages(prospective, record_failure)
        NULL
      }
    }, error = function(e) {
      failure <<- conditionMessage(e)
      NULL
    })
    if (!is.null(failure)) diagnostic <- ggchord_output_issue("error", "build_failed", failure)
    if (!is.null(results$fusion)) {
      fusion <- results$fusion$export
      for (field in intersect(names(geometry), names(fusion))) geometry[[field]] <- fusion[[field]]
      if (!is.null(selections)) {
        for (field in c("fronts", "paths", "zoom_regions", "zoom_paths")) {
          d <- geometry[[field]]
          if (!nrow(d)) next
          keep <- vapply(seq_len(nrow(d)), function(row)
            any(vapply(entries, function(entry)
              identical(entry$layer_id, d$layer_id[[row]]) &&
                identical(ggchord_output_composition_key(entry$panel_path), d$composition_id[[row]]),
              logical(1))), logical(1))
          geometry[[field]] <- d[keep, , drop = FALSE]
        }
      }
      geometry$panels <- ggchord_output_panel_paths(geometry$panels)
      geometry$fronts <- ggchord_output_panel_paths(geometry$fronts)
    } else if (!is.null(results$leaf)) {
      geometry$panels <- ggchord_output_panel_table(results$leaf$context)
    }
  }
  if (length(failures)) {
    diagnostic <- ggchord_rbind_fill(c(list(diagnostic), lapply(failures, function(issue)
      ggchord_output_issue("error", "build_failed", issue$message, issue$path,
        issue$layer_id %||% NA_character_))))
  }
  measurement <- ggchord_output_measurement(request, !is.null(measurement_device),
    failure = if (need_geometry && is.null(measurement_device)) failure,
    actual = measurement_device)
  statuses <- lapply(entries, function(entry) {
    local <- Filter(function(issue)
      length(issue$path) <= length(entry$panel_path) &&
        identical(issue$path, head(entry$panel_path, length(issue$path))) &&
        (is.null(issue$layer_id) || identical(issue$layer_id, entry$layer_id)), failures)
    ggchord_output_stage_status(entry, stages, request,
      if (length(local)) local[[1L]]$message else failure)
  })
  report <- ggchord_rbind_fill(statuses) %||% data.frame(
    layer_id = character(), stage = character(), status = character(), reason = character(),
    panel_path = I(list()))
  stage_tables <- setNames(lapply(stages, function(stage) {
    lapply(entries, function(entry) {
      list(panel_path = entry$panel_path, layer_id = entry$layer_id, data = entry$capture[[stage]])
    })
  }), stages)
  if ("geometry" %in% stages) {
    geometry$layers <- stage_tables$geometry
    geometry$metrics <- ggchord_output_geometry_metrics(results$leaves, measurement$status,
      failure, failures)
    text <- lapply(entries, function(entry) {
      value <- entry$capture$text
      if (is.null(value) || !nrow(value)) return(NULL)
      value$layer_id <- rep(entry$layer_id, nrow(value))
      value$panel_path <- rep(list(entry$panel_path), nrow(value))
      value
    })
    geometry$text <- ggchord_rbind_fill(text) %||% data.frame()
    stage_tables$geometry <- geometry
  }
  structure(c(stage_tables, list(stage_status = report, measurement = measurement,
    reference = catalogue$reference, diagnostics = diagnostic,
    provenance = lapply(entries, function(entry) {
      source <- entry$capture$source
      list(panel_path = entry$panel_path, layer_id = entry$layer_id,
        import = attr(source, "ggchord_import"), conversion = attr(source, "ggchord_conversion"),
        lineage = attr(source, "ggchord_report"))
    }))), class = c("ggchord_layout_export", "list"))
}

ggchord_output_plot <- function(plot, caller) {
  if (!inherits(plot, c("ggchord", "ggchord_fusion")) &&
      (!inherits(plot, "ggplot") || is.null(plot$ggchord_alignment)))
    ggchord_stop(caller, ": expected ggchord, Fusion or registered character plot")
}

ggchord_output_request <- function(width, height, units, device, dpi, caller) {
  units <- ggchord_match_arg(units, c("in", "mm", "cm"))
  if (xor(is.null(width), is.null(height)))
    ggchord_stop(caller, ": supply width and height together")
  positive <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  if ((!is.null(width) && (!positive(width) || !positive(height))) || !positive(dpi))
    ggchord_stop(caller, ": width, height and dpi must be positive finite numbers")
  if (!is.null(device) && !is.function(device) &&
      (!is.character(device) || length(device) != 1L || is.na(device) ||
        !device %in% c("png", "pdf", "svg")))
    ggchord_stop(caller, ": use a png, pdf or svg device name or a device function")
  size <- if (!is.null(width)) c(width = width, height = height) /
    switch(units, "in" = 1, "mm" = 25.4, "cm" = 2.54)
  list(size = size, width = width, height = height, units = units, device = device, dpi = dpi)
}

ggchord_with_output_device <- function(request, code) {
  previous <- dev.cur()
  owned <- NULL
  estimated <- is.null(request$device)
  path <- tempfile("ggchord-output-", fileext = if (is.character(request$device))
    paste0(".", request$device) else ".pdf")
  on.exit({
    if (!is.null(owned) && owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
    unlink(path)
  }, add = TRUE)
  width <- request$size[["width"]]
  height <- request$size[["height"]]
  if (estimated) {
    close <- ggchord_measurement_device(width = width, height = height)
    on.exit(close(), add = TRUE)
  } else if (is.function(request$device)) {
    formals <- names(formals(request$device))
    explicit <- list(width = width, height = height, units = "in", res = request$dpi, dpi = request$dpi)
    args <- explicit[intersect(formals, names(explicit))]
    if ("..." %in% formals) for (name in setdiff(c("width", "height"), names(args)))
      args[[name]] <- explicit[[name]]
    file_arg <- intersect(c("filename", "file"), formals)
    if (!length(file_arg) && "..." %in% formals) file_arg <- "filename"
    if (length(file_arg)) args[[file_arg[1L]]] <- path
    do.call(request$device, args)
  } else {
    ggchord_diagnostic_device(path, request$device, width, height, request$dpi)
  }
  if (dev.cur() == previous) ggchord_stop("The requested measurement device did not open")
  if (!estimated) owned <- dev.cur()
  grid::grid.newpage()
  value <- code()
  attr(value, "measurement_device") <- names(dev.cur())
  value
}

ggchord_output_measurement <- function(request, used, failure = NULL, actual = NULL) {
  status <- if (!is.null(failure)) "failed" else if (!used) "unavailable" else
    if (is.null(request$device)) "estimated" else "measured"
  reason <- if (!is.null(failure)) failure else if (is.null(request$size)) "missing_dimensions" else
    if (!used) "geometry_not_requested" else NA_character_
  list(status = status, reason = reason, width = request$width, height = request$height,
    units = request$units, dpi = request$dpi,
    device = actual %||% if (is.character(request$device)) request$device,
    font_conditions = if (status == "estimated") "Private measurement backend and its fonts were used as an estimate" else
      if (status == "measured") "Metrics from the specified backend; font fallback and embedding were not inspected")
}

ggchord_output_selection <- function(layers, caller) {
  if (is.null(layers)) return(NULL)
  if (is.character(layers)) layers <- data.frame(layer_id = layers,
    panel_path = I(rep(list(character()), length(layers))))
  ggchord_require_columns(layers, c("panel_path", "layer_id"), caller)
  if (!is.character(layers$layer_id) || anyNA(layers$layer_id) ||
      any(!nzchar(layers$layer_id)) || !is.list(layers$panel_path))
    ggchord_stop(caller, ": supply nonempty layer IDs and a list of character paths")
  lapply(seq_len(nrow(layers)), function(i) {
    path <- layers$panel_path[[i]]
    if (length(path)) ggchord_fusion_check_path(path, caller)
    list(panel_path = if (!length(path)) character() else unname(path), layer_id = layers$layer_id[[i]])
  })
}

ggchord_output_layer_ids <- function(plot) {
  ids <- lapply(plot$layers, function(layer) layer$ggchord_layer_id %||% layer$ggchord_params$id)
  declared <- unlist(ids, use.names = FALSE)
  if (anyDuplicated(declared)) ggchord_stop("Duplicate declared layer ID")
  vapply(seq_along(ids), function(i) {
    id <- ids[[i]]
    if (!is.null(id)) return(id)
    id <- paste0("ggplot-layer-", i)
    while (id %in% declared) id <- paste0(id, "-native")
    id
  }, character(1))
}

ggchord_output_copy_native <- function(layer, id) {
  copy <- ggproto(NULL, layer, ggchord_layer_id = id)
  if (is.function(layer$data)) copy$data <- ggchord_render_data_once(layer$data)
  copy
}

ggchord_output_native_stat <- function(layer, capture) {
  parent <- layer$stat
  copy <- ggproto(NULL, layer)
  copy$stat <- ggproto(NULL, parent, compute_layer = function(self, data, params, layout) {
    capture$mapped <- data
    result <- ggproto_parent(parent, self)$compute_layer(data, params, layout)
    capture$statistic <- result
    result
  })
  copy
}

ggchord_output_prepare <- function(plot, path, stages, catalogue) {
  record <- function(id, type, capture, source) {
    catalogue$layers <- c(catalogue$layers, list(list(panel_path = path, layer_id = id,
      type = type, capture = capture, source = source)))
  }
  fresh_capture <- function() {
    value <- new.env(parent = emptyenv())
    value$requested <- unique(c("source",
      if (any(stages %in% c("mapped", "statistic", "geometry"))) "mapped",
      if (any(stages %in% c("statistic", "geometry"))) "statistic", stages))
    value
  }
  if (inherits(plot, "ggchord_fusion")) {
    for (key in names(plot$panels))
      plot$panels[[key]]$plot <- ggchord_output_prepare(plot$panels[[key]]$plot,
        c(path, ggchord_fusion_path_ids(key)), stages, catalogue)
    for (layer in plot$ribbons) {
      capture <- fresh_capture()
      record(layer$id, "fusion_ribbon", capture, local({ value <- layer$data; function() value }))
      plot$ribbons[[layer$id]]$capture <- capture
    }
    for (layer in plot$zooms) {
      capture <- fresh_capture()
      source <- layer$region
      source$component_type <- rep("decoration", nrow(source))
      source$from_panel_path <- rep(list(layer$from), nrow(source))
      source$to_panel_path <- rep(list(layer$to), nrow(source))
      record(layer$id, "fusion_zoom", capture, local({ value <- source; function() value }))
      plot$zooms[[layer$id]]$capture <- capture
    }
    return(plot)
  }
  genomic <- inherits(plot, "ggchord")
  if (genomic) {
    plot <- tryCatch(ggchord_render_inputs(plot), error = function(e) {
      if (!length(path) || is.null(catalogue$on_error)) ggchord_rethrow(e)
      catalogue$on_error(path, e)
      plot
    })
    plot <- ggchord_register_layers(ggchord_render_plot_copy(plot))
    reference <- plot$ggchord$data$seq_data
  } else reference <- plot$ggchord_alignment
  if (!is.null(reference)) catalogue$reference <- c(catalogue$reference,
    list(list(panel_path = path, data = reference)))
  ids <- ggchord_output_layer_ids(plot)
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    capture <- fresh_capture()
    capture$track <- identical(layer$ggchord_type, "track")
    capture$anchor <- layer$ggchord_params$anchor_layer
    native <- is.null(layer$ggchord_type)
    source <- if (native) {
      plot_data <- plot$data
      layer <- ggchord_output_copy_native(layer, ids[i])
      local({ lyr <- layer; global <- plot_data; function() lyr$layer_data(global) })
    } else {
      layer <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
      layer$ggchord_params$stage_capture <- capture
      if (!is.null(layer$ggchord_params$source_stat))
        layer$ggchord_params$source_stat_cache <- new.env(parent = emptyenv())
      local({
        lyr <- layer
        fallback <- ggchord_layer_fallback(lyr, plot$ggchord$data)
        anchor <- lyr$ggchord_params$anchor_layer
        function() {
          if (is.null(anchor) || is.data.frame(lyr$ggchord_input_data))
            return(ggchord_layer_source(lyr, fallback))
          # An anchored callback consumes the target's source/statistic stage.
          # Resolve precisely those two layers without building text geometry.
          at <- match(anchor, ids)
          if (is.na(at)) ggchord_stop("Unknown anchor layer '", anchor, "'")
          branch <- ggchord_prepare_facet_sources(plot)
          owner <- branch$layers[[at]]
          if (!is.null(owner$ggchord_params$source_stat)) {
            owner <- reconstruct_layer(owner, owner$data, mapping = owner$mapping)
            owner$ggchord_params$source_stat_facet <- branch$ggchord$source_facet_setup
          }
          branch$layers <- list(owner, lyr)
          branch <- ggchord_bind_domain_labels(branch)
          ggchord_layer_source(branch$layers[[2L]], fallback)
        }
      })
    }
    if (native) {
      if (any(stages %in% c("mapped", "statistic"))) layer <- ggchord_output_native_stat(layer, capture)
    }
    plot$layers[[i]] <- layer
    record(ids[i], if (native) "native" else layer$ggchord_type, capture, source)
  }
  plot
}

ggchord_output_capture_source <- function(entry) {
  entry$capture$source <- tryCatch(entry$source(), error = function(e) {
    entry$capture$source_error <- conditionMessage(e)
    NULL
  })
}

ggchord_output_build_stages <- function(plot, on_error = NULL, path = character()) {
  if (inherits(plot, "ggchord_fusion")) {
    for (layer in plot$ribbons) {
      layer$capture$mapped <- tryCatch({
        ggchord_fusion_check_style_mapping(layer)
        ggchord_fusion_input(layer)
      }, error = function(e) {
        ggchord_fusion_layer_error(layer, "mapped", e, function(ignored, error, id) {
          if (is.null(on_error)) ggchord_rethrow(error)
          on_error(path, error, id)
        })
      })
      layer$capture$statistic <- layer$capture$mapped
    }
    for (id in names(plot$panels)) {
      child <- c(path, ggchord_fusion_path_ids(id))
      tryCatch(ggchord_output_build_stages(plot$panels[[id]]$plot, on_error, child),
        error = function(e) {
          if (is.null(on_error)) ggchord_rethrow(e)
          on_error(child, e)
        })
    }
  } else {
    built <- ggplot_build(plot)
    if (!is.null(plot$ggchord_alignment)) ggchord_alignment_layout(built, plot$ggchord_alignment)
    ggchord_output_capture_build(plot, built, built$plot$ggchord$ref$layout)
  }
}

ggchord_output_capture_build <- function(plot, built, layout) {
  for (layer in plot$layers) {
    capture <- layer$ggchord_params$stage_capture
    if (is.null(capture)) next
    id <- layer$ggchord_layer_id
    if ("mapped" %in% capture$requested && is.null(capture$mapped))
      capture$mapped <- capture$source
    if ("statistic" %in% capture$requested && is.null(capture$statistic)) {
      data <- if (identical(layer$ggchord_type, "ribbon"))
        layout$layer_layouts[[id]]$ribbon_stat_data else layout$track_statistics[[id]]
      capture$statistic <- ggchord_restore_instance_data(data, layout$display_instances) %||%
        capture$mapped %||% capture$source
    }
  }
}

ggchord_output_build_geometry <- function(plot, on_error = NULL) {
  if (inherits(plot, "ggchord_fusion")) {
    fusion <- ggchord_build_fusion(plot, on_error)
    grid.draw(fusion$grob)
    return(list(fusion = fusion, leaves = fusion$panel_registry))
  }
  built <- ggplot_build(plot)
  if (inherits(plot, "ggchord")) {
    table <- ggchord_render_table(built, force = TRUE)
    leaf <- list(built = attr(table, "ggchord_render_built"),
      layout = attr(table, "ggchord_render_layout"), context = attr(table, "ggchord_render_context"),
      seq_data = built$plot$ggchord$data$seq_data)
  } else {
    layout <- ggchord_alignment_layout(built, plot$ggchord_alignment)
    table <- ggplot_gtable(built)
    leaf <- list(built = built, layout = layout, context = ggchord_panel_context(table, built))
  }
  grid.draw(table)
  leaf$guides <- ggchord_diagnostic_guides(table, ggchord_render_dimensions())
  ggchord_output_capture_build(plot, leaf$built, leaf$layout)
  list(leaf = leaf, leaves = setNames(list(leaf), ""))
}

ggchord_output_composition_key <- function(path) {
  ggchord_fusion_path_key(path)
}

ggchord_output_attach_geometry <- function(entries, results) {
  for (entry in entries) {
    capture <- entry$capture
    if (entry$type %in% c("fusion_ribbon", "fusion_zoom")) {
      tables <- if (entry$type == "fusion_ribbon") c("connections", "fronts", "paths") else
        c("zoom_regions", "zoom_paths")
      key <- ggchord_output_composition_key(entry$panel_path)
      capture$geometry <- setNames(lapply(tables, function(table) {
        value <- results$fusion$export[[table]]
        if (!nrow(value)) return(value)
        value[value$composition_id == key & value$layer_id == entry$layer_id, , drop = FALSE]
      }), tables)
      next
    }
    key <- ggchord_fusion_path_key(entry$panel_path)
    leaf <- if (!nzchar(key)) results$leaf else results$leaves[[key]]
    if (is.null(leaf) || !identical(leaf$context$status, "measured")) next
    data <- leaf$layout$layer_geometry[[entry$layer_id]]
    if (isTRUE(leaf$layout$dotplot))
      data <- list(alignment = leaf$layout$alignment[[entry$layer_id]])
    built <- leaf$built
    ids <- vapply(built$plot$layers, function(layer) layer$ggchord_layer_id %||% "", character(1))
    hit <- which(ids == entry$layer_id)
    if (length(hit)) {
      capture$text <- ggchord_rbind_fill(lapply(hit, function(i) {
        d <- built$data[[i]]
        if (!is.data.frame(d) || !nrow(d)) return(NULL)
        styled <- ggchord_diagnostic_text_styles(d, built$plot$layers[[i]])
        ggchord_rbind_fill(lapply(split(seq_len(nrow(styled)), as.integer(styled$PANEL)), function(rows) {
          panel <- as.integer(styled$PANEL[rows[1L]])
          metric <- leaf$context$facet_contexts[[panel]] %||% leaf$context
          value <- ggchord_diagnostic_text(styled[rows, , drop = FALSE], built$layout$coord,
            built$layout$panel_params[[panel]], metric)
          if (nrow(value)) value$PANEL <- rep(panel, nrow(value))
          value
        }))
      }))
    }
    if (is.null(data)) data <- list(vertices = ggchord_rbind_fill(built$data[hit]) %||% data.frame())
    if (identical(entry$type, "ribbon") && length(hit)) {
      tips <- lapply(hit, function(i) {
        params <- built$plot$layers[[i]]$geom_params
        d <- built$data[[i]]
        if (!nrow(d) || identical(params$tip %||% "none", "none")) return(NULL)
        raw <- built$plot$layers[[i]]$data
        if (is.null(d$source_row) && nrow(d) == nrow(raw)) d$source_row <- raw$source_row
        ggchord_rbind_fill(lapply(split(d, interaction(d$PANEL, d$group, drop = TRUE)), function(piece) {
          panel <- as.integer(piece$PANEL[1L])
          context <- leaf$context$facet_contexts[[panel]] %||% leaf$context
          trace <- ggchord_ribbon_trace(piece, params$tip, params$tip_length,
            built$layout$coord, built$layout$panel_params[[panel]], context$panel_mm)
          ggchord_rbind_fill(lapply(names(trace$tips), function(side) {
            xy <- trace$tips[[side]]
            data.frame(x = xy[, 1L], y = xy[, 2L], PANEL = panel,
              group = piece$group[1L], source_row = piece$source_row[1L] %||% NA_integer_,
              endpoint = side, component_type = "decoration")
          }))
        }))
      })
      if (any(lengths(tips) > 0L)) data$tip <- ggchord_rbind_fill(tips)
    }
    capture$geometry <- lapply(data, ggchord_restore_instance_data,
      instances = leaf$layout$display_instances)
  }
}

ggchord_output_panel_paths <- function(data) {
  if (is.data.frame(data) && "panel_id" %in% names(data))
    data$panel_path <- lapply(data$panel_id, ggchord_fusion_path_ids)
  data
}

ggchord_output_panel_table <- function(context) {
  if (!is.null(context$facet_panels_mm)) return(context$facet_panels_mm)
  if (is.null(context$panel_bounds_mm)) return(data.frame())
  as.data.frame(as.list(context$panel_bounds_mm))
}

ggchord_output_geometry_metrics <- function(leaves, device_status, failure, failures) {
  row <- function(path, panel = NA_integer_, panel_status = "unavailable",
      metric_status = "unavailable", converged = NA, passes = NA_integer_) {
    data.frame(panel_path = I(list(ggchord_fusion_path_ids(path))), PANEL = panel,
      panel_measurement = panel_status, geometry_metric_status = metric_status,
      converged = converged, refinement_passes = passes)
  }
  if (!is.null(failure)) return(row("", metric_status = "failed"))
  if (identical(device_status, "unavailable")) return(row(""))
  measured <- lapply(seq_along(leaves), function(at) {
    path <- names(leaves)[at]
    leaf <- leaves[[at]]
    context <- leaf$context
    panels <- as.integer(leaf$built$layout$layout$PANEL)
    passes <- context$refinement_passes %||% 0L
    ggchord_rbind_fill(lapply(seq_along(panels), function(i) {
      info <- context$facet_contexts[[i]] %||% context
      status <- if (!identical(info$status, "measured")) "unavailable" else
        if (identical(device_status, "estimated")) "estimated" else
        if (!is.null(context$refinement_passes)) info$geometry_metric_status else
        if (isTRUE(leaf$layout$character)) "measured_positions" else "measured"
      row(path, panels[i], info$status, status,
        if (passes > 0L && identical(info$status, "measured")) info$converged else NA,
        as.integer(passes))
    }))
  })
  panel_failures <- Filter(function(x) is.null(x$layer_id), failures)
  failed_paths <- unique(vapply(panel_failures, function(x)
    ggchord_fusion_path_key(x$path), character(1)))
  failed <- lapply(setdiff(failed_paths, names(leaves)), function(path)
    row(path, metric_status = "failed"))
  ggchord_rbind_fill(c(measured, failed)) %||% row("")
}

ggchord_output_stage_status <- function(entry, stages, request, build_failure) {
  capture <- entry$capture
  issue <- function(stage) {
    if (stage == "geometry" && is.null(request$size)) return("missing_dimensions")
    if (stage %in% c("statistic", "geometry") && length(capture$statistic_errors))
      return(paste(capture$statistic_errors, collapse = "\n"))
    capture[[paste0(stage, "_error")]] %||% if (is.null(capture[[stage]]) && stage != "source") build_failure
  }
  rows <- lapply(stages, function(stage) {
    problem <- issue(stage)
    if (identical(problem, "missing_dimensions"))
      return(data.frame(stage = stage, status = "unavailable", reason = problem))
    if (!is.null(problem)) return(data.frame(stage = stage, status = "failed", reason = problem))
    if (!is.null(capture[[stage]])) return(data.frame(stage = stage, status = "available", reason = NA_character_))
    data.frame(stage = stage, status = "unavailable", reason = "stage_not_produced")
  })
  value <- do.call(rbind, rows)
  value$panel_path <- rep(list(entry$panel_path), nrow(value))
  value$layer_id <- rep(entry$layer_id, nrow(value))
  value
}

ggchord_output_issue <- function(severity, code, message, panel_path = character(),
    layer_id = NA_character_, component = "plot", suggestion = NA_character_) {
  data.frame(severity = severity, code = code, panel_path = I(list(panel_path)),
    layer_id = layer_id, component = component, message = message, suggestion = suggestion)
}
ggchord_sequence_display <- function(layout) {
  ids <- layout$seqs
  ref <- layout$sequence_reference
  length <- unname(ref$lens[ids])
  tab <- layout$sequence_layout
  span <- if (!isTRUE(layout$collinear)) abs(ref$ends[ids] - ref$starts[ids]) else
    if (all(c("axis_min", "axis_max") %in% names(tab))) tab$axis_max - tab$axis_min else
      tab$xmax - tab$xmin
  data.frame(accver = ids, source_length = length, display_length = unname(span),
    display_scale_factor = unname(span / length), length_mode = layout$length_mode %||% "proportional",
    sequence_unit = layout$sequence_unit %||% "bp",
    display_unit = if (isTRUE(layout$collinear)) "layout" else "radian")
}
