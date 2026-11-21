# Render once inside final panel containers; never freeze a previous-size image.
ggchord_build_fusion_core <- function(x, on_error = NULL) {
  previous_plot <- last_plot()
  on.exit(set_last_plot(previous_plot),add=TRUE)
  canvas <- c(width=convertWidth(unit(1,"npc"),"mm",TRUE),
    height=convertHeight(unit(1,"npc"),"mm",TRUE))
  origin <- grid::deviceLoc(unit(0,"npc"),unit(0,"npc"),valueOnly=TRUE)
  origin <- c(origin$x,origin$y)*25.4
  collected <- list(boxes=list())
  if (identical(x$fusion_guides,"collect")) {
    collected <- ggchord_fusion_collect_guides(x, partial = !is.null(on_error))
    x <- collected$plot
  }
  for (scope in x$scope_scales) {
    if (any(!scope$layers %in% names(x$ribbons)))
      ggchord_stop("Fusion: scale_for() must select root ribbon layers")
  }
  ribbon_mapped <- lapply(x$ribbons, function(layer) {
    tryCatch(ggchord_fusion_input(layer), error = function(e)
      ggchord_fusion_layer_error(layer, "mapped", e, on_error))
  })
  originals <- x
  x <- ggchord_fusion_subset_ribbons(originals, names(Filter(Negate(is.null), ribbon_mapped)))
  styles <- if (is.null(on_error)) ggchord_fusion_build_styles(x) else
    tryCatch(ggchord_fusion_build_styles(x), error = identity)
  if (inherits(styles, "error")) {
    # Retry only after a joint build fails. Keep the layers that still train
    # together so valid ribbons retain their shared scale and guides.
    selected <- character()
    styles <- ggchord_fusion_build_styles(ggchord_fusion_subset_ribbons(originals, selected))
    for (id in names(x$ribbons)) {
      trial <- ggchord_fusion_subset_ribbons(originals, c(selected, id))
      built <- tryCatch(ggchord_fusion_build_styles(trial), error = identity)
      if (inherits(built, "error")) {
        ggchord_fusion_layer_error(x$ribbons[[id]], "mapped", built, on_error)
      } else {
        selected <- c(selected, id)
        styles <- built
      }
    }
    x <- ggchord_fusion_subset_ribbons(originals, selected)
  }
  ribbon_mapped <- ribbon_mapped[names(x$ribbons)]
  styles <- ggchord_fusion_measure_styles(x, styles)
  styles <- ggchord_fusion_shared_guides(styles,collected$boxes,x,canvas)
  guide_size <- if (styles$side == "none") c(0, 0) else if (styles$side != "multiple")
    styles$sizes[styles$side, ] else apply(styles$sizes, 2, max)
  budget <- ggchord_fusion_guide_budget(canvas, styles$sizes)
  area <- budget$area; shift <- budget$shift
  guide_boxes <- lapply(names(styles$guides), function(side) {
    size <- styles$sizes[side, ]
    xy <- switch(side, right = c(canvas[1]-size[1]/2, shift[2]+area[2]/2),
      left = c(size[1]/2, shift[2]+area[2]/2),
      top = c(shift[1]+area[1]/2, canvas[2]-size[2]/2),
      bottom = c(shift[1]+area[1]/2, size[2]/2))
    data.frame(side = side, xmin = xy[1]-size[1]/2, xmax = xy[1]+size[1]/2,
      ymin = xy[2]-size[2]/2, ymax = xy[2]+size[2]/2, width_mm = size[1], height_mm = size[2])
  })
  guide_boxes <- ggchord_rbind_fill(guide_boxes) %||% data.frame()
  if (min(area) <= 0 || (nrow(guide_boxes) && any(guide_boxes$xmin < 0 | guide_boxes$ymin < 0 |
      guide_boxes$xmax > canvas[1] | guide_boxes$ymax > canvas[2]))) ggchord_stop(
    "Fusion: canvas cannot contain required cross-panel guides; increase width/height or set explicit guide layout")
  if (nrow(guide_boxes) > 1L) for (i in seq_len(nrow(guide_boxes))) for (j in seq_len(i-1L)) {
    a <- guide_boxes[i, ]; b <- guide_boxes[j, ]
    if (min(a$xmax, b$xmax) > max(a$xmin, b$xmin) && min(a$ymax, b$ymax) > max(a$ymin, b$ymin))
      ggchord_stop("Fusion: cross-panel guide boxes overlap under the explicit layout; enlarge the canvas or revise guide layout")
  }
  styles$report$candidate_bounds <- styles$report$bounds
  styles$report$bounds <- if (!nrow(guide_boxes)) data.frame() else data.frame(
    guide = paste0("guide-box-", guide_boxes$side), xmin_mm = guide_boxes$xmin,
    xmax_mm = guide_boxes$xmax, ymin_mm = guide_boxes$ymin, ymax_mm = guide_boxes$ymax,
    overflow = FALSE)
  rectangles <- ggchord_fusion_rectangles(x, area)
  panels <- list(); panel_rows <- list(); instances <- list(); diagnostics <- list()
  backgrounds <- list(); foregrounds <- list(); nested_results <- list(); sources <- list()
  measured_leaves <- list()
  render_events <- list(); connection_grobs <- list(); overlay_grobs <- list(); obstacle_parts <- list()
  for (id in names(rectangles)) {
    panel_path <- ggchord_fusion_path_ids(id)
    spec <- x$panels[[id]]; rect <- rectangles[[id]]
    dimensions <- c(diff(rect[1:2]),diff(rect[3:4]))*area
    center <- c(mean(rect[1:2]),mean(rect[3:4]))*area+shift
    viewport <- grid::viewport(x=center[1],y=center[2],width=dimensions[1],height=dimensions[2],
      default.units="mm",clip="off",name=paste0("fusion-",id))
    panel_error <- function(e) {
      if (is.null(on_error)) stop(e)
      on_error(panel_path, e)
      NULL
    }
    if (inherits(spec$plot, "ggchord_fusion")) {
      child_error <- if (!is.null(on_error)) function(path, error, layer_id = NULL)
        on_error(c(panel_path, path), error, layer_id)
      nested <- tryCatch(local({
        pushViewport(viewport)
        on.exit(popViewport())
        ggchord_build_fusion(spec$plot, child_error)
      }), error = panel_error)
      if (is.null(nested)) next
      offset <- center - dimensions / 2
      child <- ggchord_fusion_nested_export(nested$export, id, offset)
      nested_results[[id]] <- list(export = child)
      obstacle_parts[[id]] <- child$routing_obstacles[child$routing_obstacles$kind != "container", , drop = FALSE]
      for (key in names(nested$panel_registry)) {
        full <- paste0(id, "/", key)
        if (full %in% names(panels) || full %in% names(x$panels)) ggchord_stop("Fusion: nested panel path collision: ", full)
        panel <- nested$panel_registry[[key]]
        panel$context$panel_bounds_mm <- panel$context$panel_bounds_mm + rep(offset, each = 2L)
        panels[[full]] <- panel
      }
      container <- data.frame(panel_id = id, side = "area", coordinate = "fusion", is_container = TRUE,
        container_xmin = offset[1], container_xmax = offset[1] + dimensions[1],
        container_ymin = offset[2], container_ymax = offset[2] + dimensions[2],
        panel_xmin = offset[1], panel_xmax = offset[1] + dimensions[1],
        panel_ymin = offset[2], panel_ymax = offset[2] + dimensions[2],
        measurement = "measured_nested_composition", geometry_metric_status = "nested", converged = TRUE)
      panel_rows[[id]] <- ggchord_rbind_fill(list(container, child$panels))
      instances[[id]] <- child$instances
      diagnostics[[id]] <- list(panel_id = id, canvas_mm = dimensions, nested = TRUE,
        guide_layout = child$metadata$guide_layout, panel_layout = child$metadata$panel_layout,
        plot_decoration = child$metadata$plot_decoration, content_bounds_mm = child$metadata$content_bounds_mm,
        composition_guide_bounds_mm = child$metadata$guide_bounds_mm)
      for (key in names(child$diagnostics)) {
        diagnostic <- child$diagnostics[[key]]
        diagnostic$panel_id <- paste0(id, "/", key)
        if (!is.null(diagnostic$render_context$panel_bounds_mm))
          diagnostic$render_context$panel_bounds_mm <- diagnostic$render_context$panel_bounds_mm + rep(offset, each = 2L)
        diagnostics[[paste0(id, "/", key)]] <- diagnostic
      }
      sources[[id]] <- nested$export$source_data
      for (event in nested$render_events) {
        if (!is.null(event$panel_id)) event$panel_id <- paste0(id, "/", event$panel_id)
        event$background <- grobTree(event$background, vp = viewport)
        event$foreground <- grobTree(event$foreground, vp = viewport)
        if (!is.null(event$bounds)) event$bounds <- event$bounds + rep(offset, each = 2L)
        if (!is.null(event$container_bounds)) event$container_bounds <- event$container_bounds + rep(offset, each = 2L)
        render_events[[length(render_events) + 1L]] <- event
      }
      connection_viewport <- viewport; connection_viewport$clip <- grid::viewport(clip = "inherit")$clip
      connection_grobs <- c(connection_grobs, lapply(nested$connection_grobs, function(g) grobTree(g, vp = connection_viewport)))
      for (j in seq_along(nested$overlay_grobs)) overlay_grobs <- c(overlay_grobs,
        setNames(list(grobTree(nested$overlay_grobs[[j]], vp = viewport)),
          paste0(id, "/", names(nested$overlay_grobs)[j])))
      next
    }
    if (id %in% names(panels)) ggchord_stop("Fusion: nested panel path collision: ", id)
    # Clone build-local ggproto state without serializing quosure/data-function
    # environments: a teaching environment can contain many unrelated plots.
    character_panel <- !is.null(spec$plot$ggchord_alignment)
    genomic_panel <- inherits(spec$plot, "ggchord")
    plot <- if (genomic_panel) ggchord_render_plot_copy(spec$plot) else spec$plot
    plot$layers <- lapply(plot$layers,function(layer) {
      copy <- ggproto(NULL,layer)
      copy$geom <- ggproto(NULL,layer$geom)
      copy$stat <- ggproto(NULL,layer$stat)
      copy$position <- ggproto(NULL,layer$position)
      copy
    })
    table <- tryCatch(local({
      pushViewport(viewport)
      on.exit(popViewport())
      ggchord_with_render_context(list(available_canvas_mm=dimensions), {
        if (!genomic_panel) {
          built <- ggplot_build(plot)
          layout <- if (character_panel) ggchord_alignment_layout(built, plot$ggchord_alignment)
          table <- ggplot_gtable(built)
          context <- ggchord_panel_context(table, built)
          context$geometry_metric_status <- "measured_positions; glyph_fit_at_draw"
          context$converged <- TRUE
          attr(table, "ggchord_render_context") <- context
          attr(table, "ggchord_render_layout") <- layout
          attr(table, "ggchord_render_built") <- built
          table
        } else ggplotGrob(plot)
      })
    }), error = panel_error)
    if (is.null(table)) next
    # Guide bounds are measured while the actual child viewport is active.
    local({
      pushViewport(viewport)
      on.exit(popViewport())
      obstacle_parts[[id]] <<- ggchord_fusion_guide_obstacles(table, id, origin)
    })
    context <- attr(table,"ggchord_render_context")
    if (!identical(context$status, "measured")) {
      panel_error(simpleError(paste0("Fusion panel '", id, "': actual panel measurement unavailable")))
      next
    }
    if (is.null(context$panel_bounds_mm) && length(context$facet_contexts)) {
      b <- do.call(rbind, lapply(context$facet_contexts, `[[`, "panel_bounds_mm"))
      context$panel_bounds_mm <- c(min(b[, 1L]), max(b[, 2L]), min(b[, 3L]), max(b[, 4L]))
    }
    context$device_panel_bounds_mm <- context$panel_bounds_mm
    context$panel_bounds_mm <- context$panel_bounds_mm - rep(origin, each = 2L)
    layout <- attr(table,"ggchord_render_layout")
    built <- attr(table,"ggchord_render_built")
    ggchord_output_capture_build(plot, built, layout)
    panels[[id]] <- list(layout = layout, built = built, context = context,
      seq_data = built$plot$ggchord$data$seq_data %||% spec$plot$ggchord$data$seq_data,
      guides = ggchord_diagnostic_guides(table, dimensions))
    measured_leaves[[id]] <- panels[[id]]
    if (!is.null(layout))
      obstacle_parts[[paste0(id, "/text")]] <- ggchord_fusion_text_obstacles(panels[[id]], id)
    sources[[id]] <- spec$plot$ggchord_alignment %||% spec$plot$ggchord$sources
    backgrounds[[id]] <- grobTree(ggchord_fusion_partition(table,TRUE),vp=viewport)
    foregrounds[[id]] <- grobTree(ggchord_fusion_partition(table,FALSE),vp=viewport)
    bounds <- context$panel_bounds_mm
    render_events[[length(render_events) + 1L]] <- list(panel_id = id, background = backgrounds[[id]],
      foreground = foregrounds[[id]], bounds = bounds,
      container_bounds = c(center[1]-dimensions[1]/2, center[1]+dimensions[1]/2,
        center[2]-dimensions[2]/2, center[2]+dimensions[2]/2))
    panel_rows[[id]] <- data.frame(panel_id=id,side="area",
      coordinate=if(character_panel) "alignment" else if (!genomic_panel) "cartesian" else
        if (isTRUE(layout$dotplot)) "dotplot" else if(isTRUE(layout$sankey_projected)) "sankey" else
        if(isTRUE(layout$radial_projected)) "radial" else if(isTRUE(spec$plot$coordinates$ggchord_circular)) "circular" else
        if(isTRUE(layout$collinear)) "collinear" else "chord",
      container_xmin=center[1]-dimensions[1]/2,container_xmax=center[1]+dimensions[1]/2,
      container_ymin=center[2]-dimensions[2]/2,container_ymax=center[2]+dimensions[2]/2,
      panel_xmin=bounds[1],panel_xmax=bounds[2],panel_ymin=bounds[3],panel_ymax=bounds[4],
      measurement=context$status,geometry_metric_status=context$geometry_metric_status,
      converged=context$converged,stringsAsFactors=FALSE)
    d <- layout$display_instances
    if (is.null(d)) d <- data.frame(instance_id = layout$seqs %||% character(),
      accver = layout$seqs %||% character(), key = layout$seqs %||% character())
    d$panel_id <- rep(id, nrow(d))
    instances[[id]] <- d
    # Report actual gtable guide bounds without claiming arbitrary label collision coverage.
    guide_at <- which(grepl("^guide-box",table$layout$name) &
      !vapply(table$grobs,inherits,logical(1),"zeroGrob"))
    guide_bounds <- lapply(guide_at,function(at) {
      g <- table$grobs[[at]]
      c(width=convertWidth(grobWidth(g),"mm",TRUE),height=convertHeight(grobHeight(g),"mm",TRUE))
    })
    oversized <- any(vapply(guide_bounds,function(size) any(size>dimensions),logical(1)))
    if (oversized) panel_error(simpleError(
      paste0("Fusion panel '", id, "': guides exceed its container")))
    diagnostics[[id]] <- list(panel_id=id,canvas_mm=dimensions,render_context=context,
      adaptive_layout=context$adaptive_layout,guide_layout=context$guide_layout,
      label_fit=ggchord_label_fit_diagnostics(layout),label_visibility=ggchord_label_visibility(layout),
      guide_sizes_mm=guide_bounds,guide_overflow=oversized,collision_status="not globally measured")
  }
  containers <- lapply(render_events, function(e) {
    if (is.null(e$panel_id) || is.null(e$container_bounds)) return(NULL)
    b <- e$container_bounds
    data.frame(panel_id = e$panel_id, kind = "container", xmin = b[1], xmax = b[2],
      ymin = b[3], ymax = b[4], measurement = "resolved_container")
  })
  if (nrow(guide_boxes)) obstacle_parts[["composition-guide"]] <- data.frame(panel_id = "",
    kind = "guide", guide_boxes[c("xmin", "xmax", "ymin", "ymax")], measurement = "composition_guide")
  obstacles <- unique(ggchord_rbind_fill(c(containers, obstacle_parts))) %||% data.frame()
  connection_rows <- list()
  fronts <- list()
  paths <- list()
  ribbons <- list()
  for (i in seq_along(x$ribbons)) {
    layer <- x$ribbons[[i]]
    input <- ribbon_mapped[[layer$id]]
    if (!is.null(layer$capture)) {
      layer$capture$mapped <- input
      layer$capture$statistic <- input
    }
    incomplete <- ggchord_fusion_incomplete(input, layer$required)
    if (any(incomplete) && !layer$na.rm)
      warning("Fusion layer '", layer$id, "': omitted ", sum(incomplete), " incomplete relationship row(s)", call. = FALSE)
    for (row in seq_len(nrow(input))) {
      record <- input[row, , drop = FALSE]
      report <- record
      report$qpanel_path <- list(report$qpanel_path[[1L]])
      report$spanel_path <- list(report$spanel_path[[1L]])
      report$layer_id <- layer$id
      report$composition_id <- ""
      report$source_row <- row
      report$component_type <- "relationship"
      if (!"source_members" %in% names(report)) report$source_members <- list(row)
      if (incomplete[row]) {
        report$retained <- FALSE
        report$reason <- "missing_input"
        connection_rows[[length(connection_rows) + 1L]] <- report
        next
      }
      connection <- tryCatch({
        q <- ggchord_fusion_endpoint(panels, record, "q", layer$id, row)
        s <- ggchord_fusion_endpoint(panels, record, "s", layer$id, row)
        if (q$panel == s$panel) ggchord_stop("Fusion: links require distinct panels")
        if (q$unit != s$unit) ggchord_stop("Fusion: bp/aa links require an explicit source projection")
        ggchord_check_character_origin(q, s, "fusion_ribbon()")
        ggchord_check_interval_unit(record, q$unit, "fusion_ribbon()")
        qf <- ggchord_fusion_front(panels[[q$panel]], q)
        sf <- ggchord_fusion_front(panels[[s$panel]], s)
        path <- ggchord_fusion_path(qf, sf, if (layer$route == "straight") 0 else layer$curvature)
        style <- styles$data[[i]][row, , drop = FALSE]
        normals <- attr(path, "endpoint_normals")
        row_fronts <- lapply(c("q", "s"), function(prefix) {
          endpoint <- if (prefix == "q") q else s
          front <- if (prefix == "q") qf else sf
          front$endpoint <- prefix
          front$panel_id <- endpoint$panel
          front$instance_id <- endpoint$instance
          front$accver <- endpoint$entity
          front$layer_id <- layer$id
          front$composition_id <- ""
          front$source_row <- row
          front$sequence_unit <- endpoint$unit
          front$control_normal_x <- normals[prefix, 1L]
          front$control_normal_y <- normals[prefix, 2L]
          front
        })
        d <- data.frame(x = path[, 1L], y = path[, 2L], fill = style$fill,
          colour = style$colour, alpha = style$alpha, linewidth = style$linewidth,
          linetype = style$linetype, layer_id = layer$id, composition_id = "",
          source_row = row, units = "mm", component_type = "relationship")
        grob <- grid::pathGrob(d$x, d$y,
          default.units = "mm", rule = "evenodd", gp = gpar(
            fill = alpha(style$fill, style$alpha),
            col = alpha(style$colour, style$alpha), lwd = style$linewidth * ggplot2::.pt,
            lty = style$linetype, lineend = layer$fixed$lineend %||% "butt"))
        report$route <- layer$route
        report$qpanel <- q$panel
        report$spanel <- s$panel
        report$qinstance <- q$instance
        report$sinstance <- s$instance
        report$sequence_unit <- q$unit
        report$retained <- TRUE
        report$reason <- NA_character_
        list(fronts = row_fronts, path = d, grob = grob, report = report)
      }, error = function(e) ggchord_fusion_layer_error(layer, "geometry", e, on_error))
      if (is.null(connection)) {
        report$retained <- FALSE
        report$reason <- "geometry_failed"
      } else {
        fronts <- c(fronts, connection$fronts)
        paths[[length(paths) + 1L]] <- connection$path
        ribbons[[length(ribbons) + 1L]] <- connection$grob
        report <- connection$report
      }
      connection_rows[[length(connection_rows) + 1L]] <- report
    }
  }
  zooms <- do.call(c, unname(lapply(x$zooms, function(layer) {
    tryCatch(ggchord_render_zoom(layer, panels), error = function(e) {
      if (is.null(on_error)) stop(e)
      if (!is.null(layer$capture)) layer$capture$geometry_error <- conditionMessage(e)
      on_error(character(), e, layer$id)
      NULL
    })
  })))
  guide_grobs <- lapply(seq_len(nrow(guide_boxes)), function(i) {
    b <- guide_boxes[i, ]
    grobTree(styles$guides[[b$side]], vp = grid::viewport(
      x = mean(c(b$xmin, b$xmax)), y = mean(c(b$ymin, b$ymax)),
      width = b$width_mm, height = b$height_mm, default.units = "mm"))
  })
  guide <- do.call(grobTree, guide_grobs)
  overlaps <- list()
  ids <- names(rectangles)
  for (i in seq_along(ids)) for (j in seq_len(i - 1L)) {
    a <- rectangles[[i]]; b <- rectangles[[j]]
    box <- c(max(a[1], b[1]), min(a[2], b[2]), max(a[3], b[3]), min(a[4], b[4]))
    if (box[1] >= box[2] || box[3] >= box[4]) next
    box <- box * rep(area, each = 2L) + rep(shift, each = 2L)
    overlaps[[length(overlaps) + 1L]] <- data.frame(under = ids[j], over = ids[i],
      xmin = box[1], xmax = box[2], ymin = box[3], ymax = box[4], units = "mm",
      policy = "explicit_allow; later panel paints over earlier panel")
  }
  connection_grobs <- c(ribbons, lapply(zooms, `[[`, "grob"), connection_grobs)
  overlay_grobs <- c(overlay_grobs, do.call(c, unname(lapply(zooms, `[[`, "overlay"))))
  render_events[[length(render_events) + 1L]] <- list(panel_id = NULL,
    background = nullGrob(), foreground = guide, bounds = NULL)
  overlap_mode <- isTRUE(x$allow_overlap) || any(vapply(nested_results,
    function(n) nrow(n$export$overlaps) > 0L, logical(1)))
  grob <- ggchord_fusion_render_events(render_events, connection_grobs, overlay_grobs, canvas, overlap_mode)
  bind <- function(rows) if(length(rows)) ggchord_rbind_fill(rows) else data.frame()
  child_tables <- function(name) lapply(nested_results, function(n) n$export[[name]])
  combine <- function(rows, name) bind(c(rows, child_tables(name)))
  export <- list(
    panels = bind(panel_rows), instances = bind(instances), overlaps = combine(overlaps, "overlaps"),
    connections = combine(connection_rows, "connections"), fronts = combine(fronts, "fronts"),
    paths = combine(paths, "paths"), zoom_regions = combine(lapply(zooms, `[[`, "region"), "zoom_regions"),
    zoom_paths = combine(lapply(zooms, `[[`, "path"), "zoom_paths"), diagnostics = diagnostics,
    source_layouts = lapply(panels, `[[`, "layout"), source_data = sources,
    routing_obstacles = obstacles,
    metadata = list(canvas_mm = canvas, panel_layout = x$fusion_areas,
      coordinate_space = "composition", units = "mm",
      guide_side = styles$side, guide_size_mm = guide_size, guide_sizes_mm = styles$sizes,
      guide_bounds_mm = guide_boxes, guide_layout = styles$report, measurement = "measured",
      allow_overlap = isTRUE(x$allow_overlap), device = names(dev.cur()),
      font_conditions = "Selected device metrics; fallback and embedding have not been inspected"))
  list(grob = grob, export = export, panel_registry = panels, measured_leaves = measured_leaves,
    ribbon_mapped = ribbon_mapped, nested = nested_results,
    render_events = render_events, connection_grobs = connection_grobs, overlay_grobs = overlay_grobs)
}
