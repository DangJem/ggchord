# Render once inside final panel containers; never freeze a previous-size image.
ggchord_build_fusion <- function(x) {
  previous_plot <- ggplot2::last_plot()
  on.exit(ggplot2::set_last_plot(previous_plot),add=TRUE)
  canvas <- c(width=grid::convertWidth(grid::unit(1,"npc"),"mm",TRUE),
    height=grid::convertHeight(grid::unit(1,"npc"),"mm",TRUE))
  origin <- grid::deviceLoc(grid::unit(0,"npc"),grid::unit(0,"npc"),valueOnly=TRUE)
  origin <- c(origin$x,origin$y)*25.4
  styles <- ggchord_fusion_styles(x)
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
  if (min(area)<15 || (nrow(guide_boxes) && any(guide_boxes$xmin < 0 | guide_boxes$ymin < 0 |
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
  render_events <- list(); connection_grobs <- list(); overlay_grobs <- list(); obstacle_parts <- list()
  for (id in names(x$panels)) {
    spec <- x$panels[[id]]; rect <- rectangles[[id]]
    dimensions <- c(diff(rect[1:2]),diff(rect[3:4]))*area
    center <- c(mean(rect[1:2]),mean(rect[3:4]))*area+shift
    viewport <- grid::viewport(x=center[1],y=center[2],width=dimensions[1],height=dimensions[2],
      default.units="mm",clip="off",name=paste0("fusion-",id))
    if (inherits(spec$plot, "ggchord_fusion")) {
      nested <- local({
        grid::pushViewport(viewport)
        on.exit(grid::popViewport())
        ggchord_build_fusion(spec$plot)
      })
      offset <- center - dimensions / 2
      child <- ggchord_fusion_nested_export(nested$export, id, offset)
      nested_results[[id]] <- list(export = child, permissions = nested$source_permissions)
      obstacle_parts[[id]] <- child$routing_obstacles[child$routing_obstacles$kind != "container", , drop = FALSE]
      for (key in names(nested$panel_registry)) {
        full <- paste0(id, "/", key)
        if (full %in% names(panels) || full %in% names(x$panels)) ggchord_stop("Fusion: nested panel path collision: ", full)
        panel <- nested$panel_registry[[key]]
        panel$context$panel_bounds_mm <- panel$context$panel_bounds_mm + rep(offset, each = 2L)
        panels[[full]] <- panel
      }
      container <- data.frame(panel_id = id, side = spec$side, coordinate = "fusion", is_container = TRUE,
        container_xmin = offset[1], container_xmax = offset[1] + dimensions[1],
        container_ymin = offset[2], container_ymax = offset[2] + dimensions[2],
        panel_xmin = offset[1], panel_xmax = offset[1] + dimensions[1],
        panel_ymin = offset[2], panel_ymax = offset[2] + dimensions[2],
        measurement = "measured_nested_composition", geometry_metric_status = "nested", converged = TRUE)
      panel_rows[[id]] <- ggchord_rbind_fill(list(container, child$panels))
      instances[[id]] <- child$instances
      diagnostics[[id]] <- list(panel_id = id, canvas_mm = dimensions, nested = TRUE,
        guide_layout = child$metadata$guide_layout, panel_layout = child$metadata$panel_layout,
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
        event$background <- grid::grobTree(event$background, vp = viewport)
        event$foreground <- grid::grobTree(event$foreground, vp = viewport)
        if (!is.null(event$bounds)) event$bounds <- event$bounds + rep(offset, each = 2L)
        if (!is.null(event$container_bounds)) event$container_bounds <- event$container_bounds + rep(offset, each = 2L)
        render_events[[length(render_events) + 1L]] <- event
      }
      connection_viewport <- viewport; connection_viewport$clip <- NA
      connection_grobs <- c(connection_grobs, lapply(nested$connection_grobs, function(g) grid::grobTree(g, vp = connection_viewport)))
      for (j in seq_along(nested$overlay_grobs)) overlay_grobs <- c(overlay_grobs,
        stats::setNames(list(grid::grobTree(nested$overlay_grobs[[j]], vp = viewport)),
          paste0(id, "/", names(nested$overlay_grobs)[j])))
      next
    }
    if (id %in% names(panels)) ggchord_stop("Fusion: nested panel path collision: ", id)
    # Clone build-local ggproto state without serializing quosure/data-function
    # environments: a teaching environment can contain many unrelated plots.
    character_panel <- !is.null(spec$plot$ggchord_alignment)
    plot <- if (character_panel) spec$plot else ggchord_render_plot_copy(spec$plot)
    plot$layers <- lapply(plot$layers,function(layer) {
      copy <- ggplot2::ggproto(NULL,layer)
      copy$geom <- ggplot2::ggproto(NULL,layer$geom)
      copy$stat <- ggplot2::ggproto(NULL,layer$stat)
      copy$position <- ggplot2::ggproto(NULL,layer$position)
      copy
    })
    table <- local({
      grid::pushViewport(viewport)
      on.exit(grid::popViewport())
      ggchord_with_render_context(list(available_canvas_mm=dimensions), {
        if (character_panel) {
          built <- ggplot2::ggplot_build(plot)
          layout <- ggchord_alignment_layout(built, plot$ggchord_alignment)
          table <- ggplot2::ggplot_gtable(built)
          context <- ggchord_panel_context(table, built)
          context$geometry_metric_status <- "measured_positions; glyph_fit_at_draw"
          context$converged <- TRUE
          attr(table, "ggchord_render_context") <- context
          attr(table, "ggchord_render_layout") <- layout
          attr(table, "ggchord_render_built") <- built
          table
        } else ggplot2::ggplotGrob(plot)
      })
    })
    # Guide bounds are measured while the actual child viewport is active.
    local({
      grid::pushViewport(viewport)
      on.exit(grid::popViewport())
      obstacle_parts[[id]] <<- ggchord_fusion_guide_obstacles(table, id, origin)
    })
    context <- attr(table,"ggchord_render_context")
    if (!identical(context$status,"measured")) ggchord_stop("Fusion panel '",id,"': actual panel measurement unavailable")
    context$device_panel_bounds_mm <- context$panel_bounds_mm
    context$panel_bounds_mm <- context$panel_bounds_mm-rep(origin,each=2L)
    layout <- attr(table,"ggchord_render_layout")
    built <- attr(table,"ggchord_render_built")
    panels[[id]] <- list(layout=layout,built=built,context=context,
      seq_data=spec$plot$ggchord$data$seq_data)
    obstacle_parts[[paste0(id, "/text")]] <- ggchord_fusion_text_obstacles(panels[[id]], id)
    sources[[id]] <- spec$plot$ggchord_alignment %||% spec$plot$ggchord$sources
    backgrounds[[id]] <- grid::grobTree(ggchord_fusion_partition(table,TRUE),vp=viewport)
    foregrounds[[id]] <- grid::grobTree(ggchord_fusion_partition(table,FALSE),vp=viewport)
    bounds <- context$panel_bounds_mm
    render_events[[length(render_events) + 1L]] <- list(panel_id = id, background = backgrounds[[id]],
      foreground = foregrounds[[id]], bounds = bounds,
      container_bounds = c(center[1]-dimensions[1]/2, center[1]+dimensions[1]/2,
        center[2]-dimensions[2]/2, center[2]+dimensions[2]/2))
    panel_rows[[id]] <- data.frame(panel_id=id,side=spec$side,
      coordinate=if(character_panel) "alignment" else if(isTRUE(layout$sankey_projected)) "sankey" else
        if(isTRUE(layout$radial_projected)) "radial" else if(isTRUE(spec$plot$coordinates$ggchord_circular)) "circular" else
        if(isTRUE(layout$collinear)) "collinear" else "chord",
      container_xmin=center[1]-dimensions[1]/2,container_xmax=center[1]+dimensions[1]/2,
      container_ymin=center[2]-dimensions[2]/2,container_ymax=center[2]+dimensions[2]/2,
      panel_xmin=bounds[1],panel_xmax=bounds[2],panel_ymin=bounds[3],panel_ymax=bounds[4],
      measurement=context$status,geometry_metric_status=context$geometry_metric_status,
      converged=context$converged,stringsAsFactors=FALSE)
    d <- layout$display_instances
    if (is.null(d)) d <- data.frame(instance_id=layout$seqs,accver=layout$seqs,key=layout$seqs)
    d$panel_id <- id; instances[[id]] <- d
    # Report actual gtable guide bounds without claiming arbitrary label collision coverage.
    guide_at <- which(grepl("^guide-box",table$layout$name) &
      !vapply(table$grobs,inherits,logical(1),"zeroGrob"))
    guide_bounds <- lapply(guide_at,function(at) {
      g <- table$grobs[[at]]
      c(width=grid::convertWidth(grid::grobWidth(g),"mm",TRUE),height=grid::convertHeight(grid::grobHeight(g),"mm",TRUE))
    })
    oversized <- any(vapply(guide_bounds,function(size) any(size>dimensions),logical(1)))
    if (oversized) ggchord_stop("Fusion panel '",id,"': guides exceed its container; increase size or configure local guides")
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
  # Use the shared Ribbon renderer with a single full-canvas Cartesian frame.
  frame <- ggplot2::ggplot() + ggplot2::coord_cartesian(xlim=c(0,canvas[1]),ylim=c(0,canvas[2]),expand=FALSE)
  frame <- ggplot2::ggplot_build(frame)
  connection_rows <- list(); fronts <- list(); paths <- list(); ribbons <- list(); seen <- list()
  corridors <- list(); corridor_rows <- list()
  for (i in seq_along(x$ribbons)) {
    layer <- x$ribbons[[i]]; input <- layer$data
    layer_id <- paste0("fusion-ribbon-",i)
    for (name in names(layer$selectors)) if (!is.null(layer$selectors[[name]])) input[[name]] <- rep(layer$selectors[[name]],nrow(input))
    unknown <- ggchord_alignment_unknown(input,"fusion_ribbon()")
    if (any(unknown)) warning("Fusion: omitted relationships with unknown direction; see connections",call.=FALSE)
    for (row in seq_len(nrow(input))) {
      record <- input[row,,drop=FALSE]
      report <- cbind(data.frame(layer_id=layer_id,source_row=row),record)
      if (unknown[row]) {
        report$retained <- FALSE; report$reason <- "unknown_direction"
        connection_rows[[length(connection_rows)+1L]] <- report
        next
      }
      q <- ggchord_fusion_endpoint(panels,record,"q",layer_id,row)
      s <- ggchord_fusion_endpoint(panels,record,"s",layer_id,row)
      if (q$panel == s$panel) ggchord_stop("Fusion ",layer_id," row ",row,": cross-panel links require distinct panels")
      if (q$unit != s$unit) ggchord_stop("Fusion ",layer_id," row ",row,": bp/aa links require an explicit source projection")
      ggchord_check_interval_unit(record,q$unit,"fusion_ribbon()")
      if ("source_hsp_id" %in% names(record)) {
        key <- as.character(record$source_hsp_id)
        if (is.na(key) || !nzchar(key)) ggchord_stop("Fusion: source_hsp_id must be nonempty when supplied")
        if (!is.null(seen[[key]]) && !(layer$allow_repeated && all(unlist(seen[[key]]))))
          ggchord_stop("Fusion: repeated source_hsp_id '",key,"' requires allow_repeated=TRUE on all participating layers")
        seen[[key]] <- c(seen[[key]],list(layer$allow_repeated))
      }
      qf <- ggchord_fusion_front(panels[[q$panel]],q); sf <- ggchord_fusion_front(panels[[s$panel]],s)
      path <- ggchord_fusion_path(qf,sf,layer$curvature)
      style <- styles$data[[i]][row,,drop=FALSE]
      relevant <- ggchord_route_obstacles(obstacles, c(q$panel, s$panel))
      stroke <- if (is.na(style$colour) || !is.finite(style$linewidth)) 0 else style$linewidth/2
      relevant <- ggchord_route_inflate(relevant, stroke)
      group <- if (identical(layer$route, "bundle")) record$route_group else NA_character_
      # Length-prefixed keys avoid collisions with slashes or arbitrary user IDs.
      key <- if (!is.na(group)) paste0(nchar(c(layer_id, q$panel, s$panel, group)), ":",
        c(layer_id, q$panel, s$panel, group), collapse = "") else NULL
      corridor <- if (!is.null(key)) corridors[[key]] else NULL
      routed <- ggchord_route_ribbon(qf, sf, path, relevant, canvas,
        layer$route %||% "direct", layer$route_width %||% .8, corridor = corridor, tip = layer$tip)
      if (!is.null(key)) {
        corridors[[key]] <- routed$corridor
        corridor_rows[[length(corridor_rows)+1L]] <- data.frame(layer_id = layer_id,
          route_group = group, source_row = row, qpanel = q$panel, spanel = s$panel,
          vertex = seq_len(nrow(routed$corridor)), x = routed$corridor[, 1], y = routed$corridor[, 2],
          coordinate_space = "composition", units = "mm", semantic = "shared_drawing_corridor")
      }
      path <- routed$path
      normals <- attr(path,"endpoint_normals")
      for (prefix in c("q","s")) {
        endpoint <- if(prefix=="q") q else s
        front <- if(prefix=="q") qf else sf
        front$endpoint <- prefix; front$panel_id <- endpoint$panel; front$instance_id <- endpoint$instance
        front$accver <- endpoint$entity; front$layer_id <- layer_id; front$source_row <- row
        front$sequence_unit <- endpoint$unit
        front$control_normal_x <- normals[prefix,1]; front$control_normal_y <- normals[prefix,2]
        fronts[[length(fronts)+1L]] <- front
      }
      d <- data.frame(x=path[,1],y=path[,2],group=1L,fill=style$fill,colour=style$colour,
        alpha=style$alpha,linewidth=style$linewidth,linetype=style$linetype,
        flow_qn=nrow(qf),flow_sn=nrow(sf),flow_tip=layer$tip,flow_from=NA_character_,flow_to=NA_character_)
      if (layer$gradient) {
        endpoint_colour <- function(endpoint,front) {
          panel <- panels[[endpoint$panel]]
          if (isTRUE(panel$layout$character)) ggchord_stop("Fusion: character endpoints require explicit fill; sequence-colour gradient has no backbone colour")
          colours <- ggchord_sequence_display_colours(panel$built$plot,panel$layout)
          ids <- if(is.null(panel$layout$display_instances)) endpoint$key else endpoint$instance
          at <- which(colours$instance == ids)
          if (!length(at)) ggchord_stop("Fusion: cannot resolve endpoint sequence colour")
          at <- at[which.min((colours$x[at]-mean(front$local_x))^2+(colours$y[at]-mean(front$local_y))^2)]
          colours$colour[at]
        }
        d$flow_from <- endpoint_colour(q,qf); d$flow_to <- endpoint_colour(s,sf)
      }
      ribbons[[length(ribbons)+1L]] <- GeomRibbonFlow$draw_panel(d,frame$layout$panel_params[[1]],frame$layout$coord,
        flow=layer$gradient,tip=layer$tip)
      d$layer_id <- layer_id; d$source_row <- row; paths[[length(paths)+1L]] <- d
      report$route <- layer$route %||% "direct"; report$route_status <- routed$status
      report$obstacle_collisions <- routed$collisions; report$tip <- layer$tip; report$route_stroke_mm <- stroke
      report$qpanel <- q$panel; report$spanel <- s$panel
      report$q_instance <- q$instance; report$s_instance <- s$instance
      report$sequence_unit <- q$unit; report$retained <- TRUE; report$reason <- "retained"
      connection_rows[[length(connection_rows)+1L]] <- report
    }
  }
  for (nested in nested_results) for (key in names(nested$permissions)) {
    permissions <- c(seen[[key]], nested$permissions[[key]])
    if (length(permissions) > 1L && !all(unlist(permissions)))
      ggchord_stop("Fusion: repeated source_hsp_id '", key, "' across nested panels requires allow_repeated=TRUE on all participating layers")
    seen[[key]] <- permissions
  }
  focus_settings <- x$focus
  if (!is.null(focus_settings)) { focus_settings$obstacles <- obstacles; focus_settings$canvas <- canvas }
  focus <- ggchord_build_focus(focus_settings, panels)
  for (name in c("regions", "fronts", "leaders", "leader_paths")) focus[[name]]$composition_id <- rep("root", nrow(focus[[name]]))
  guide_grobs <- lapply(seq_len(nrow(guide_boxes)), function(i) {
    b <- guide_boxes[i, ]
    grid::grobTree(styles$guides[[b$side]], vp = grid::viewport(
      x = mean(c(b$xmin, b$xmax)), y = mean(c(b$ymin, b$ymax)),
      width = b$width_mm, height = b$height_mm, default.units = "mm"))
  })
  guide <- do.call(grid::grobTree, guide_grobs)
  overlaps <- list()
  ids <- names(x$panels)
  for (i in seq_along(ids)) for (j in seq_len(i - 1L)) {
    a <- rectangles[[i]]; b <- rectangles[[j]]
    box <- c(max(a[1], b[1]), min(a[2], b[2]), max(a[3], b[3]), min(a[4], b[4]))
    if (box[1] >= box[2] || box[3] >= box[4]) next
    box <- box * rep(area, each = 2L) + rep(shift, each = 2L)
    overlaps[[length(overlaps) + 1L]] <- data.frame(under = ids[j], over = ids[i],
      xmin = box[1], xmax = box[2], ymin = box[3], ymax = box[4], units = "mm",
      policy = "explicit_allow; later panel paints over earlier panel")
  }
  connection_grobs <- c(ribbons, focus$underlay, connection_grobs)
  overlay_grobs <- c(overlay_grobs, focus$overlay)
  render_events[[length(render_events) + 1L]] <- list(panel_id = NULL,
    background = grid::nullGrob(), foreground = guide, bounds = NULL)
  overlap_mode <- isTRUE(x$allow_overlap) || any(vapply(nested_results,
    function(n) nrow(n$export$overlaps) > 0L, logical(1)))
  if (overlap_mode) {
    containers <- Filter(Negate(is.null), lapply(render_events, `[[`, "container_bounds"))
    clip_regions <- function(regions) lapply(regions, function(bounds)
      ggchord_fusion_clip_connections(connection_grobs, bounds, canvas))
    children <- clip_regions(ggchord_fusion_rect_subtract(c(0,canvas[1],0,canvas[2]), containers))
    for (i in seq_along(render_events)) {
      event <- render_events[[i]]
      children <- c(children, list(event$background))
      if (!is.null(event$panel_id)) {
        later <- if (i < length(render_events)) Filter(Negate(is.null),
          lapply(render_events[seq.int(i+1L,length(render_events))], `[[`, "container_bounds")) else list()
        children <- c(children, clip_regions(ggchord_fusion_rect_subtract(event$bounds, later)))
      }
      children <- c(children, list(event$foreground),
        unname(overlay_grobs[names(overlay_grobs) == event$panel_id]))
    }
  } else children <- c(lapply(render_events, `[[`, "background"), connection_grobs,
    lapply(render_events, `[[`, "foreground"), unname(overlay_grobs))
  grob <- do.call(grid::grobTree, c(unname(children), list(vp = grid::viewport(xscale=c(0,1),yscale=c(0,1)))))
  bind <- function(rows) if(length(rows)) ggchord_rbind_fill(rows) else data.frame()
  child_tables <- function(name) lapply(nested_results, function(n) n$export[[name]])
  combine <- function(rows, name) bind(c(rows, child_tables(name)))
  focus_regions <- bind(c(list(focus$regions), child_tables("focus_regions")))
  focus_fronts <- bind(c(list(focus$fronts), child_tables("focus_fronts")))
  leaders <- bind(c(list(focus$leaders), child_tables("leaders")))
  leader_paths <- bind(c(list(focus$leader_paths), child_tables("leader_paths")))
  for (id in names(panel_rows)) {
    d <- panel_rows[[id]]; at <- which(d$panel_id == id); spec <- x$panels[[id]]
    for (name in c("requested_angle", "requested_distance")) if (!name %in% names(d)) d[[name]] <- NA_real_
    d$requested_angle[at] <- spec$angle %||% NA_real_
    for (name in c("auto_angle", "auto_size_factor")) {
      if (!name %in% names(d)) d[[name]] <- NA_real_
      d[[name]][at] <- attr(rectangles[[id]], name) %||% NA_real_
    }
    d$requested_distance[at] <- spec$distance %||% NA_real_
    if (!"requested_box" %in% names(d)) d$requested_box <- I(rep(list(NULL), nrow(d)))
    d$requested_box[at] <- list(spec$box)
    panel_rows[[id]] <- d
  }
  export <- list(panels=bind(panel_rows),instances=bind(instances),overlaps=combine(overlaps,"overlaps"),
    connections=combine(connection_rows,"connections"),fronts=combine(fronts,"fronts"),paths=combine(paths,"paths"),
    focus_region_data=focus$input,focus_regions=focus_regions,focus_fronts=focus_fronts,leaders=leaders,leader_paths=leader_paths,
    nested_focus_region_data=lapply(nested_results,function(n) n$export$focus_region_data),
    focus_region_data_by_composition=c(list(root=focus$input), unlist(lapply(names(nested_results), function(id) {
      data <- nested_results[[id]]$export$focus_region_data_by_composition
      stats::setNames(data, paste0(id, "/", names(data)))
    }), recursive=FALSE)),
    routing_obstacles=obstacles,routing_corridors=combine(corridor_rows,"routing_corridors"),diagnostics=diagnostics,source_layouts=lapply(panels,`[[`,"layout"),source_data=sources,
    relationship_data=c(lapply(x$ribbons,`[[`,"data"),unlist(lapply(nested_results,function(n) n$export$relationship_data),recursive=FALSE)),
    metadata=list(canvas_mm=canvas,panel_layout=attr(rectangles,"automatic_layout") %||% list(status="fixed"),coordinate_space="composition",units="mm",layout=x$fusion_layout,
      guide_side=styles$side,guide_size_mm=guide_size,guide_sizes_mm=styles$sizes,
      guide_bounds_mm=guide_boxes,guide_layout=styles$report,measurement="measured",allow_overlap=isTRUE(x$allow_overlap),
      device=names(grDevices::dev.cur()),font_conditions="Selected device metrics; fallback and embedding not certified"))
  # A parent can add overlapping leaves after a child was measured. Recheck
  # child geometry in the complete enclosing canvas rather than exporting a
  # stale clear claim. Exact routes fail explicitly if the parent blocks them.
  for (i in seq_len(nrow(export$connections))) {
    connection <- export$connections[i, , drop = FALSE]
    if (!isTRUE(connection$retained)) next
    path <- export$paths[export$paths$layer_id == connection$layer_id &
      export$paths$source_row == connection$source_row, c("x", "y"), drop = FALSE]
    relevant <- ggchord_route_inflate(ggchord_route_obstacles(obstacles, c(connection$qpanel, connection$spanel)), connection$route_stroke_mm)
    qn <- sum(export$fronts$layer_id == connection$layer_id & export$fronts$source_row == connection$source_row & export$fronts$endpoint == "q")
    sn <- sum(export$fronts$layer_id == connection$layer_id & export$fronts$source_row == connection$source_row & export$fronts$endpoint == "s")
    polygon <- ggchord_ribbon_visible_polygon(as.matrix(path), qn, sn, connection$tip)
    collisions <- sum(ggchord_route_polygon_hits(polygon, relevant))
    if (collisions && connection$route %in% c("avoid", "bundle"))
      ggchord_stop("Fusion: enclosing layout obstructs a measured route in ", connection$layer_id,
        "; move the overlapping parent panels or use route = 'direct' with explicit diagnostics")
    export$connections$obstacle_collisions[i] <- collisions
    if (collisions) export$connections$route_status[i] <- "direct_collision"
  }
  for (i in seq_len(nrow(export$leaders))) {
    leader <- export$leaders[i, , drop = FALSE]
    selected <- export$focus_regions$composition_id == leader$composition_id &
      export$focus_regions$region_row == leader$region_row & export$focus_regions$region_piece == leader$region_piece
    region <- export$focus_regions[selected, , drop = FALSE]
    selected <- export$leader_paths$composition_id == leader$composition_id &
      export$leader_paths$region_row == leader$region_row & export$leader_paths$region_piece == leader$region_piece &
      export$leader_paths$boundary == leader$boundary
    line <- as.matrix(export$leader_paths[selected, c("x", "y"), drop = FALSE])
    relevant <- ggchord_route_inflate(ggchord_route_obstacles(obstacles, c(region$source_panel, region$detail_panel)), region$leader_linewidth/2)
    hits <- rep(FALSE, nrow(relevant))
    for (j in seq_len(nrow(line)-1L)) hits <- hits | ggchord_route_segment_hits(line[j, ], line[j+1L, ], relevant)
    if (any(hits) && leader$route == "avoid")
      ggchord_stop("Fusion: enclosing layout obstructs a focus leader; move the overlapping panels or inspect route = 'direct'")
    export$leaders$obstacle_collisions[i] <- sum(hits)
    if (any(hits)) export$leaders$route_status[i] <- "direct_collision"
  }
  for (i in seq_len(nrow(export$focus_regions))) {
    region <- export$focus_regions[i, , drop = FALSE]
    selected <- export$leaders$composition_id == region$composition_id &
      export$leaders$region_row == region$region_row & export$leaders$region_piece == region$region_piece
    states <- export$leaders$route_status[selected]
    export$focus_regions$collision_status[i] <- if (length(states)) paste(unique(states), collapse = ";") else "no_leaders"
  }
  if ("source_hsp_id" %in% names(export$connections)) {
    d <- export$connections
    ids <- unique(as.character(d$source_hsp_id[!is.na(d$source_hsp_id) & d$retained]))
    export$source_observations <- ggchord_rbind_fill(lapply(ids, function(id) data.frame(
      source_hsp_id = id, displayed_connections = sum(d$source_hsp_id %in% id & d$retained),
      identity_basis = "explicit source_hsp_id; display copies are not new observations"))) %||% data.frame()
  } else export$source_observations <- data.frame()
  list(grob=grob,export=export,panel_registry=panels,source_permissions=seen,
    render_events=render_events,connection_grobs=connection_grobs,overlay_grobs=overlay_grobs)
}
