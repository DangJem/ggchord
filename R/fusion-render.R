# Render once inside final panel containers; never freeze a previous-size image.
ggchord_build_fusion <- function(x) {
  canvas <- c(width=grid::convertWidth(grid::unit(1,"npc"),"mm",TRUE),
    height=grid::convertHeight(grid::unit(1,"npc"),"mm",TRUE))
  origin <- grid::deviceLoc(grid::unit(0,"npc"),grid::unit(0,"npc"),valueOnly=TRUE)
  origin <- c(origin$x,origin$y)*25.4
  styles <- ggchord_fusion_styles(x)
  guide_size <- if (styles$side=="none") c(0,0) else c(
    grid::convertWidth(grid::grobWidth(styles$guide),"mm",TRUE),
    grid::convertHeight(grid::grobHeight(styles$guide),"mm",TRUE))
  area <- canvas; shift <- c(0,0)
  if (styles$side %in% c("left","right")) {
    area[1] <- canvas[1]-guide_size[1]-4
    if (styles$side=="left") shift[1] <- guide_size[1]+4
  }
  if (styles$side %in% c("top","bottom")) {
    area[2] <- canvas[2]-guide_size[2]-4
    if (styles$side=="bottom") shift[2] <- guide_size[2]+4
  }
  if (min(area)<15 || any(guide_size > canvas)) ggchord_stop(
    "Fusion: canvas cannot contain required cross-panel guides; increase width/height or set explicit guide layout")
  rectangles <- ggchord_fusion_rectangles(x, area)
  panels <- list(); panel_rows <- list(); instances <- list(); diagnostics <- list()
  backgrounds <- list(); foregrounds <- list()
  for (id in names(x$panels)) {
    spec <- x$panels[[id]]; rect <- rectangles[[id]]
    dimensions <- c(diff(rect[1:2]),diff(rect[3:4]))*area
    center <- c(mean(rect[1:2]),mean(rect[3:4]))*area+shift
    viewport <- grid::viewport(x=center[1],y=center[2],width=dimensions[1],height=dimensions[2],
      default.units="mm",clip="off",name=paste0("fusion-",id))
    # Deep clone layers as well: extension geoms may carry their own state.
    plot <- unserialize(serialize(ggchord_render_plot_copy(spec$plot),NULL))
    table <- local({
      grid::pushViewport(viewport)
      on.exit(grid::popViewport())
      ggchord_with_render_context(list(available_canvas_mm=dimensions),ggplot2::ggplotGrob(plot))
    })
    context <- attr(table,"ggchord_render_context")
    if (!identical(context$status,"measured")) ggchord_stop("Fusion panel '",id,"': actual panel measurement unavailable")
    context$device_panel_bounds_mm <- context$panel_bounds_mm
    context$panel_bounds_mm <- context$panel_bounds_mm-rep(origin,each=2L)
    layout <- attr(table,"ggchord_render_layout")
    built <- attr(table,"ggchord_render_built")
    panels[[id]] <- list(layout=layout,built=built,context=context,
      seq_data=spec$plot$ggchord$data$seq_data)
    backgrounds[[id]] <- grid::grobTree(ggchord_fusion_partition(table,TRUE),vp=viewport)
    foregrounds[[id]] <- grid::grobTree(ggchord_fusion_partition(table,FALSE),vp=viewport)
    bounds <- context$panel_bounds_mm
    panel_rows[[id]] <- data.frame(panel_id=id,side=spec$side,
      coordinate=if(isTRUE(spec$plot$coordinates$ggchord_circular)) "circular" else
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
      label_fit=ggchord_label_fit_diagnostics(layout),label_visibility=ggchord_label_visibility(layout),
      guide_sizes_mm=guide_bounds,guide_overflow=oversized,collision_status="not globally measured")
  }
  # Use the shared Ribbon renderer with a single full-canvas Cartesian frame.
  frame <- ggplot2::ggplot() + ggplot2::coord_cartesian(xlim=c(0,canvas[1]),ylim=c(0,canvas[2]),expand=FALSE)
  frame <- ggplot2::ggplot_build(frame)
  connection_rows <- list(); fronts <- list(); paths <- list(); ribbons <- list(); seen <- list()
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
      style <- styles$data[[i]][row,,drop=FALSE]
      d <- data.frame(x=path[,1],y=path[,2],group=1L,fill=style$fill,colour=style$colour,
        alpha=style$alpha,linewidth=style$linewidth,linetype=style$linetype,
        flow_qn=nrow(qf),flow_sn=nrow(sf),flow_tip=layer$tip,flow_from=NA_character_,flow_to=NA_character_)
      if (layer$gradient) {
        endpoint_colour <- function(endpoint,front) {
          panel <- panels[[endpoint$panel]]
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
      report$qpanel <- q$panel; report$spanel <- s$panel
      report$q_instance <- q$instance; report$s_instance <- s$instance
      report$sequence_unit <- q$unit; report$retained <- TRUE; report$reason <- "retained"
      connection_rows[[length(connection_rows)+1L]] <- report
    }
  }
  guide <- styles$guide
  if (styles$side!="none") {
    xy <- switch(styles$side,right=c(canvas[1]-guide_size[1]/2,canvas[2]/2),
      left=c(guide_size[1]/2,canvas[2]/2),top=c(canvas[1]/2,canvas[2]-guide_size[2]/2),
      bottom=c(canvas[1]/2,guide_size[2]/2))
    guide <- grid::grobTree(guide,vp=grid::viewport(x=xy[1],y=xy[2],width=guide_size[1],height=guide_size[2],default.units="mm"))
  }
  grob <- do.call(grid::grobTree,c(unname(backgrounds),ribbons,unname(foregrounds),
    list(guide, vp=grid::viewport(xscale=c(0,1),yscale=c(0,1)))))
  bind <- function(rows) if(length(rows)) ggchord_rbind_fill(rows) else data.frame()
  list(grob=grob,export=list(panels=bind(panel_rows),instances=bind(instances),
    connections=bind(connection_rows),fronts=bind(fronts),paths=bind(paths),
    diagnostics=diagnostics,source_layouts=lapply(panels,`[[`,"layout"),
    source_data=lapply(x$panels,function(p)p$plot$ggchord$sources),
    relationship_data=lapply(x$ribbons,`[[`,"data"),
    metadata=list(canvas_mm=canvas,coordinate_space="composition",units="mm",layout=x$fusion_layout,
      guide_side=styles$side,guide_size_mm=guide_size,measurement="measured",
      device=names(grDevices::dev.cur()),font_conditions="Selected device metrics; fallback and embedding not certified")))
}
