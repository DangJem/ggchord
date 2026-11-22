#' Project reference relationships onto two aligned sets of paths
#'
#' @inheritParams coord_collinear
#' @return A coordinate system for a ggchord plot.
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
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' ggchord(data = seqs) +
#'   geom_seq() +
#'   geom_gene(data = genes) +
#'   geom_link_ribbon(data = links) +
#'   layout_sequences(column_gap = 25) +
#'   coord_sankey()
coord_sankey <- function(axis = "y", ratio = NULL, xlim = NULL, ylim = NULL,
    expand = FALSE, clip = "on", fit = "geometry") {
  ggchord_check_argument_names(sys.call(), sys.function())
  base <- coord_collinear(axis = axis, ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = clip, fit = fit)
  ggproto("CoordSankey", base, ggchord_sankey = TRUE,
    sankey_self_links = TRUE, sankey_axis = axis, collinear_seq_axis = "sankey")
}

ggchord_sankey_points <- function(x, y, id, frames) {
  at <- match(id, frames$key)
  normal <- frames$column_x[at] + frames$normal_sign[at] * (y - frames$source_y[at])
  along <- frames$top[at] - (x - frames$axis_min[at])
  if (isTRUE(attr(frames, "horizontal"))) cbind(x = -along, y = -normal) else
    cbind(x = normal, y = along)

}

ggchord_project_sankey_data <- function(data, frames, id = NULL) {
  if (!is.data.frame(data) || !nrow(data)) return(data)
  ids <- data$accver %||% rep(id, nrow(data))
  if (is.null(ids) || !length(ids)) return(data)
  for (pair in list(c("x", "y"), c("x0", "y0"), c("x1", "y1"), c("xend", "yend"),
      c("text_x", "text_y"), c("label_x", "label_y"), c("anchor_x", "anchor_y"))) {
    if (!all(pair %in% names(data))) next
    xy <- ggchord_sankey_points(data[[pair[1]]], data[[pair[2]]], ids, frames)
    data[[pair[1]]] <- xy[,1]; data[[pair[2]]] <- xy[,2]
  }
  angles <- intersect(c("angle", "text_angle", "label_angle"), names(data))
  if (isTRUE(attr(frames, "horizontal"))) return(data)
  for (field in angles) {
    angle <- data[[field]] * pi/180
    degrees <- atan2(-cos(angle), sin(angle)) * 180/pi
    above <- which(degrees > 90)
    degrees[above] <- degrees[above] - 180
    below <- which(degrees < -90)
    degrees[below] <- degrees[below] + 180
    data[[field]] <- degrees
  }
  data
}

ggchord_sankey_route <- function(a, b, route, horizontal) {
  if (identical(route$route, "straight")) return(rbind(a, b))
  normal <- if (horizontal) 2L else 1L
  delta <- (b[normal] - a[normal]) * ggchord_connection_bend(route)
  first <- a; first[normal] <- first[normal] + delta
  second <- b; second[normal] <- second[normal] - delta
  bezier_pts(a, b, first, second, n = 50)
}

ggchord_project_sankey_layout <- function(layout, plot) {
  if (isTRUE(layout$sankey_projected)) return(layout)
  instances <- layout$display_instances
  ref <- layout$sequence_reference
  tab <- layout$sequence_layout
  left <- which(instances$side == "left")
  span <- setNames(tab$xmax-tab$xmin, tab$accver)
  gaps <- if (length(left) > 1L) abs(diff(unname(ref$y[instances$key[left]]))) else numeric()
  top <- c(0, -cumsum(head(unname(span[instances$key[left]]), -1) + gaps))
  frames <- instances
  horizontal <- identical(plot$coordinates$sankey_axis, "x")
  attr(frames, "horizontal") <- horizontal
  frames$top <- rep(top, 2)
  gap <- plot$coordinates$sankey_column_gap_mm * ggchord_layout_units_per_mm(if (horizontal) 2L else 1L)
  frames$column_x <- ifelse(frames$side == "left", 0, gap)
  frames$normal_sign <- if (horizontal) -1 else 1
  facing <- setNames(ifelse(frames$side == "left", 1, -1) * frames$normal_sign,
    frames$key)
  frames$source_y <- unname(ref$y[frames$key])
  frames$axis_min <- tab$xmin[match(frames$key,tab$accver)]
  normal_mm <- ggchord_units_per_mm(ggchord_path_normal_axis(plot$coordinates))
  for (owner in names(layout$layer_geometry)) for (component in names(layout$layer_geometry[[owner]])) {
    d <- layout$layer_geometry[[owner]][[component]]
    if (!component %in% c("ribbon", "link")) {
      projected <- ggchord_project_sankey_data(d,frames)
      if (!horizontal && component == "seq_label" && is.data.frame(projected) && nrow(projected)) {
        automatic <- !(projected$.ggchord_vjust_explicit %||% rep(FALSE, nrow(projected)))
        # Projection changes a left-of-row name into an above-column name.
        # Anchor its lower edge there, clear of the first coordinate tick.
        projected$vjust[automatic] <- 0
      }
      layout$layer_geometry[[owner]][[component]] <- projected
      next
    }
    input <- layout$layer_inputs[[owner]][[component]]
    if (!nrow(d) || is.null(input)) next
    layer <- plot$layers[[which(vapply(plot$layers,function(l) identical(l$ggchord_layer_id,owner),logical(1)))[1]]]
    params <- layer$ggchord_params
    if (!identical(params$link_branch %||% "none", "none"))
      ggchord_stop("coord_sankey(): route = 'shared' is unsupported; use route = 'curve'")
    configured_gap <- params$ribbon_gap %||% params$link_gap
    gaps_by_id <- process_sequence_param(configured_gap %||% .035, layout$seqs, "link_gap", .035)
    mapping <- if (isTRUE(layer$inherit.aes)) plot$mapping else NULL
    ribbon_gap_mm <- ggchord_ribbon_gap_mm(params, layout$seqs, layer, mapping)
    along_mm <- ggchord_units_per_mm(if (horizontal) "x" else "y")
    avoid <- if (is.null(configured_gap)) params$link_avoid %||% "none" else "none"
    obstacles <- ggchord_obstacle_index(layout$obstacles, avoid)
    offsets <- ggchord_connection_offsets(params, layout$seqs,
      ggchord_path_normal_axis(plot$coordinates))
    front <- function(id, bounds, samples) {
      direction <- facing[[id]]
      base <- if (is.null(configured_gap))
        ggchord_ribbon_base_gap(layout$ribbon_clearance, id, direction, normal_mm,
          ribbon_gap_mm[[id]]) else
        gaps_by_id[[id]]
      gap_mm <- if (component == "ribbon") ribbon_gap_mm[[id]] else .30
      profile <- ggchord_connection_front(bounds, obstacles[[id]], direction, avoid,
        base, normal_mm, ggchord_collinear_scale(id, ref) * normal_mm / along_mm,
        n = samples, gap_mm = gap_mm, strength = params$avoid_strength %||% .5)
      ggchord_sankey_points(ggchord_collinear_x(profile$position, id, ref),
        unname(ref$y[id]) + direction * profile$gap + offsets[[id]], id, frames)
    }
    pieces <- list()
    for (i in seq_len(nrow(input))) {
      template <- d[d$source_row == i, , drop=FALSE]
      if (!nrow(template)) next
      q <- as.character(input$qaccver[i]); s <- as.character(input$saccver[i])
      if (component == "ribbon") {
        qxy <- front(q, ggchord_interval_front(input$qstart[i], input$qend[i]), 2L)
        sxy <- front(s, ggchord_interval_front(input$sstart[i], input$send[i]), 2L)
        edge <- function(a, b) ggchord_sankey_route(a, b, params$route, horizontal)
        a <- edge(qxy[1,],sxy[1,]); b <- edge(qxy[nrow(qxy),],sxy[nrow(sxy),])
        xy <- rbind(qxy,b,sxy[nrow(sxy):1,,drop=FALSE],a[nrow(a):1,,drop=FALSE])
      } else {
        a <- front(q, rep(input$qpos[i], 2L), 1L)[1L, ]
        b <- front(s, rep(input$spos[i], 2L), 1L)[1L, ]
        xy <- ggchord_sankey_route(a, b, params$route, horizontal)

      }
      piece <- template[rep(1L,nrow(xy)),,drop=FALSE]
      piece$x <- xy[,1]; piece$y <- xy[,2]; piece$group <- i
      piece$.ggchord_adjacent <- TRUE; piece$linetype_val <- "solid"
      if(component=="ribbon") {piece$.q_n<-nrow(qxy);piece$.s_n<-nrow(sxy)}
      pieces[[length(pieces)+1L]] <- piece
    }
    layout$layer_geometry[[owner]][[component]] <- ggchord_rbind_fill(pieces) %||% d[FALSE,,drop=FALSE]
  }
  layout$seq_arcs <- lapply(layout$seq_arcs,ggchord_project_sankey_data,frames=frames)
  for (field in c("axis_lines","axis_ticks","gene_polys","gene_labels","gene_label_segments",
      "region_polys","seq_labels_df","sequence_scales","track_axes"))
    layout[[field]] <- ggchord_project_sankey_data(layout[[field]],frames)
  for (pair in list(c("ribbon_polys","ribbon"),c("link_lines","link"))) {
    parts <- lapply(layout$layer_geometry,function(owner)owner[[pair[2]]])
    layout[[pair[1]]] <- ggchord_rbind_fill(parts)
  }
  tab$axis_min <- tab$xmin; tab$axis_max <- tab$xmax
  at <- match(tab$accver,frames$key)
  tab$ymax <- frames$top[at]; tab$ymin <- tab$ymax - unname(span[tab$accver])
  tab$xmin <- tab$xmax <- frames$column_x[at]; tab$y <- (tab$ymax+tab$ymin)/2
  layout$sequence_layout <- tab
  geometry <- c(unlist(layout$layer_geometry,recursive=FALSE), list(layout$axis_lines,layout$axis_ticks))
  xx <- yy <- numeric()
  for(d in geometry) if(is.data.frame(d)) for(pair in list(c("x","y"),c("x0","y0"),c("x1","y1"),c("text_x","text_y"),c("label_x","label_y")))
    if (all(pair %in% names(d))) {xx<-c(xx,d[[pair[1]]]);yy<-c(yy,d[[pair[2]]])}
  layout$extremes <- list(x_min=min(xx,na.rm=TRUE),x_max=max(xx,na.rm=TRUE),
    y_min=min(yy,na.rm=TRUE),y_max=max(yy,na.rm=TRUE))
  layout$sankey_frames <- frames
  layout$seq_axis <- if (horizontal) "horizontal" else "vertical"
  layout$sankey_projected <- TRUE
  layout
}
