#' Mirrored two-column genomic interval relationships
#'
#' Each source sequence appears once in each column. Input query intervals
#' connect from the left to subject intervals on the right; input reciprocal
#' records remain distinct. This is not a multistage flow or a conservation law.
#' Source lengths and both HSP endpoints determine the displayed interval widths.
#' @param column_gap Positive separation of the two backbones in layout units.
#' @param self_links Show relationships whose query and subject entity are equal.
#' @inheritParams coord_collinear
#' @return A ggchord coordinate object.
#' @export
#' @examples
#' # Synthetic intervals illustrate the mirrored layout.
#' s <- data.frame(accver = c("A", "B"), length = c(1000, 600))
#' h <- data.frame(qaccver = "A", saccver = "B",
#'   qstart = 100, qend = 300, sstart = 50, send = 250)
#' ggchord(s, h) + geom_seq() + geom_link_ribbon() +
#'   geom_seq_label() + coord_sankey()
coord_sankey <- function(column_gap = 3, self_links = FALSE, ratio = 1,
    xlim = NULL, ylim = NULL, expand = FALSE, clip = "off",
    fit = c("labels", "geometry", "manual")) {
  if (!is.numeric(column_gap) || length(column_gap) != 1L || !is.finite(column_gap) || column_gap <= 0)
    ggchord_stop("coord_sankey(): column_gap must be positive")
  if (!is.logical(self_links) || length(self_links) != 1L || is.na(self_links)) ggchord_stop("self_links must be TRUE or FALSE")
  base_coord <- coord_collinear(ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = clip, fit = fit)
  ggplot2::ggproto("CoordSankey", base_coord, ggchord_sankey = TRUE,
    sankey_column_gap = column_gap, sankey_self_links = self_links,
    collinear_seq_axis = "sankey")
}

ggchord_sankey_points <- function(x, y, id, frames) {
  at <- match(id, frames$key)
  cbind(x = frames$column_x[at] + frames$normal_sign[at] * (y - frames$source_y[at]),
    y = frames$top[at] - (x - frames$axis_min[at]))
}

ggchord_project_sankey_data <- function(data, frames, id = NULL) {
  if (!is.data.frame(data) || !nrow(data)) return(data)
  ids <- data$accver %||% rep(id, nrow(data))
  if (is.null(ids) || !length(ids)) return(data)
  for (pair in list(c("x", "y"), c("x0", "y0"), c("x1", "y1"),
      c("text_x", "text_y"), c("label_x", "label_y"), c("anchor_x", "anchor_y"))) {
    if (!all(pair %in% names(data))) next
    xy <- ggchord_sankey_points(data[[pair[1]]], data[[pair[2]]], ids, frames)
    data[[pair[1]]] <- xy[,1]; data[[pair[2]]] <- xy[,2]
  }
  at <- match(ids, frames$key)
  for (field in intersect(c("angle", "text_angle", "label_angle"), names(data))) {
    angle <- data[[field]] * pi/180
    degrees <- atan2(-cos(angle), frames$normal_sign[at] * sin(angle)) * 180/pi
    degrees[degrees > 90] <- degrees[degrees > 90] - 180
    degrees[degrees < -90] <- degrees[degrees < -90] + 180
    data[[field]] <- degrees
  }
  data
}

ggchord_project_sankey_layout <- function(layout, plot) {
  if (isTRUE(layout$sankey_projected)) return(layout)
  instances <- layout$display_instances
  ref <- layout$sequence_reference
  tab <- layout$sequence_layout
  left <- which(instances$side == "left")
  span <- stats::setNames(tab$xmax-tab$xmin, tab$accver)
  gaps <- if (length(left) > 1L) abs(diff(unname(ref$y[instances$key[left]]))) else numeric()
  top <- c(0, -cumsum(head(unname(span[instances$key[left]]), -1) + gaps))
  frames <- instances
  frames$top <- rep(top, 2)
  frames$column_x <- ifelse(frames$side == "left", 0, plot$coordinates$sankey_column_gap)
  frames$normal_sign <- ifelse(frames$side == "left", -1, 1)
  frames$source_y <- unname(ref$y[frames$key])
  frames$axis_min <- tab$xmin[match(frames$key,tab$accver)]
  for (owner in names(layout$layer_geometry)) for (component in names(layout$layer_geometry[[owner]])) {
    d <- layout$layer_geometry[[owner]][[component]]
    if (!component %in% c("ribbon", "link")) {
      layout$layer_geometry[[owner]][[component]] <- ggchord_project_sankey_data(d,frames)
      next
    }
    input <- layout$layer_inputs[[owner]][[component]]
    if (!nrow(d) || is.null(input)) next
    layer <- plot$layers[[which(vapply(plot$layers,function(l) identical(l$ggchord_layer_id,owner),logical(1)))[1]]]
    params <- layer$ggchord_params
    if (!identical(params$link_branch %||% "none", "none"))
      ggchord_stop("coord_sankey(): shared branches are outside the first mirrored layout; use link_branch='none'")
    configured_gap <- params$ribbon_gap %||% params$link_gap
    gaps_by_id <- process_sequence_param(configured_gap %||% .035, layout$seqs, "link_gap", .035)
    pieces <- list()
    for (i in seq_len(nrow(input))) {
      template <- d[d$source_row == i, , drop=FALSE]
      if (!nrow(template)) next
      q <- as.character(input$qaccver[i]); s <- as.character(input$saccver[i])
      if (component == "ribbon") {
        front <- function(id,start,end) {
          positions <- seq(start,end,length.out=16)
          gap <- ggchord_collinear_gap_profile(positions,id,-1,layout$obstacles,
            avoid=if (is.null(configured_gap)) params$link_avoid %||% "none" else "none",
            close_gap=unname(gaps_by_id[id]))
          ggchord_sankey_points(ggchord_collinear_x(positions,id,ref),unname(ref$y[id])-gap,id,frames)
        }
        qxy <- front(q,input$qstart[i],input$qend[i]); sxy <- front(s,input$sstart[i],input$send[i])
        edge <- function(a,b,which_edge) {
          control <- params$ribbon_ctrl_point
          cp <- if (is.null(control)) list(c(a[1]+(b[1]-a[1])*.45,a[2]),
            c(b[1]-(b[1]-a[1])*.45,b[2])) else {
              specified <- ggchord_link_controls(control,i,(a+b)/2)
              rep(specified[which_edge],2)
            }
          bezier_pts(a,b,cp[[1]],cp[[2]],n=50)
        }
        a <- edge(qxy[1,],sxy[1,],1); b <- edge(qxy[nrow(qxy),],sxy[nrow(sxy),],2)
        xy <- rbind(qxy,b,sxy[nrow(sxy):1,,drop=FALSE],a[nrow(a):1,,drop=FALSE])
      } else {
        a <- ggchord_sankey_points(ggchord_collinear_x(input$qpos[i],q,ref),ref$y[q]-gaps_by_id[q],q,frames)[1,]
        b <- ggchord_sankey_points(ggchord_collinear_x(input$spos[i],s,ref),ref$y[s]-gaps_by_id[s],s,frames)[1,]
        xy <- if (identical(layer$ggchord_params$link_type,"straight")) rbind(a,b) else
          { cp <- if (is.null(params$link_ctrl_point)) list(c(a[1]+(b[1]-a[1])*.45,a[2]),
            c(b[1]-(b[1]-a[1])*.45,b[2])) else ggchord_link_controls(params$link_ctrl_point,i,(a+b)/2)
            bezier_pts(a,b,cp[[1]],cp[[2]],n=50) }
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
    if(all(pair %in% names(d))) {xx<-c(xx,d[[pair[1]]]);yy<-c(yy,d[[pair[2]]])}
  layout$extremes <- list(x_min=min(xx,na.rm=TRUE),x_max=max(xx,na.rm=TRUE),
    y_min=min(yy,na.rm=TRUE),y_max=max(yy,na.rm=TRUE))
  layout$sankey_frames <- frames
  layout$seq_axis <- "vertical"; layout$sankey_projected <- TRUE
  layout
}
