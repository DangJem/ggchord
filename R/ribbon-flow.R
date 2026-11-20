# One decorative tip geometry for rendering and physical route validation.
ggchord_ribbon_tip_edges <- function(edge1, edge2, front, tip) {
  if (tip == "subject") {
    n <- nrow(edge1); center <- (edge1+edge2)/2; span <- edge2[n, ]-edge1[n, ]
    if (sum(span^2) < 1e-20 || sum((center[n, ]-center[1L, ])^2) < 1e-20) tip <- "none" else {
      shoulder <- max(2L, floor(.88*(n-1L))+1L)
      middle <- (nrow(front)+1)/2
      marker <- colMeans(front[unique(c(floor(middle), ceiling(middle))), , drop = FALSE])
      base1 <- center[shoulder, ]-span/2; base2 <- center[shoulder, ]+span/2
      u <- seq(0, 1, length.out = n-shoulder+1L)
      edge1[shoulder:n, ] <- outer(1-u, base1)+outer(u, marker)
      edge2[shoulder:n, ] <- outer(1-u, base2)+outer(u, marker)
    }
  }
  list(edge1 = edge1, edge2 = edge2, tip = tip)
}

ggchord_ribbon_visible_polygon <- function(path, nq, ns, tip) {
  if (tip == "none") return(path)
  ne <- (nrow(path)-nq-ns)/2
  subject <- path[nq+ne+seq_len(ns), , drop = FALSE]
  shape <- ggchord_ribbon_tip_edges(path[nrow(path)-seq_len(ne)+1L, , drop = FALSE],
    path[nq+seq_len(ne), , drop = FALSE], subject, tip)
  rbind(path[seq_len(nq), , drop = FALSE], shape$edge2,
    if (shape$tip == "none") subject, shape$edge1[ne:1L, , drop = FALSE])
}

# Keep the projected interval fronts and compute tip decoration separately.
ggchord_ribbon_trace <- function(data, tip, tip_length, coord, panel, metric) {
  if (!all(c("flow_qn", "flow_sn") %in% names(data))) return(NULL)
  nq <- as.integer(data$flow_qn[1L])
  ns <- as.integer(data$flow_sn[1L])
  ne <- (nrow(data) - nq - ns) / 2
  if (anyNA(c(nq, ns, ne)) || nq < 1 || ns < 1 || ne < 2 || ne != floor(ne)) return(NULL)
  xy <- as.matrix(data[c("x", "y")])
  query <- xy[seq_len(nq), , drop = FALSE]
  subject <- xy[nq + ne + seq_len(ns), , drop = FALSE]
  a <- xy[nrow(xy) - seq_len(ne) + 1L, , drop = FALSE]
  b <- xy[nq + seq_len(ne), , drop = FALSE]
  tips <- list()
  from <- tip %in% c("from", "both") && !identical(data$.arrow_first[1L], FALSE)
  to <- tip %in% c("to", "both") && !identical(data$.arrow_last[1L], FALSE)
  if (from || to) {
    center <- (a + b) / 2
    screen <- coord$transform(data.frame(x = center[, 1L], y = center[, 2L]), panel)
    arc <- c(0, cumsum(sqrt((diff(screen$x) * metric[1L])^2 +
      (diff(screen$y) * metric[2L])^2)))
    total <- tail(arc, 1)
    if (is.finite(total) && total > 0) {
      tip_span <- min(tip_length, total / 2)
      cuts <- c(if (from) tip_span, if (to) total - tip_span)
      samples <- sort(unique(c(arc, cuts)))
      resample <- function(edge) vapply(seq_len(2L), function(axis)
        stats::approx(arc, edge[, axis], xout = samples, ties = "ordered")$y,
        numeric(length(samples)))
      a <- resample(a)
      b <- resample(b)
      add_tip <- function(front, at, distance, side) {
        cut <- if (side == "from") tail(at, 1L) else at[1L]
        marker <- colMeans(front)
        a[at, ] <<- outer(1 - distance, a[cut, ]) + outer(distance, marker)
        b[at, ] <<- outer(1 - distance, b[cut, ]) + outer(distance, marker)
        tips[[side]] <<- rbind(a[at, , drop = FALSE], b[rev(at), , drop = FALSE])
      }
      if (from) {
        at <- which(samples <= tip_span)
        add_tip(query, at, 1 - samples[at] / tip_span, "from")
        query <- query[FALSE, , drop = FALSE]
      }
      if (to) {
        at <- which(samples >= total - tip_span)
        add_tip(subject, at, (samples[at] - (total - tip_span)) / tip_span, "to")
        subject <- subject[FALSE, , drop = FALSE]
      }
    }
  }
  list(polygon = rbind(query, b, subject, a[rev(seq_len(nrow(a))), , drop = FALSE]),
    edge1 = a, edge2 = b, tips = tips)
}

# Preserve complete fronts in exported data and add separate decoration at draw time.
GeomRibbonFlow <- ggplot2::ggproto("GeomRibbonFlow", ggplot2::GeomPolygon,
  default_aes = aes(fill = "grey70", colour = NA, alpha = .6, linewidth = .3,
    linetype = 1, fill_from = NA_character_, fill_to = NA_character_,
    flow_qn = 0, flow_sn = 0),
  extra_params = c("na.rm", "flow", "tip", "tip_length"),
  draw_panel = function(data, panel_params, coord, na.rm = FALSE, rule = "evenodd",
      flow = FALSE, tip = "none", tip_length = 2) {
    if (!flow && tip == "none")
      return(ggplot2::GeomPolygon$draw_panel(data, panel_params, coord, rule = rule))
    groups <- split(data, data$group)
    device_pixels <- if (grDevices::dev.cur() == 1L) 1e6 else prod(grDevices::dev.size("px"))
    pixel_budget <- min(8e6, max(2e6, device_pixels * 4)) / max(1, length(groups))
    metric <- c(grid::convertWidth(grid::unit(1, "npc"), "mm", TRUE),
      grid::convertHeight(grid::unit(1, "npc"), "mm", TRUE))
    grobs <- lapply(groups, function(d) {
      trace <- ggchord_ribbon_trace(d, tip, tip_length, coord, panel_params, metric)
      if (is.null(trace))
        return(ggplot2::GeomPolygon$draw_panel(d, panel_params, coord, rule = rule))
      xy <- trace$polygon
      border <- d[rep(1L, nrow(xy)), , drop = FALSE]
      border$x <- xy[, 1L]; border$y <- xy[, 2L]
      if (!flow) return(ggplot2::GeomPolygon$draw_panel(border, panel_params, coord, rule = rule))
      endpoints <- c(d$fill_from[1L], d$fill_to[1L])
      if (anyNA(endpoints)) return(grid::nullGrob())
      opacity <- grDevices::col2rgb(endpoints, alpha = TRUE)[4L, ] / 255
      n <- max(2L, nrow(trace$edge1) - 1L)
      colors <- scales::alpha(grDevices::colorRampPalette(endpoints, space = "Lab")(n),
        seq(opacity[1L], opacity[2L], length.out = n))
      fill <- ggchord_ribbon_mesh_fill(trace$edge1, trace$edge2, xy, colors, d$alpha[1L],
        coord, panel_params, pixel_budget)
      border$fill <- NA
      grid::grobTree(fill,
        ggplot2::GeomPolygon$draw_panel(border, panel_params, coord, rule = rule))
    })
    do.call(grid::grobTree, grobs)
  },
  draw_key = function(data, params, size) {
    if (isTRUE(params$flow)) {
      endpoints <- c(data$fill_from, data$fill_to)
      data$fill <- c(endpoints[!is.na(endpoints)], "grey70")[1L]
    }
    ggplot2::draw_key_polygon(data, params, size)
  })

ggchord_ribbon_flow <- function(layer, geometry, layout, plot) {
  if (!identical(layer$ggchord_type, "ribbon") || !nrow(geometry))
    return(list(layer = layer, data = geometry))
  flow <- identical(layer$ggchord_params$fill_mode, "endpoint")
  mapping <- layer$ggchord_input_mapping
  channels <- c(names(mapping), names(layer$aes_params))
  if (!flow && any(c("fill_from", "fill_to") %in% channels))
    ggchord_stop("Endpoint colours require fill_mode = 'endpoint'")
  if (flow) {
    if ("fill" %in% channels) ggchord_stop("Endpoint ribbons use fill_from and fill_to")
    scopes <- ggchord_scale_scopes(plot)
    for (channel in c("fill_from", "fill_to")) {
      if (channel %in% names(layer$aes_params)) next
      scoped <- any(vapply(scopes, function(scope) {
        channel %in% scope$aesthetics && any(vapply(plot$layers[scope$selected],
          function(candidate) identical(candidate$ggchord_layer_id, layer$ggchord_layer_id), logical(1)))
      }, logical(1)))
      if (!channel %in% names(mapping) || (!plot$scales$has_scale(channel) && !scoped))
        ggchord_stop("Endpoint ribbons require a fixed colour or an explicit scale for ", channel)
    }
  }
  layer$mapping$flow_qn <- aes(flow_qn = I(.q_n))$flow_qn
  layer$mapping$flow_sn <- aes(flow_sn = I(.s_n))$flow_sn
  list(layer = layer, data = geometry)
}

# Paint the mesh once so translucent strips cannot accumulate edge seams.
# Geometry/outlines remain the original paths; only the continuous fill is raster.
ggchord_ribbon_mesh_fill <- function(edge1,edge2,outline,colours,alpha,coord,panel_params,pixel_budget=8e6) {
  transform <- function(xy) {
    d <- coord$transform(data.frame(x=xy[,1],y=xy[,2]),panel_params)
    cbind(x=d$x,y=d$y)
  }
  a<-transform(edge1);b<-transform(edge2);polygon<-transform(outline)
  limits<-c(range(polygon[,1]),range(polygon[,2]))
  canvas<-grDevices::dev.size('in');pixels<-grDevices::dev.size('px')
  resolution<-pmax(150,pixels/canvas)
  panel_inches<-c(grid::convertWidth(grid::unit(1,'npc'),'in',TRUE),
    grid::convertHeight(grid::unit(1,'npc'),'in',TRUE))
  dimensions<-pmax(2L,pmin(1600L,ceiling(c(diff(limits[1:2]),diff(limits[3:4]))*panel_inches*resolution)))
  reduced <- prod(dimensions)>pixel_budget
  if(reduced)dimensions <- pmax(2L,floor(dimensions*sqrt(pixel_budget/prod(dimensions))))
  nx<-dimensions[1];ny<-dimensions[2]
  xx<-seq(limits[1],limits[2],length.out=nx+1L);xx<-(head(xx,-1)+tail(xx,-1))/2
  yy<-seq(limits[3],limits[4],length.out=ny+1L);yy<-(head(yy,-1)+tail(yy,-1))/2
  # Outside the two sampled side edges only the true curved endpoint cap can
  # remain. Its colour is the corresponding endpoint, not screen-x colour.
  q<-colMeans(rbind(a[1,],b[1,]));s<-colMeans(rbind(tail(a,1),tail(b,1)))
  values<-outer(yy,xx,function(y,x) as.numeric((x-s[1])^2+(y-s[2])^2 < (x-q[1])^2+(y-q[2])^2))
  t<-seq(0,1,length.out=nrow(a))
  trace<-function(p,v) {
    ix<-which(xx>=min(p[,1]) & xx<=max(p[,1]));iy<-which(yy>=min(p[,2]) & yy<=max(p[,2]))
    if(!length(ix)||!length(iy))return()
    det<-(p[2,2]-p[3,2])*(p[1,1]-p[3,1])+(p[3,1]-p[2,1])*(p[1,2]-p[3,2])
    if(abs(det)<1e-16)return()
    x<-rep(xx[ix],each=length(iy));y<-rep(yy[iy],times=length(ix))
    w1<-((p[2,2]-p[3,2])*(x-p[3,1])+(p[3,1]-p[2,1])*(y-p[3,2]))/det
    w2<-((p[3,2]-p[1,2])*(x-p[3,1])+(p[1,1]-p[3,1])*(y-p[3,2]))/det
    w3<-1-w1-w2;inside<-w1>=-1e-8 & w2>=-1e-8 & w3>=-1e-8
    at<-cbind(rep(iy,times=length(ix)),rep(ix,each=length(iy)))
    values[at[inside,,drop=FALSE]]<<-(w1*v[1]+w2*v[2]+w3*v[3])[inside]
  }
  for(i in seq_len(nrow(a)-1L)) {
    trace(rbind(a[i,],a[i+1,],b[i+1,]),c(t[i],t[i+1],t[i+1]))
    trace(rbind(a[i,],b[i+1,],b[i,]),c(t[i],t[i+1],t[i]))
  }
  mask<-matrix(FALSE,ny,nx)
  closed<-rbind(polygon,polygon[1,,drop=FALSE])
  for(j in seq_len(ny)) {
    y<-yy[j];cross<-which((head(closed[,2],-1)>y)!=(tail(closed[,2],-1)>y))
    if(!length(cross))next
    cuts<-sort(closed[cross,1]+(y-closed[cross,2])*
      (closed[cross+1,1]-closed[cross,1])/(closed[cross+1,2]-closed[cross,2]))
    mask[j,]<-findInterval(xx,cuts)%%2L==1L
  }
  ramp<-if(all(is.na(colours)))rep(NA_character_,256) else
    grDevices::colorRampPalette(colours,space='Lab')(256)
  index<-pmax(1L,pmin(256L,round(values*255)+1L))
  opacity<-if(all(is.na(colours))) rep(0,256) else stats::approx(seq(0,1,length.out=length(colours)),
    grDevices::col2rgb(colours,alpha=TRUE)[4,]/255,xout=seq(0,1,length.out=256))$y
  raster<-matrix(scales::alpha(ramp[index],(if(is.na(alpha))1 else alpha)*opacity[index]),ny,nx)
  raster[!mask]<-NA_character_
  grob <- grid::rasterGrob(as.raster(raster[ny:1,,drop=FALSE]),
    x=mean(limits[1:2]),y=mean(limits[3:4]),width=diff(limits[1:2]),height=diff(limits[3:4]),
    default.units='npc',interpolate=TRUE)
  device_dimensions <- pmax(2L,ceiling(c(diff(limits[1:2]),diff(limits[3:4]))*panel_inches*(pixels/canvas)))
  attr(grob,"ggchord_reduced_fill") <- any(dimensions < .8*device_dimensions)
  grob
}
