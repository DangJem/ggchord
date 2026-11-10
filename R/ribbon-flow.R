# Ribbon colour follows the sampled Bezier parameter, never screen x/y.
GeomRibbonFlow <- ggplot2::ggproto("GeomRibbonFlow", ggplot2::GeomPolygon,
  default_aes = utils::modifyList(ggplot2::GeomPolygon$default_aes,
    ggplot2::aes(flow_from=NA_character_, flow_to=NA_character_, flow_qn=0, flow_sn=0, flow_tip=NA_character_)),
  extra_params = c("na.rm", "flow", "tip"),
  draw_panel = function(data, panel_params, coord, na.rm=FALSE, rule="evenodd",
      flow=FALSE, tip="none") {
    if (!flow && tip == "none") return(ggplot2::GeomPolygon$draw_panel(data,panel_params,coord,rule=rule))
    groups <- split(data,data$group)
    device_pixels <- if(grDevices::dev.cur()==1L) 1e6 else prod(grDevices::dev.size("px"))
    pixel_budget <- min(8e6,max(2e6,device_pixels*4))/max(1,length(groups))
    reduced <- 0L
    grobs <- lapply(groups,function(d) {
      if(!is.null(d$flow_tip) && !is.na(d$flow_tip[1]))tip <- as.character(d$flow_tip[1])
      nq <- as.integer(d$flow_qn[1]); ns <- as.integer(d$flow_sn[1])
      ne <- (nrow(d)-nq-ns)/2
      if (nq < 1 || ns < 1 || ne < 2 || ne != floor(ne))
        return(ggplot2::GeomPolygon$draw_panel(d,panel_params,coord,rule=rule))
      # Geometry is query front, second edge, reversed subject front,
      # reversed first edge. Both edges below run query -> subject.
      edge2 <- d[nq+seq_len(ne),c("x","y")]
      edge1 <- d[nrow(d)-seq_len(ne)+1L,c("x","y")]
      xy1 <- as.matrix(edge1); xy2 <- as.matrix(edge2)
      t <- seq(0,1,length.out=ne)
      subject_rows <- nq+ne+seq_len(ns)
      if (tip == "subject") {
        centre <- (xy1+xy2)/2
        span <- xy2[ne,]-xy1[ne,]
        if (sum(span^2)<1e-20 || sum((centre[ne,]-centre[1,])^2)<1e-20) {
          warning("Ribbon subject tip has a degenerate interval/channel; keeping the full front",call.=FALSE)
          tip <- "none"
        } else {
          shoulder <- max(2L,floor(.88*(ne-1L))+1L)
          front <- as.matrix(d[subject_rows,c("x","y")])
          middle <- (nrow(front)+1)/2
          marker <- colMeans(front[unique(c(floor(middle),ceiling(middle))),,drop=FALSE])
          base1 <- centre[shoulder,]-span/2
          base2 <- centre[shoulder,]+span/2
          u <- seq(0,1,length.out=ne-shoulder+1L)
          xy1[shoulder:ne,] <- outer(1-u,base1)+outer(u,marker)
          xy2[shoulder:ne,] <- outer(1-u,base2)+outer(u,marker)
        }
      }
      colors <- if (flow) {
        if (is.na(d$flow_from[1]) || is.na(d$flow_to[1])) rep(NA_character_,ne-1L) else
          { endpoints <- c(d$flow_from[1],d$flow_to[1])
            opacity <- grDevices::col2rgb(endpoints,alpha=TRUE)[4,]/255
            scales::alpha(grDevices::colorRampPalette(endpoints,space="Lab")(ne-1L),
              seq(opacity[1],opacity[2],length.out=ne-1L)) }
      } else rep(d$fill[1],ne-1L)
      subject_rows <- nq+ne+seq_len(ns)
      xy <- rbind(as.matrix(d[seq_len(nq),c("x","y")]),xy2,
        if(tip=="none") as.matrix(d[subject_rows,c("x","y")]),xy1[ne:1,,drop=FALSE])
      fill <- ggchord_ribbon_mesh_fill(xy1,xy2,xy,colors,d$alpha[1],coord,panel_params,pixel_budget)
      if(isTRUE(attr(fill,"ggchord_reduced_fill")))reduced <<- reduced+1L
      border <- d[rep(1L,nrow(xy)),,drop=FALSE]
      border$x<-xy[,1];border$y<-xy[,2];border$fill<-NA
      grid::grobTree(fill,
        ggplot2::GeomPolygon$draw_panel(border,panel_params,coord,rule=rule))
    })
    if(reduced>0L)warning(sprintf("Gradient fill resolution reduced for %d ribbon(s) to bound device memory; source intervals and vector outlines are retained",reduced),call.=FALSE)
    do.call(grid::grobTree,grobs)
  })

ggchord_ribbon_flow <- function(layer, geometry, layout, plot, sequence_colours=NULL) {
  if (!identical(layer$ggchord_type,"ribbon") || !nrow(geometry)) return(list(layer=layer,data=geometry))
  flow <- ggchord_ribbon_uses_flow(layer,plot)
  if (flow) {
    if (is.null(sequence_colours)) sequence_colours <- ggchord_sequence_display_colours(plot,layout)
    endpoint_colour <- function(query) {
      ids <- if(is.null(layout$display_instances)) {
        if(query) geometry$qaccver else geometry$saccver
      } else if(query) geometry$q_instance else geometry$s_instance
      fallback <- unname((layout$legend_seq_colors %||% layout$seq_colors)[if(query) geometry$qaccver else geometry$saccver])
      for (rows in split(seq_len(nrow(geometry)),geometry$group)) {
        d<-geometry[rows,,drop=FALSE]
        front <- if(query) seq_len(d$.q_n[1]) else d$.q_n[1]+(nrow(d)-d$.q_n[1]-d$.s_n[1])/2+seq_len(d$.s_n[1])
        centre <- colMeans(d[front,c("x","y"),drop=FALSE])
        candidates <- which(sequence_colours$instance==ids[rows[1]])
        if(length(candidates)) {
          at <- candidates[which.min((sequence_colours$x[candidates]-centre[1])^2+(sequence_colours$y[candidates]-centre[2])^2)]
          fallback[rows] <- sequence_colours$colour[at]
        }
      }
      fallback
    }
    geometry$.flow_from <- endpoint_colour(TRUE)
    geometry$.flow_to <- endpoint_colour(FALSE)
    layer$mapping$ribbon_fill <- NULL
    layer$mapping$flow_from <- ggplot2::aes(flow_from=I(.flow_from))$flow_from
    layer$mapping$flow_to <- ggplot2::aes(flow_to=I(.flow_to))$flow_to
  }
  geometry$.flow_tip <- layer$ggchord_params$tip %||% "none"
  layer$mapping$flow_tip <- ggplot2::aes(flow_tip=I(.flow_tip))$flow_tip
  layer$mapping$flow_qn <- ggplot2::aes(flow_qn=I(.q_n))$flow_qn
  layer$mapping$flow_sn <- ggplot2::aes(flow_sn=I(.s_n))$flow_sn
  layer$geom_params$flow <- flow
  layer$geom_params$tip <- layer$ggchord_params$tip %||% "none"
  list(layer=layer,data=geometry)
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
  triangle<-function(p,v) {
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
    triangle(rbind(a[i,],a[i+1,],b[i+1,]),c(t[i],t[i+1],t[i+1]))
    triangle(rbind(a[i,],b[i+1,],b[i,]),c(t[i],t[i+1],t[i]))
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

# Resolve the actual sequence aesthetic through the same geom, staged mapping
# and scales as the visible layer. The small bare ggplot cannot recurse into
# ggchord layout. Last sequence layer is the visible owner where layers overlap.
ggchord_sequence_display_colours <- function(plot,layout) {
  layers <- Filter(function(layer) identical(layer$ggchord_type,"seq"),plot$layers)
  empty <- data.frame(instance=character(),x=numeric(),y=numeric(),colour=character())
  if(!length(layers)) layers <- list(geom_seq(show.legend=FALSE))
  layer <- tail(layers,1)[[1]]
  geometry <- extract_ggchord_layer_data(layer,layout)
  if(!nrow(geometry))return(empty)
  if (!".component" %in% names(geometry)) geometry$.component <- "path"
  if (!"group" %in% names(geometry)) geometry$group <- match(geometry$accver,unique(geometry$accver))
  layer$ggchord_resolved_input <- geometry
  layer <- reconstruct_layer(layer,geometry)
  native <- ggplot2::ggplot(mapping=plot$mapping)+layer
  explicit <- plot$scales$get_scales("seq_colour")
  if(!is.null(explicit))native$scales$add(explicit$clone())
  defaults <- list(scale_seq_colour_manual(values=layout$legend_seq_colors %||% layout$seq_colors,guide="none"))
  defaults <- ggchord_infer_visual_scales(native,layout,defaults)
  native <- attach_ggchord_scales(native,defaults)
  built <- ggplot2::ggplot_build(native)$data[[1]]
  data.frame(instance=if(is.null(layout$display_instances)) as.character(geometry$accver) else geometry$display_instance_id,
    x=geometry$x,y=geometry$y,colour=built$seq_colour,stringsAsFactors=FALSE)
}

ggchord_ribbon_uses_flow <- function(layer,plot) {
  if(!identical(layer$ggchord_type,"ribbon"))return(FALSE)
  user <- layer$ggchord_input_mapping
  if(isTRUE(layer$inherit.aes))user <- c(user,plot$mapping[setdiff(names(plot$mapping),names(user))])
  is.null(layer$aes_params$ribbon_fill) && !"ribbon_fill" %in% names(user) &&
    is.null(plot$scales$get_scales("ribbon_fill"))
}
