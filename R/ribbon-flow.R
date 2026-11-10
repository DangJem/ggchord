# Ribbon colour follows the sampled Bezier parameter, never screen x/y.
GeomRibbonFlow <- ggplot2::ggproto("GeomRibbonFlow", ggplot2::GeomPolygon,
  default_aes = utils::modifyList(ggplot2::GeomPolygon$default_aes,
    ggplot2::aes(flow_from=NA_character_, flow_to=NA_character_, flow_qn=0, flow_sn=0)),
  extra_params = c("na.rm", "flow", "tip"),
  draw_panel = function(data, panel_params, coord, na.rm=FALSE, rule="evenodd",
      flow=FALSE, tip="none") {
    if (!flow && tip == "none") return(ggplot2::GeomPolygon$draw_panel(data,panel_params,coord,rule=rule))
    grobs <- lapply(split(data,data$group),function(d) {
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
      if (tip == "subject") {
        centre <- (xy1+xy2)/2
        taper <- pmin(1,(1-t)/.12)
        xy1 <- centre+(xy1-centre)*taper
        xy2 <- centre+(xy2-centre)*taper
      }
      colors <- if (flow && !is.na(d$flow_from[1]) && !is.na(d$flow_to[1]))
        grDevices::colorRampPalette(c(d$flow_from[1],d$flow_to[1]),space="Lab")(ne-1L) else
        rep(d$fill[1],ne-1L)
      strips <- lapply(seq_len(ne-1L),function(i) {
        piece <- d[rep(1L,4L),,drop=FALSE]
        xy <- rbind(xy1[i,],xy1[i+1L,],xy2[i+1L,],xy2[i,])
        piece$x<-xy[,1];piece$y<-xy[,2];piece$group<-i
        piece$fill<-colors[i];piece$colour<-colors[i];piece$linewidth<-.04
        piece
      })
      fill <- ggplot2::GeomPolygon$draw_panel(do.call(rbind,strips),panel_params,coord,rule=rule)
      border <- d[rep(1L,2L*ne),,drop=FALSE]
      xy <- rbind(xy1,xy2[ne:1,,drop=FALSE]);border$x<-xy[,1];border$y<-xy[,2];border$fill<-NA
      grid::grobTree(fill,ggplot2::GeomPolygon$draw_panel(border,panel_params,coord,rule=rule))
    })
    do.call(grid::grobTree,grobs)
  })

ggchord_ribbon_flow <- function(layer, geometry, layout, plot) {
  if (!identical(layer$ggchord_type,"ribbon") || !nrow(geometry)) return(list(layer=layer,data=geometry))
  user <- layer$ggchord_input_mapping
  if (isTRUE(layer$inherit.aes)) user <- c(user,plot$mapping[setdiff(names(plot$mapping),names(user))])
  flow <- is.null(layer$aes_params$ribbon_fill) && !"ribbon_fill" %in% names(user) &&
    is.null(plot$scales$get_scales("ribbon_fill"))
  if (flow) {
    colours <- layout$legend_seq_colors %||% layout$seq_colors
    seq_layer <- Filter(function(l) identical(l$ggchord_type,"seq"), plot$layers)
    fixed_colour <- if (length(seq_layer)) tail(seq_layer,1)[[1]]$aes_params$seq_colour else NULL
    seq_scale <- plot$scales$get_scales("seq_colour")
    if (!is.null(seq_scale)) {
      seq_scale <- seq_scale$clone();seq_scale$train(unique(c(geometry$qaccver,geometry$saccver)))
      colour_for <- function(ids) seq_scale$map(ids)
    } else colour_for <- function(ids) unname(colours[ids])
    if (!is.null(fixed_colour)) colour_for <- function(ids) rep_len(fixed_colour,length(ids))
    geometry$.flow_from <- colour_for(geometry$qaccver)
    geometry$.flow_to <- colour_for(geometry$saccver)
    layer$mapping$ribbon_fill <- NULL
    layer$mapping$flow_from <- ggplot2::aes(flow_from=I(.flow_from))$flow_from
    layer$mapping$flow_to <- ggplot2::aes(flow_to=I(.flow_to))$flow_to
  }
  layer$mapping$flow_qn <- ggplot2::aes(flow_qn=I(.q_n))$flow_qn
  layer$mapping$flow_sn <- ggplot2::aes(flow_sn=I(.s_n))$flow_sn
  layer$geom_params$flow <- flow
  layer$geom_params$tip <- layer$ggchord_params$tip %||% "subject"
  list(layer=layer,data=geometry)
}
