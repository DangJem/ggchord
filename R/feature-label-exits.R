# Try coordinated straight/one-elbow exits without moving text or source roots.
# A whole conflicting owner group is accepted together, or left unchanged.
ggchord_clear_simple_feature_exits <- function(registry, layout) {
  if (!isTRUE(layout$circular)) return(registry)
  boxes <- ggchord_annotation_text_obstacles(registry, layout)
  leaders <- ggchord_annotation_leader_segments(registry)
  if (!nrow(boxes) || !nrow(leaders)) return(registry)
  units <- layout$text_units_per_inch %||% .35
  key <- function(d) paste(d$owner, d$component, d$source, sep='\r')
  route <- function(d) data.frame(x0=d$x1, y0=d$y1, x1=d$x2, y1=d$y2)
  own_boxes <- function(d) boxes$owner %in% d$owner[1L] &
    boxes$source_row %in% d$source[1L]
  crosses <- function(d, others) nrow(others) && any(vapply(seq_len(nrow(d)),
    function(i) any(ggchord_route_segment_crossings(c(d$x1[i],d$y1[i]),
      c(d$x2[i],d$y2[i]),others)),logical(1)))
  for (owner in names(registry)) {
    geometry <- registry[[owner]]$gene_label_repel
    if (!is.data.frame(geometry) || !all(c('.component','feature_label_mode','label_layout') %in% names(geometry))) next
    labels <- geometry[geometry$.component %in% 'text' &
      geometry$feature_label_mode %in% 'external' & geometry$label_layout %in% 'feature',,drop=FALSE]
    eligible <- leaders$owner %in% owner & leaders$component %in% 'gene_label_repel' &
      leaders$source %in% labels$source_row
    keys <- key(leaders)
    blocked <- Filter(function(k) {
      d <- leaders[keys %in% k,,drop=FALSE]
      ggchord_inner_leader_hits(route(d),boxes[!own_boxes(d),,drop=FALSE])
    },unique(keys[eligible]))
    if (!length(blocked)) next
    accepted <- list()
    pending <- blocked
    trial <- leaders
    complete <- TRUE
    for (k in blocked) {
      d <- trial[key(trial) %in% k,,drop=FALSE]
      starts <- which(!vapply(seq_len(nrow(d)),function(i)
        any((d$x2-d$x1[i])^2+(d$y2-d$y1[i])^2 < 1e-16),logical(1)))
      ends <- which(!vapply(seq_len(nrow(d)),function(i)
        any((d$x1-d$x2[i])^2+(d$y1-d$y2[i])^2 < 1e-16),logical(1)))
      if (length(starts)!=1L || length(ends)!=1L) {complete<-FALSE;break}
      a <- c(d$x1[starts],d$y1[starts]); b <- c(d$x2[ends],d$y2[ends])
      radius <- sqrt(sum(a^2))
      if (!is.finite(radius) || radius<=0) {complete<-FALSE;break}
      others <- trial[!key(trial) %in% pending,,drop=FALSE]
      obstacles <- boxes[!own_boxes(d),,drop=FALSE]
      hit <- vapply(seq_len(nrow(obstacles)),function(j)
        ggchord_inner_leader_hits(route(d),obstacles[j,,drop=FALSE]),logical(1))
      corners <- obstacles[hit,,drop=FALSE]
      corners$ow <- corners$ow+.4*units/25.4
      corners$oh <- corners$oh+.4*units/25.4
      direction <- a/radius; tangent <- c(-direction[2L],direction[1L])
      grid <- expand.grid(radial=c(-.1,-.05,0,.025,.05,.075,.1,.15,.2,.3,.5),
        tangent=c(0,.025,-.025,.05,-.05,.1,-.1,.15,-.15,.25,-.25,.4,-.4,.6,-.6,.8,-.8))
      scale <- max(1,radius)
      points <- rbind(b,a,outer(grid$radial*scale,direction)+
        outer(grid$tangent*scale,tangent)+matrix(a,nrow(grid),2,byrow=TRUE),
        ggchord_boxes_corners(corners))
      cost <- sqrt(rowSums((points-matrix(a,nrow(points),2,byrow=TRUE))^2))+
        sqrt(rowSums((points-matrix(b,nrow(points),2,byrow=TRUE))^2))
      limit <- sum(sqrt((d$x2-d$x1)^2+(d$y2-d$y1)^2))*1.3+2*units/25.4
      chosen <- NULL
      for (i in order(cost)) {
        if (!is.finite(cost[i]) || cost[i]>limit) next
        elbow <- points[i,]
        candidate <- d[rep(1L,2L),,drop=FALSE]
        candidate$x1 <- c(a[1],elbow[1]);candidate$y1 <- c(a[2],elbow[2])
        candidate$x2 <- c(elbow[1],b[1]);candidate$y2 <- c(elbow[2],b[2])
        candidate <- candidate[(candidate$x2-candidate$x1)^2+
          (candidate$y2-candidate$y1)^2>1e-16,,drop=FALSE]
        if (!nrow(candidate) || ggchord_inner_leader_hits(route(candidate),obstacles) ||
            crosses(candidate,others)) next
        chosen <- candidate;break
      }
      if (is.null(chosen)) {complete<-FALSE;break}
      accepted[[k]] <- chosen
      trial <- rbind(trial[!key(trial) %in% k,,drop=FALSE],chosen)
      pending <- setdiff(pending,k)
    }
    if (!complete) next
    # Recheck the entire solution, including all previously pending routes.
    valid <- all(vapply(accepted,function(d)
      !ggchord_inner_leader_hits(route(d),boxes[!own_boxes(d),,drop=FALSE]) &&
      !crosses(d,trial[!key(trial) %in% key(d),,drop=FALSE]),logical(1)))
    if (!valid) next
    for (d in accepted) {
      source <- d$source[1L]
      at <- geometry$.component %in% 'segment' & geometry$source_row %in% source
      template <- geometry[which(at)[1L],,drop=FALSE]
      geometry <- geometry[!at,,drop=FALSE]
      replacement <- template[rep(1L,nrow(d)),,drop=FALSE]
      replacement$x <- d$x1;replacement$y <- d$y1
      replacement$xend <- d$x2;replacement$yend <- d$y2
      for (pair in list(c('x0','x'),c('y0','y'),c('x1','xend'),c('y1','yend')))
        if (pair[1L] %in% names(replacement)) replacement[[pair[1L]]] <- replacement[[pair[2L]]]
      replacement$leader_text_status <- 'simple_exit_clearance'
      geometry <- ggchord_rbind_fill(list(geometry,replacement))
    }
    registry[[owner]]$gene_label_repel <- geometry
    leaders <- trial
  }
  registry
}
