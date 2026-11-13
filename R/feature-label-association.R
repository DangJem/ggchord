# Fit the glyph envelopes actually drawn, rather than a whole-word width alone.
ggchord_feature_contains_label <- function(boxes, polygons, source_row) {
  own <- Filter(function(p) identical(p$source_row[1L], source_row), polygons)
  if (!length(own) || !nrow(boxes)) return(FALSE)
  all(vapply(seq_len(nrow(boxes)), function(i) {
    corners <- ggchord_box_corners(boxes[i, , drop = FALSE])
    probes <- rbind(corners, (corners + corners[c(2,3,4,1), ])/2,
      c(boxes$cx[i], boxes$cy[i]))
    any(vapply(own, function(p) all(vapply(seq_len(nrow(probes)), function(j)
      ggchord_point_in_polygon(probes[j,1], probes[j,2], p), logical(1))), logical(1)))
  }, logical(1)))
}

# Intersect the unchanged source-angle ray with its actual drawn contour.
# The side is a display attachment, never a new source position or interval.
ggchord_feature_label_attachment <- function(candidate, polygons, outward, units_per_inch) {
  if (isTRUE(candidate$repeat_count > 1L)) return(candidate)
  radius <- sqrt(candidate$anchor_x^2 + candidate$anchor_y^2)
  if (!is.finite(radius) || radius < 1e-8) return(candidate)
  direction <- c(candidate$anchor_x, candidate$anchor_y)/radius
  cross <- function(x,y) x[1]*y[2]-x[2]*y[1]
  distances <- numeric()
  for (p in polygons) {
    if (!identical(p$source_row[1L], candidate$source_row)) next
    xy <- cbind(p$x, p$y); xy <- rbind(xy, xy[1L, ])
    a <- xy[-nrow(xy), , drop = FALSE]
    edge <- xy[-1L, , drop = FALSE]-a
    den <- direction[1]*edge[,2]-direction[2]*edge[,1]
    at <- which(abs(den) > 1e-12)
    t <- (a[at,1]*edge[at,2]-a[at,2]*edge[at,1])/den[at]
    u <- (a[at,1]*direction[2]-a[at,2]*direction[1])/den[at]
    distances <- c(distances, t[u >= -1e-8 & u <= 1+1e-8 & t > 0])
  }
  if (!length(distances)) return(candidate)
  # A tiny visible-side clearance avoids a numerical re-entry at the boundary.
  edge <- if (outward) max(distances) else min(distances)
  edge <- edge + if (outward) .03*units_per_inch/25.4 else -.03*units_per_inch/25.4
  candidate$anchor_x <- edge*direction[1]
  candidate$anchor_y <- edge*direction[2]
  candidate
}

ggchord_feature_label_shifted <- function(candidate, height) {
  anchor_angle <- atan2(candidate$anchor_y, candidate$anchor_x)
  angle <- atan2(candidate$text_y, candidate$text_x)-anchor_angle
  abs(atan2(sin(angle), cos(angle))) * sqrt(candidate$text_x^2+candidate$text_y^2) > height*.5
}
# Preserve straight source-to-label association. Clip glyph envelopes instead
# of making a visibility-graph detour; transparent feature fills stay honest.
ggchord_simple_feature_leaders <- function(registry, layout) {
  registry <- ggchord_clear_simple_feature_exits(registry, layout)
  boxes <- ggchord_annotation_text_obstacles(registry, layout)
  units <- layout$text_units_per_inch %||% .35
  source_maps <- ggchord_annotation_source_maps(layout)
  attachment_key <- function(owner, source, components) {
    keys <- unique(source_maps$source_id[source_maps$layer_id %in% owner &
      source_maps$input_row %in% source & source_maps$component %in% components])
    if (length(keys)==1L) keys else paste(owner,source,sep="/")
  }
  feature_registry <- registry
  for (id in names(feature_registry)) {
    p <- feature_registry[[id]]$gene_poly
    if (!is.data.frame(p) || !nrow(p)) next
    p$.attachment_id <- vapply(p$source_row,function(source)
      attachment_key(id,source,"gene_poly"),character(1))
    feature_registry[[id]]$gene_poly <- p
  }
  for (owner in names(registry)) {
    geometry <- registry[[owner]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry) ||
        !"label_layout" %in% names(geometry) ||
        !any(geometry$label_layout %in% "feature")) next
    segments <- which(geometry$.component %in% "segment")
    if (!length(segments)) next
    parts <- list(geometry[-segments, , drop = FALSE])
    texts <- geometry[geometry$.component %in% "text", , drop = FALSE]
    for (source in unique(geometry$source_row[segments])) {
      row <- texts[texts$source_row %in% source, , drop = FALSE]
      if (nrow(row) != 1L || is.na(row$label) || !nzchar(row$label)) next
      if (row$feature_label_mode %in% "external") {
        # The shared perimeter has already ordered and attached this route.
        # Reconnecting its endpoints can undo that ordering and cross callouts.
        own <- boxes$owner == owner & !is.na(boxes$source_row) & boxes$source_row==source
        parts[[length(parts)+1L]] <- ggchord_feature_route_occlusion(
          geometry[segments[geometry$source_row[segments]==source],,drop=FALSE],
          row,boxes[!own,,drop=FALSE])
        next
      }
      own <- boxes$owner == owner & !is.na(boxes$source_row) & boxes$source_row == source
      row$.attachment_id <- attachment_key(owner,source,
        c("gene_text","gene_text_repel","gene_label_repel"))
      if (ggchord_feature_short_attachment(row,boxes[own,,drop=FALSE],
          feature_registry,units)) next
      route <- data.frame(x0 = row$anchor_x, y0 = row$anchor_y,
        x1 = row$x, y1 = row$y)
      own_cuts <- ggchord_inner_leader_intervals(route, boxes[own, , drop = FALSE])
      if (nrow(own_cuts)) {
        stop <- max(0, min(own_cuts$start) - .01*units /
          max(sqrt((route$x1-route$x0)^2+(route$y1-route$y0)^2), 1e-8))
        route$x1 <- route$x0 + stop*(route$x1-route$x0)
        route$y1 <- route$y0 + stop*(route$y1-route$y0)
      }
      policy <- row$.leader_overlap %||% "fade"
      fade <- row$.leader_overlap_alpha %||% .18
      if ((route$x1-route$x0)^2+(route$y1-route$y0)^2 < 1e-16) next
      source_segments <- segments[geometry$source_row[segments] == source]
      base_alpha <- geometry$alpha[source_segments]
      base_alpha <- if (any(is.finite(base_alpha))) max(base_alpha, na.rm=TRUE) else NA_real_
      cuts <- if (identical(policy, "show")) data.frame(start=numeric(),end=numeric()) else
        ggchord_inner_leader_intervals(route, boxes[!own, , drop = FALSE])
      times <- sort(unique(c(0, 1, cuts$start, cuts$end)))
      times <- times[times >= 0 & times <= 1]
      for (i in seq_len(length(times)-1L)) {
        a <- times[i]; b <- times[i+1L]; middle <- (a+b)/2
        covered <- any(middle > cuts$start & middle < cuts$end)
        if (b-a < 1e-10 || (covered && identical(policy, "clip"))) next
        d <- geometry[segments[which(geometry$source_row[segments] == source)[1L]], , drop = FALSE]
        d$x <- route$x0 + a*(route$x1-route$x0)
        d$y <- route$y0 + a*(route$y1-route$y0)
        d$xend <- route$x0 + b*(route$x1-route$x0)
        d$yend <- route$y0 + b*(route$y1-route$y0)
        for (pair in list(c('x0','x'), c('y0','y'), c('x1','xend'), c('y1','yend')))
          if (pair[1L] %in% names(d)) d[[pair[1L]]] <- d[[pair[2L]]]
        d$occluded <- covered
        d$alpha <- if (covered && identical(policy, "fade")) base_alpha * fade else base_alpha
        d$leader_text_status <- 'direct_single_segment'
        parts[[length(parts)+1L]] <- d
      }
    }
    registry[[owner]]$gene_label_repel <- ggchord_rbind_fill(parts)
  }
  registry
}

# Omit an isolated, centred, near-source attachment. Tangential displacement,
# intervening features, and competing nearby sources still require a leader.
ggchord_feature_short_attachment <- function(row, boxes, registry, units) {
  if (!row$feature_label_mode %in% "adjacent" ||
      isTRUE(row$repeat_count>1L) || !nrow(boxes)) return(FALSE)
  if (ggchord_feature_label_shifted(row,max(boxes$oh))) return(FALSE)
  route <- data.frame(x0=row$anchor_x,y0=row$anchor_y,x1=row$x,y1=row$y)
  cuts <- ggchord_inner_leader_intervals(route,boxes)
  fraction <- if(nrow(cuts)) min(cuts$start) else 1
  distance <- sqrt((row$x-row$anchor_x)^2+(row$y-row$anchor_y)^2)
  if (distance*fraction > 1.5*units/25.4) return(FALSE)
  # Association remains explicit when another feature is just as close to
  # the name. Compare contour distance, not distant polygon vertices alone.
  for (owner in registry) {
    p <- owner$gene_poly
    if (!is.data.frame(p) || !nrow(p)) next
    p <- p[p$.component %in% "polygon",,drop=FALSE]
    for (part in split(p,p$group,drop=TRUE)) {
      if (!nrow(part)) next
      # Shared feature/label source IDs are local; caller supplies own key.
      if (identical(as.character(part$.attachment_id[1L]),row$.attachment_id)) next
      xy<-cbind(part$x,part$y);edge<-xy[c(2:nrow(xy),1L),,drop=FALSE]-xy
      relative<-matrix(c(row$x,row$y),nrow(xy),2,byrow=TRUE)-xy
      along<-pmin(1,pmax(0,rowSums(relative*edge)/pmax(rowSums(edge^2),1e-20)))
      closest<-xy+edge*along
      if (min(sqrt((closest[,1L]-row$x)^2+(closest[,2L]-row$y)^2)) <=
          distance*1.15) return(FALSE)
    }
  }
  TRUE
}

# Preserve an already ordered external route and apply the requested text
# overlap appearance to each original leg. Splitting a leg adds no bends.
ggchord_feature_route_occlusion <- function(rows, text, boxes) {
  policy <- as.character(text$.leader_overlap %||% "fade")[1L]
  if (!nrow(rows) || !nrow(boxes) || identical(policy,"show")) return(rows)
  fade <- as.numeric(text$.leader_overlap_alpha %||% .18)[1L]
  nominal <- if(any(is.finite(rows$alpha))) max(rows$alpha,na.rm=TRUE) else 1
  parts <- list()
  for(i in seq_len(nrow(rows))) {
    row <- rows[i,,drop=FALSE]
    route <- data.frame(x0=row$x,y0=row$y,x1=row$xend,y1=row$yend)
    cuts <- ggchord_inner_leader_intervals(route,boxes)
    times <- sort(unique(c(0,1,cuts$start,cuts$end)))
    for(k in seq_len(length(times)-1L)) {
      a<-times[k];b<-times[k+1L];middle<-(a+b)/2
      covered<-any(middle>cuts$start & middle<cuts$end)
      if(b-a<1e-10 || (covered && identical(policy,"clip")))next
      part<-row
      part$x<-route$x0+a*(route$x1-route$x0)
      part$y<-route$y0+a*(route$y1-route$y0)
      part$xend<-route$x0+b*(route$x1-route$x0)
      part$yend<-route$y0+b*(route$y1-route$y0)
      for(pair in list(c('x0','x'),c('y0','y'),c('x1','xend'),c('y1','yend')))
        if(pair[1L] %in% names(part))part[[pair[1L]]]<-part[[pair[2L]]]
      part$occluded<-covered
      part$alpha<-nominal*if(covered && identical(policy,"fade"))fade else 1
      parts[[length(parts)+1L]]<-part
    }
  }
  ggchord_rbind_fill(parts) %||% rows[0,,drop=FALSE]
}
