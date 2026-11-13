# Conservative polygon interiors, not exact ink/occlusion. Endpoints may touch.
association_segment_inside <- function(a, b, polygon) {
  if (nrow(polygon) < 3L) return(FALSE)
  if (max(a[1], b[1]) < min(polygon$x) || min(a[1], b[1]) > max(polygon$x) ||
      max(a[2], b[2]) < min(polygon$y) || min(a[2], b[2]) > max(polygon$y)) return(FALSE)
  cross <- function(x, y) x[1]*y[2]-x[2]*y[1]
  delta <- b-a; times <- c(0, 1)
  xy <- cbind(polygon$x, polygon$y); xy <- rbind(xy, xy[1L, ])
  for (j in seq_len(nrow(xy)-1L)) {
    edge <- xy[j+1L, ]-xy[j, ]; den <- cross(delta, edge)
    if (abs(den) < 1e-12) next
    t <- cross(xy[j, ]-a, edge)/den
    u <- cross(xy[j, ]-a, delta)/den
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) times <- c(times, t)
  }
  times <- sort(unique(times)); middle <- (head(times, -1L)+tail(times, -1L))/2
  any(vapply(middle, function(t) ggchord_point_in_polygon(
    a[1]+t*delta[1], a[2]+t*delta[2], polygon), logical(1)))
}

association_metrics <- function(layout) {
  maps <- ggchord_annotation_source_maps(layout)
  identity <- function(owner, row, components) {
    values <- unique(maps$source_id[maps$layer_id == owner & maps$input_row == row &
      maps$component %in% components])
    if (length(values) == 1L) values else NA_character_
  }
  polygons <- list(); labels <- list(); crossings <- list()
  for (owner in names(layout$layer_geometry)) {
    d <- layout$layer_geometry[[owner]]$gene_poly
    if (!is.data.frame(d) || !nrow(d)) next
    d <- d[is.na(d$.component) | d$.component == 'polygon', , drop = FALSE]
    for (group in unique(d$group)) {
      p <- d[d$group == group, , drop = FALSE]
      p$association_id <- identity(owner, p$source_row[1L], 'gene_poly')
      p$polygon_id <- paste(owner, group, sep = '/')
      polygons[[length(polygons)+1L]] <- p
    }
  }
  for (owner in names(layout$layer_geometry)) {
    d <- layout$layer_geometry[[owner]]$gene_label_repel
    if (!is.data.frame(d) || !nrow(d)) next
    texts <- d[d$.component == 'text' & !is.na(d$label), , drop = FALSE]
    for (i in seq_len(nrow(texts))) {
      t <- texts[i, ]; source <- t$source_row
      source_id <- identity(owner, source, c('gene_text', 'gene_text_repel', 'gene_label_repel'))
      own <- Filter(function(p) !is.na(source_id) && identical(p$association_id[1L], source_id), polygons)
      glyphs <- d[d$.component == 'arc_text' & d$source_row == source, , drop = FALSE]
      if (!nrow(glyphs)) glyphs <- t
      glyphs$text <- glyphs$label; glyphs$text_x <- glyphs$x; glyphs$text_y <- glyphs$y
      glyphs$text_angle <- glyphs$angle
      boxes <- ggchord_text_boxes(glyphs, units_per_inch = layout$text_units_per_inch, box_padding = 0)
      corners <- ggchord_boxes_corners(boxes)
      contained <- length(own) > 0 && all(vapply(seq_len(nrow(corners)), function(j)
        any(vapply(own, function(p) ggchord_point_in_polygon(corners[j,1], corners[j,2], p), logical(1))), logical(1)))
      routes <- d[d$.component == 'segment' & d$source_row == source, , drop = FALSE]
      shift <- abs(atan2(sin(atan2(t$y,t$x)-atan2(t$anchor_y,t$anchor_x)),
        cos(atan2(t$y,t$x)-atan2(t$anchor_y,t$anchor_x)))) * sqrt(t$x^2+t$y^2)
      tangent_error <- if (nrow(glyphs) && isTRUE(t$.draw_as_arc)) max(abs(
        ((glyphs$angle-atan2(glyphs$y,glyphs$x)*180/pi-90+90)%%180)-90)) else NA_real_
      # Count supporting lines, merging collinear fragments caused by glyph
      # clipping. This observes route complexity without locking coordinates.
      supports <- list()
      for (j in seq_len(nrow(routes))) {
        a <- c(routes$x[j],routes$y[j]); v <- c(routes$xend[j],routes$yend[j])-a
        length <- sqrt(sum(v^2))
        if (length < 1e-10) next
        v <- v/length
        same <- any(vapply(supports, function(q)
          abs(v[1]*q$v[2]-v[2]*q$v[1]) < 1e-6 &&
          abs((a[1]-q$a[1])*q$v[2]-(a[2]-q$a[2])*q$v[1]) < 1e-7, logical(1)))
        if (!same) supports[[length(supports)+1L]] <- list(a=a,v=v)
      }
      # Check the actual drawable owner's contour, not the label sub-layout's
      # recomputed surrogate. A small physical tolerance admits edge clearance.
      shared <- "repeat_count" %in% names(t) && isTRUE(t$repeat_count>1L)
      anchor_on_source <- if (shared || !length(own)) NA else any(vapply(own, function(p) {
        a <- cbind(p$x,p$y);b<-a[c(2:nrow(a),1L),,drop=FALSE]
        edge<-b-a;point<-c(t$anchor_x,t$anchor_y)
        delta<-matrix(point,nrow(a),2,byrow=TRUE)-a
        fraction<-pmin(1,pmax(0,rowSums(delta*edge)/pmax(rowSums(edge^2),1e-20)))
        nearest<-a+edge*fraction
        ggchord_point_in_polygon(point[1L],point[2L],p) ||
          min(sqrt(rowSums((nearest-matrix(point,nrow(a),2,byrow=TRUE))^2))) <=
            .1*(layout$text_units_per_inch %||% .35)/25.4
      },logical(1)))
      labels[[length(labels)+1L]] <- data.frame(owner = owner, source_row = source,
        source_id = source_id, label = t$label, mode = t$feature_label_mode,
        contained = contained, anchor_on_source = anchor_on_source,
        tangent_error_deg = tangent_error,
        arc_text = isTRUE(t$.draw_as_arc),
        supporting_segments = length(supports), tangential_shift = shift, has_leader = nrow(routes) > 0L,
        ambiguous_adjacent = t$feature_label_mode == 'adjacent' && !nrow(routes) && shift > max(boxes$oh)*.5)
      for (j in seq_len(nrow(routes))) for (p in polygons) if (association_segment_inside(
          c(routes$x[j], routes$y[j]), c(routes$xend[j], routes$yend[j]), p)) {
        crossings[[length(crossings)+1L]] <- data.frame(owner = owner, source_row = source,
          label = t$label, polygon_id = p$polygon_id[1L], feature_source_id = p$association_id[1L],
          own = !is.na(source_id) && identical(source_id, p$association_id[1L]))
      }
    }
  }
  labels <- ggchord_rbind_fill(labels) %||% data.frame()
  crossings <- unique(ggchord_rbind_fill(crossings) %||% data.frame())
  list(labels = labels, leader_feature_crossings = crossings,
    summary = list(inside_labels_outside_own_feature = sum(labels$mode == 'inside' & !labels$contained),
      anchors_off_own_feature = sum(labels$anchor_on_source %in% FALSE),
      internal_glyphs_off_tangent = sum(labels$tangent_error_deg>1,na.rm=TRUE),
      ambiguous_adjacent_without_leader = sum(labels$ambiguous_adjacent),
      leaders_over_one_bend = sum(labels$supporting_segments > 2L),
      internal_labels_without_arc = sum(labels$mode %in% c('inside','adjacent') & !labels$arc_text),
      leader_feature_crossings = nrow(crossings)))
}

association_text_quality <- function(table) {
  layout <- attr(table,'ggchord_render_layout')
  context <- attr(table,'ggchord_render_context')
  built <- attr(table,'ggchord_render_built')
  text <- ggchord_rbind_fill(lapply(seq_along(built$data),function(i) {
    d <- built$data[[i]]; layer <- built$plot$layers[[i]]
    if (is.data.frame(layer$data) && nrow(layer$data)==nrow(d) &&
        'source_row' %in% names(layer$data)) d$source_row <- layer$data$source_row
    t <- ggchord_diagnostic_text(ggchord_diagnostic_text_styles(d,layer),
      built$plot$coordinates,built$layout$panel_params[[1L]],context,
      parse=layer$geom_params$parse %||% FALSE,size_unit=layer$geom_params$size.unit %||% 'mm')
    if (!nrow(t)) return(NULL)
    t$layer_id <- layer$ggchord_layer_id %||% paste0('native-',i)
    t
  })) %||% data.frame()
  ggchord_label_quality(layout,text,built,context)
}
