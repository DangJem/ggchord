# Display-only grouping. Biological rows and feature geometry stay independent.
ggchord_prepare_repeat_labels <- function(gl, polygons, units_per_inch) {
  if (".repeat_angles" %in% names(gl)) return(gl)
  if (!"label_repeat" %in% names(gl) || !any(gl$label_repeat %in% "bracket")) return(gl)
  gl$source_rows <- I(lapply(gl$source_row, function(x) as.integer(x)))
  gl$repeat_count <- 1L
  gl$.repeat_angles <- I(rep(list(numeric()), nrow(gl)))
  gl$.repeat_radius <- NA_real_
  eligible <- gl$label_repeat %in% "bracket" & !is.na(gl$text) & nzchar(gl$text)
  angle <- atan2(gl$anchor_y, gl$anchor_x)
  spans <- ((gl$interval_end-gl$interval_start) %% gl$sequence_length) /
    gl$sequence_length * 2*pi
  key <- paste(gl$text, gl$repeat_type, gl$repeat_direction,
    gl$feature_label_fill, gl$family, gl$fontface, sep="\r")
  for (band in split(seq_len(nrow(gl)), paste(gl$accver,gl$lane,sep="\r"))) {
    polar <- (angle[band]+2*pi) %% (2*pi)
    ord <- order(polar); band <- band[ord]; polar <- polar[ord]
    cut <- which.max(diff(c(polar,polar[1L]+2*pi)))
    start <- if(cut==length(band)) 1L else cut+1L
    order <- c(seq.int(start,length(band)),if(start>1L) seq_len(start-1L) else integer())
    band <- band[order]; polar <- polar[order]
    polar <- polar+cumsum(c(0,as.numeric(diff(polar)<0)))*2*pi
    runs <- list(); current <- integer()
    for(k in seq_along(band)) {
      i <- band[k]
      connected <- length(current)>0L && eligible[i] && eligible[tail(band[current],1)] &&
        key[i]==key[tail(band[current],1)] && spans[i]>0 && spans[i]<.25 &&
        polar[k]-polar[tail(current,1)] >= (spans[i]+spans[tail(band[current],1)])/2 &&
        polar[k]-polar[tail(current,1)] <= 3*max(spans[i],spans[tail(band[current],1)],.02) &&
        polar[k]-polar[current[1L]] < pi/2
      if (!connected && length(current)) {runs[[length(runs)+1L]]<-current;current<-integer()}
      current<-c(current,k)
    }
    if(length(current)) runs[[length(runs)+1L]]<-current
    for(run in runs) {
      members <- band[run]
      if(length(members)<2L || !all(eligible[members]) || any(spans[members] >= .25)) next
      representative <- members[1L]
      radius <- vapply(members,function(i) {
        own<-polygons[polygons$source_row %in% gl$source_row[i] &
          polygons$.component %in% "polygon",,drop=FALSE]
        if(nrow(own)) min(sqrt(own$x^2+own$y^2)) else
          sqrt(gl$anchor_x[i]^2+gl$anchor_y[i]^2)
      },numeric(1))
      bracket_radius <- min(radius)-.3*units_per_inch/25.4
      extent <- range(c(polar[run]-spans[members]/2,polar[run]+spans[members]/2))
      centre <- mean(extent)
      gl$source_rows[[representative]] <- gl$source_row[members]
      gl$repeat_count[representative] <- length(members)
      gl$.repeat_angles[[representative]] <- extent
      gl$.repeat_radius[representative] <- bracket_radius
      gl$anchor_x[representative] <- bracket_radius*cos(centre)
      gl$anchor_y[representative] <- bracket_radius*sin(centre)
      gl$text[representative] <- paste0(gl$text[representative]," \u00d7 ",length(members))
      gl$feature_label_inside[members] <- FALSE
      gl$text[members[-1L]] <- NA_character_
      gl$repeat_count[members[-1L]] <- 0L
    }
  }
  gl
}

ggchord_repeat_bracket_geometry <- function(text) {
  if (!nrow(text) || !"repeat_count" %in% names(text)) return(data.frame())
  rows <- which(text$repeat_count>1L & !is.na(text$label) & nzchar(text$label))
  paths <- list()
  for(i in rows) {
    angle <- text$.repeat_angles[[i]]; radius <- text$.repeat_radius[i]
    if(length(angle)!=2L || !is.finite(radius)) next
    # A curved range bracket, distinct from the label's simple leader.
    tick <- min(.012, radius*.02)
    arc <- seq(angle[1L],angle[2L],length.out=48L)
    angles <- c(angle[1L],arc,angle[2L])
    radii <- c(radius+tick,rep(radius,length(arc)),radius+tick)
    d <- text[rep(i,length(angles)),,drop=FALSE]
    d$x <- radii*cos(angles); d$y <- radii*sin(angles)
    d$xend <- NA_real_;d$yend <- NA_real_
    d$label <- NA_character_;d$.component <- "repeat_bracket"
    d$group <- -abs(text$group[i]);d$alpha<-1;d$linetype<-"solid"
    paths[[length(paths)+1L]]<-d
  }
  ggchord_rbind_fill(paths) %||% data.frame()
}
