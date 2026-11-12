# Optional physical perimeter placement; explicit allocations never move.
ggchord_fusion_auto_rectangles <- function(x, rects, canvas) {
  ids <- names(x$panels)[vapply(x$panels, function(p) p$side == 'auto', logical(1))]
  if (!length(ids)) return(rects)
  fixed <- Filter(Negate(is.null), rects)
  center <- rects$center
  anchor <- c(mean(center[1:2]), mean(center[3:4]))*canvas
  gap <- 1; margin <- 2
  count <- length(ids); step <- 360/count
  # A minimum default container size; explicit fonts and sizes stay untouched.
  default <- max(15, min(.18*min(canvas), .7*2*sum(canvas)/count))
  fits <- function(a, boxes) {
    a <- a*rep(canvas, each = 2L)
    all(vapply(boxes, function(b) {
      b <- b*rep(canvas, each = 2L)
      min(a[2]+gap, b[2]) <= max(a[1]-gap, b[1]) ||
        min(a[4]+gap, b[4]) <= max(a[3]-gap, b[3])
    }, logical(1)))
  }
  phases <- unique(c(0, step*c(.25, -.25, .5, -.5, .75, -.75), seq(0, 360, length.out = 17L)[-17L]))
  progress <- seq(0, 360, length.out = max(64L, count*16L)+1L)
  progress <- progress[-length(progress)]
  for (factor in c(1, .9, .8, .7)) for (phase in phases) {
    trial <- rects; placed <- fixed; last <- -Inf; valid <- TRUE
    for (i in seq_along(ids)) {
      id <- ids[i]; explicit <- x$panels[[id]]$size
      size <- explicit %||% (rep(max(15, default*factor), 2L)/canvas)
      half <- size*canvas/2
      if (any(size*canvas+2*margin > canvas)) { valid <- FALSE; break }
      candidates <- progress[progress > last]
      candidates <- candidates[order(abs(candidates-(i-1L)*step), candidates)]
      selected <- NULL
      for (angle in candidates) {
        direction <- c(cos((90-phase-angle)*pi/180), sin((90-phase-angle)*pi/180))
        boundary <- ifelse(direction > 0, canvas-margin-half, margin+half)
        distance <- (boundary-anchor)/direction
        distance[abs(direction) < 1e-10] <- Inf
        radius <- min(distance)
        if (!is.finite(radius) || radius <= 0) next
        xy <- anchor+direction*radius
        box <- c(xy[1]-half[1], xy[1]+half[1], xy[2]-half[2], xy[2]+half[2])/rep(canvas, each = 2L)
        if (any(box < 0 | box > 1) || !fits(box, placed)) next
        attr(box, 'auto_angle') <- (90-phase-angle) %% 360
        attr(box, 'auto_size_factor') <- if (is.null(explicit)) max(15, default*factor)/default else 1
        selected <- box; last <- angle; break
      }
      if (is.null(selected)) { valid <- FALSE; break }
      trial[[id]] <- selected; placed[[id]] <- selected
    }
    if (!valid) next
    attr(trial, 'automatic_layout') <- list(status = 'feasible', panel_ids = ids,
      strategy = 'bounded_clockwise_perimeter', start_phase = phase,
      default_size_factor = min(vapply(trial[ids], function(b) attr(b, 'auto_size_factor'), numeric(1))), minimum_default_container_mm = 15, gap_mm = gap,
      fixed_allocations_preserved = TRUE, placement_rescales_fonts = FALSE)
    return(trial)
  }
  ggchord_stop('Fusion: automatic perimeter panels cannot fit without overlap under the supplied constraints; enlarge the canvas, reduce fixed sizes or provide explicit boxes')
}
