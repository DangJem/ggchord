# Internal provider protocol. Positive normal points toward increasing track
# radius in map_to_curve(), not necessarily away from the plot origin.
ggchord_empty_obstacles <- function() data.frame(accver = character(), start = numeric(),
  end = numeric(), side = character(), normal_min = numeric(), normal_max = numeric(),
  priority = numeric(), source_layer = character())

ggchord_entity_obstacles <- function(layer, layout, input) {
  polys <- layout$gene_polys
  if (is.null(input) || !nrow(input) || !nrow(polys)) return(ggchord_empty_obstacles())
  refs <- layout$sequence_reference
  parts <- lapply(split(polys, polys$group), function(d) {
    row <- input[d$source_row[1], , drop = FALSE]
    id <- as.character(row$accver)
    ref <- refs$refs[[id]]
    # Project the actual rendered polygon, including lollipop heads and stems.
    xy <- ggchord_rotate_points(as.matrix(d[c("x", "y")]), -layout$rotation)
    angles <- ref$angles
    base <- map_to_curve_many(angles, refs$radius[id], ref)
    outward <- map_to_curve_many(angles, refs$radius[id] + 1, ref) - base
    nearest <- vapply(seq_len(nrow(xy)), function(i) which.min(rowSums((base - matrix(xy[i,], nrow(base), 2, byrow = TRUE))^2)), integer(1))
    offset <- rowSums((xy - base[nearest,,drop=FALSE]) * outward[nearest,,drop=FALSE])
    step <- max(sqrt(rowSums(diff(base)^2)))
    f <- (angles[nearest] - refs$starts[id]) / (refs$ends[id] - refs$starts[id])
    if (refs$orientation[id] != 1) f <- 1 - f
    pos <- f * refs$lens[id] + 1
    pad <- refs$lens[id] / (length(angles) - 1)
    lo <- min(offset) - step; hi <- max(offset) + step
    data.frame(accver = id, start = max(1, min(pos) - pad),
      end = min(refs$lens[id], max(pos) + pad),
      side = if (lo >= 0) "positive" else if (hi <= 0) "negative" else "both",
      normal_min = lo, normal_max = hi, priority = 1,
      source_layer = layer$ggchord_layer_id, stringsAsFactors = FALSE)
  })
  do.call(rbind, parts)
}

ggchord_collect_obstacles <- function(plot, cache) {
  providers <- which(vapply(plot$layers, function(l) is.function(l$ggchord_obstacle_provider), logical(1)))
  if (!length(providers)) return(ggchord_empty_obstacles())
  seq_layers <- which(vapply(plot$layers, function(l) identical(l$ggchord_type, "seq"), logical(1)))
  parts <- lapply(providers, function(i) {
    layer <- plot$layers[[i]]
    sub <- plot
    sub$layers <- plot$layers[sort(unique(c(seq_layers, i)))]
    sub$ggchord$obstacles <- ggchord_empty_obstacles()
    layout <- compute_chord_geometry_single(sub, cache)
    input <- ggchord_resolve_layer_input(layer, plot$ggchord$data$gene_data)
    out <- layer$ggchord_obstacle_provider(layer, layout, input)
    ggchord_require_columns(out, names(ggchord_empty_obstacles()), "obstacle provider")
    numeric <- c("start", "end", "normal_min", "normal_max", "priority")
    if (any(!vapply(out[numeric], is.numeric, logical(1))) || any(!is.finite(as.matrix(out[numeric]))) ||
        any(out$start > out$end | out$normal_min > out$normal_max) ||
        any(!out$side %in% c("positive", "negative", "both"))) ggchord_stop("Invalid obstacle provider ranges")
    out
  })
  do.call(rbind, parts)
}

ggchord_gap_profile <- function(pos, id, obstacles, avoid = "smooth", close_gap = .035) {
  out <- rep(close_gap, length(pos))
  if (avoid == "none" || is.null(obstacles) || !nrow(obstacles) || !length(pos)) return(out)
  obs <- obstacles[obstacles$accver == id & obstacles$normal_max > 0 & obstacles$side != "negative",,drop=FALSE]
  if (!nrow(obs)) return(out)
  if (avoid == "uniform") {
    hit <- obs$start <= max(pos) & obs$end >= min(pos)
    if (any(hit)) out[] <- max(close_gap, obs$normal_max[hit] + .025)
    return(out)
  }
  transition <- max(diff(range(pos)) * .08, 1)
  for (i in seq_len(nrow(obs))) {
    distance <- pmax(obs$start[i] - pos, pos - obs$end[i], 0)
    weight <- ifelse(distance < transition, (1 + cos(pi * pmin(distance / transition, 1))) / 2, 0)
    out <- pmax(out, close_gap + pmax(0, obs$normal_max[i] + .025 - close_gap) * weight)
  }
  out
}

# Detect only new local front failures, not legitimate crossed alignments.
ggchord_front_invalid <- function(xy, baseline) {
  a <- diff(xy); b <- diff(baseline)
  if (any(rowSums(a * b) < -1e-10)) return(TRUE)
  n <- nrow(xy)
  if (n < 4) return(FALSE)
  cross <- function(a,b) a[1]*b[2]-a[2]*b[1]
  for (i in seq_len(n-3L)) for (j in seq.int(i+2L,n-1L)) {
    u <- xy[i+1,]-xy[i,]; v <- xy[j+1,]-xy[j,]
    den <- cross(u,v)
    if (abs(den) < 1e-12) next
    w <- xy[j,]-xy[i,]; t <- cross(w,v)/den; s <- cross(w,u)/den
    if (t > 1e-8 && t < 1-1e-8 && s > 1e-8 && s < 1-1e-8) return(TRUE)
  }
  FALSE
}
