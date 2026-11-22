# Internal provider protocol. Positive normal points toward increasing track
# radius in map_to_curve(), not necessarily away from the plot origin.
ggchord_empty_obstacles <- function() data.frame(accver = character(), start = numeric(),
  end = numeric(), side = character(), normal_min = numeric(), normal_max = numeric(),
  priority = numeric(), source_layer = character())

ggchord_outline_half_mm <- function(layer, mapping = NULL) {
  colour <- layer$aes_params$colour
  mapped <- layer$ggchord_input_mapping$colour %||% mapping$colour
  if (is.null(colour) && is.null(mapped)) colour <- layer$geom$default_aes$colour
  if (length(colour) && all(is.na(colour) | colour == "transparent")) return(0)
  width <- layer$aes_params$linewidth %||% layer$geom$default_aes$linewidth
  if (is.numeric(width)) max(0, width[is.finite(width)]) / 2 else 0
}

ggchord_ribbon_gap_mm <- function(params, ids, layer = NULL, mapping = NULL) {
  gap <- process_sequence_param(params$front_gap_mm, ids, "gap", .15)
  if (!is.null(layer) && !is.null(params$front_gap_mm))
    gap <- gap + ggchord_outline_half_mm(layer, mapping)
  gap
}

ggchord_backbone_clearance <- function(plot, ids) {
  positive <- negative <- setNames(numeric(length(ids)), ids)
  instances <- plot$ggchord$display_instances
  for (layer in plot$layers) {
    if (!identical(layer$ggchord_type, "seq")) next
    params <- layer$ggchord_params
    thickness <- params$seq_thickness_mm
    offset <- params$feature_position$offset %||% 0
    if (!is.null(instances)) {
      thickness <- ggchord_instance_param(thickness, instances)
      offset <- ggchord_instance_param(offset, instances)
    }
    half <- if (identical(params$seq_style, "single")) 0 else
      process_sequence_param(thickness, ids, "thickness", 1) / 2
    shift <- ggchord_resolve_position_offset(offset, ids, "offset", "geom_seq()")
    shift <- vapply(shift, function(x) unname(x[1L]), numeric(1))
    mapping <- if (isTRUE(layer$inherit.aes)) plot$mapping else NULL
    edge <- half + ggchord_outline_half_mm(layer, mapping)
    positive[] <- pmax(positive, shift + edge)
    negative[] <- pmax(negative, edge - shift)
  }
  list(positive = positive, negative = negative)
}

# Keep only changes of slope in each outline. Duplicate along-axis positions
# represent the vertical edges of a glyph and retain both extremes.
ggchord_obstacle_profile <- function(pos, offset) {
  tolerance <- max(1, abs(pos)) * 1e-12
  at <- order(pos)
  groups <- cumsum(c(TRUE, diff(pos[at]) > tolerance))
  stations <- pos[at][!duplicated(groups)]
  pos[at] <- stations[groups]
  if (length(stations) < 2L) return(NULL)
  edges <- c(seq_along(pos)[-1L], 1L)
  envelope <- vapply(stations, function(x) {
    hits <- x >= pmin(pos, pos[edges]) & x <= pmax(pos, pos[edges])
    t <- (x - pos[hits]) / (pos[edges[hits]] - pos[hits])
    values <- offset[hits] + t * (offset[edges[hits]] - offset[hits])
    range(values[is.finite(values)], offset[pos == x])
  }, numeric(2))
  turns <- function(y) which(abs(diff(diff(y) / diff(stations))) > 1e-7) + 1L
  keep <- sort(unique(c(1L, turns(envelope[1L, ]), turns(envelope[2L, ]), length(stations))))
  cbind(pos = stations[keep], low = envelope[1L, keep], high = envelope[2L, keep])
}

# Recover a continuous station along the rendered backbone. Projecting only
# onto the sample vertices can merge distinct shoulders of a short arrow.
ggchord_obstacle_station <- function(xy, ref, radius, path) {
  if (isTRUE(ref$circle)) {
    angles <- atan2(xy[, 2L], xy[, 1L])
    middle <- mean(range(ref$angles))
    angles <- angles + 2 * pi * round((middle - angles) / (2 * pi))
    return(cbind(angle = angles,
      offset = 2 * ref$r0 - radius - sqrt(rowSums(xy^2))))
  }
  base <- path$base
  nearby <- base[path$coarse, , drop = FALSE]
  t(vapply(seq_len(nrow(xy)), function(i) {
    point <- xy[i, ]
    distance <- (nearby[, 1L] - point[1L])^2 + (nearby[, 2L] - point[2L])^2
    reach <- (sqrt(min(distance)) + 10 * path$step)^2
    local <- path$coarse[distance <= reach]
    at <- seq.int(max(1L, min(local) - 10L), min(nrow(base) - 1L, max(local) + 10L))
    delta <- path$delta[at, , drop = FALSE]
    relative <- cbind(point[1L] - base[at, 1L], point[2L] - base[at, 2L])
    fraction <- pmin(1, pmax(0, rowSums(relative * delta) / path$length2[at]))
    rest <- relative - fraction * delta
    best <- which.min(rowSums(rest^2))
    j <- at[best]
    f <- fraction[best]
    normal <- (1 - f) * path$normal[j, ] + f * path$normal[j + 1L, ]
    c(angle = (1 - f) * ref$angles[j] + f * ref$angles[j + 1L],
      offset = sum(rest[best, ] * normal))
  }, numeric(2)))
}

ggchord_entity_obstacles <- function(layer, layout, input) {
  if (isTRUE(layout$collinear)) {
    return(ggchord_collinear_entity_obstacles(layer, layout, input))
  }
  polys <- layout$gene_polys
  if (is.null(input) || !nrow(input) || !nrow(polys)) return(ggchord_empty_obstacles())
  refs <- layout$sequence_reference
  paths <- lapply(names(refs$refs), function(id) {
    ref <- refs$refs[[id]]
    if (isTRUE(ref$circle)) return(NULL)
    base <- map_to_curve_many(ref$angles, refs$radius[id], ref)
    outward <- map_to_curve_many(ref$angles, refs$radius[id] + 1, ref) - base
    coarse <- unique(c(seq.int(1L, nrow(base), by = 20L), nrow(base)))
    delta <- diff(base)
    length2 <- rowSums(delta^2)
    list(base = base, normal = outward, coarse = coarse, delta = delta,
      length2 = pmax(length2, .Machine$double.eps), step = sqrt(max(length2)))
  })
  names(paths) <- names(refs$refs)
  stroke_mm <- ggchord_outline_half_mm(layer)
  parts <- lapply(split(polys, polys$group), function(d) {
    row <- input[d$source_row[1], , drop = FALSE]
    id <- as.character(row$accver)
    ref <- refs$refs[[id]]
    # Project the actual rendered polygon, including lollipop heads and stems.
    xy <- ggchord_rotate_points(as.matrix(d[c("x", "y")]), -layout$rotation)
    station <- ggchord_obstacle_station(xy, ref, refs$radius[id], paths[[id]])
    offset <- station[, "offset"]
    f <- (station[, "angle"] - refs$starts[id]) / (refs$ends[id] - refs$starts[id])
    if (refs$orientation[id] != 1) f <- 1 - f
    size <- unname(refs$lens[id])
    pos <- unname(f * size)
    if (isTRUE(ref$periodic)) {
      pos <- pos - size * c(0, cumsum(round(diff(pos) / size)))
      pos <- pos - floor(min(pos) / size) * size
    }
    lo <- min(offset); hi <- max(offset)
    profile <- ggchord_obstacle_profile(pos, offset)
    extent <- range(pos)
    shifts <- 0
    # Keep a seam-crossing contour local to the two ends of a full circle.
    if (isTRUE(ref$periodic)) {
      if (extent[1L] < 0) shifts <- c(shifts, size)
      if (extent[2L] > size) shifts <- c(shifts, -size)
    }
    ggchord_rbind_fill(lapply(shifts, function(shift) {
      local <- profile
      if (!is.null(local)) local[, "pos"] <- local[, "pos"] + shift
      ends <- pmin(size, pmax(0, extent + shift))
      data.frame(accver = id, start = ends[1L], end = ends[2L],
        side = if (lo >= 0) "positive" else if (hi <= 0) "negative" else "both",
        normal_min = lo, normal_max = hi, priority = 1,
        source_layer = layer$ggchord_layer_id, stroke_mm = stroke_mm,
        profile = I(list(local)), stringsAsFactors = FALSE)
    }))
  })
  ggchord_rbind_fill(parts)
}

ggchord_obstacle_index <- function(obstacles, avoid) {
  if (identical(avoid, "none") || !NROW(obstacles)) return(list())
  split(obstacles, obstacles$accver)
}

ggchord_ribbon_base_gap <- function(clearance, id, direction, normal_mm, gap_mm = .15) {
  side <- if (direction >= 0) "positive" else "negative"
  ((clearance[[side]][[id]] %||% .45) + gap_mm) * normal_mm
}

# The outside of a disk swept along each outline segment. In the local
# (distance, normal) frame, a segment of slope m reaches its highest offset
# at the source station shifted by padding * m / sqrt(1 + m^2).
ggchord_front_outline <- function(profile, direction, along, padding) {
  stations <- profile[, "pos"]
  x <- (stations - stations[1L]) * along
  y <- if (direction > 0) profile[, "high"] else -profile[, "low"]
  slope <- diff(y) / diff(x)
  turn <- c(pi / 2, atan(slope), -pi / 2)
  fractions <- seq(0, 1, length.out = 7L)
  angles <- outer(turn[-length(turn)], 1 - fractions) + outer(turn[-1L], fractions)
  knots <- (x - padding * sin(angles)) / along + stations[1L]
  list(stations = stations, x = x, y = y, slope = slope,
    shift = padding * sin(turn[-c(1L, length(turn))]), knots = as.vector(knots))
}

ggchord_front_outline_at <- function(outline, positions, along, padding) {
  x <- outline$x
  y <- outline$y
  u <- (positions - outline$stations[1L]) * along
  peak <- rep(-Inf, length(u))
  for (i in seq_along(outline$slope)) {
    at <- which(u >= x[i] - padding & u <= x[i + 1L] + padding)
    if (!length(at)) next
    source <- pmax(x[i], pmin(u[at] + outline$shift[i], x[i + 1L]))
    height <- y[i] + (source - x[i]) * outline$slope[i] +
      sqrt(pmax(0, padding^2 - (source - u[at])^2))
    peak[at] <- pmax(peak[at], height)
  }
  peak
}

# The front is sampled at the backbone's usual cadence plus the knots of
# nearby obstacles. Ramps are defined in physical length, independently
# of the alignment span, and only touch the side facing the destination.
ggchord_connection_front <- function(bounds, obstacles, direction, avoid,
    base, normal_mm, along_per_position, n = 2L, gap_mm = .30, strength = .5) {
  positions <- seq(bounds[1L], bounds[2L], length.out = n)
  gap <- rep(base, n)
  if (!NROW(obstacles) || identical(avoid, "none") || direction == 0)
    return(list(position = positions, gap = gap))
  clearance <- if (direction > 0) obstacles$normal_max else -obstacles$normal_min
  stroke <- obstacles$stroke_mm %||% rep(0, NROW(obstacles))
  stroke[!is.finite(stroke)] <- 0
  # Small extra margin bounds the error from connecting samples of a corner arc.
  padding <- (gap_mm + stroke + .01) * normal_mm
  height <- pmax(0, clearance + padding - base)
  along <- max(along_per_position, .Machine$double.eps)
  transition <- (padding + pmin(2.4, pmax(.25, height / normal_mm)) *
    (.45 + 1.1 * strength) * normal_mm) / along
  ends <- range(bounds)
  keep <- which(clearance > 0 & height > 0 &
    obstacles$start - transition < ends[2L] & obstacles$end + transition > ends[1L])
  if (!length(keep)) return(list(position = positions, gap = gap))
  if (identical(avoid, "uniform")) {
    active <- keep[obstacles$start[keep] <= ends[2L] & obstacles$end[keep] >= ends[1L]]
    return(list(position = positions,
      gap = rep(max(base, clearance[active] + padding[active]), n)))
  }
  ramps <- c(.25, .50, .75, 1)
  points <- c(obstacles$start[keep], obstacles$end[keep],
    as.vector(outer(obstacles$start[keep], ramps, function(x, y) x - transition[keep] * y)),
    as.vector(outer(obstacles$end[keep], ramps, function(x, y) x + transition[keep] * y)))
  outlines <- lapply(keep, function(i) {
    profile <- obstacles$profile[[i]]
    if (is.matrix(profile)) ggchord_front_outline(profile, direction, along, padding[i])
  })
  points <- c(points, unlist(lapply(outlines, `[[`, "knots"), use.names = FALSE))
  if (n > 1L) positions <- sort(unique(c(positions,
    points[points > ends[1L] & points < ends[2L]])), decreasing = diff(bounds) < 0)
  bump <- numeric(length(positions))
  for (j in seq_along(keep)) {
    i <- keep[j]
    distance <- pmax(obstacles$start[i] - positions, positions - obstacles$end[i], 0)
    at <- which(distance < transition[i])
    if (!length(at)) next
    drop <- distance[at] / transition[i]
    weight <- 1 - drop^3 * (10 - 15 * drop + 6 * drop^2)
    outline <- outlines[[j]]
    if (is.null(outline)) {
      rise <- height[i] * weight
    } else {
      at_position <- positions[at]
      peak <- ggchord_front_outline_at(outline, at_position, along, padding[i])
      extension <- ifelse(at_position <= outline$stations[1L], outline$y[1L],
        tail(outline$y, 1L)) + padding[i] - base
      outside <- distance[at] > 0
      rise <- pmax(0, peak - base, ifelse(outside, extension * weight, 0))
    }
    # A compact smooth maximum rounds neighbouring turns while staying
    # outside both envelopes. It contributes nothing off the obstacles.
    previous <- bump[at]
    width <- pmin((.08 + .24 * strength) * normal_mm, previous, rise)
    blend <- pmax(0, width - abs(previous - rise))^2 / (4 * pmax(width, .Machine$double.eps))
    bump[at] <- pmax(previous, rise) + blend
  }
  list(position = positions, gap = base + bump)
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
    input <- ggchord_resolve_layer_input(layer, ggchord_annotation_fallback(plot$ggchord$data))
    out <- layer$ggchord_obstacle_provider(layer, layout, input)
    ggchord_require_columns(out, names(ggchord_empty_obstacles()), "obstacle provider")
    numeric <- c("start", "end", "normal_min", "normal_max", "priority")
    if (any(!vapply(out[numeric], is.numeric, logical(1))) || any(!is.finite(as.matrix(out[numeric]))) ||
        any(out$start > out$end | out$normal_min > out$normal_max) ||
        any(!out$side %in% c("positive", "negative", "both"))) ggchord_stop("Invalid obstacle provider ranges")
    out
  })
  ggchord_rbind_fill(parts)
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
  start <- xy[-n, , drop = FALSE]
  tip <- xy[-1L, , drop = FALSE]
  low <- pmin(start, tip)
  high <- pmax(start, tip)
  for (i in seq_len(n - 3L)) {
    j <- which(seq_len(n - 1L) > i + 1L &
      high[, 1L] >= low[i, 1L] & low[, 1L] <= high[i, 1L] &
      high[, 2L] >= low[i, 2L] & low[, 2L] <= high[i, 2L])
    if (!length(j)) next
    u <- a[i, ]; v <- a[j, , drop = FALSE]
    den <- u[1L] * v[, 2L] - u[2L] * v[, 1L]
    keep <- abs(den) >= 1e-12
    if (!any(keep)) next
    w <- sweep(start[j[keep], , drop = FALSE], 2L, start[i, ], "-")
    v <- v[keep, , drop = FALSE]; den <- den[keep]
    t <- (w[, 1L] * v[, 2L] - w[, 2L] * v[, 1L]) / den
    s <- (w[, 1L] * u[2L] - w[, 2L] * u[1L]) / den
    if (any(t > 1e-8 & t < 1 - 1e-8 & s > 1e-8 & s < 1 - 1e-8)) return(TRUE)
  }
  FALSE
}
