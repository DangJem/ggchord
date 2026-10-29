# Small explicit band allocator. It shares source frames, not pseudo radii.
ggchord_track_dimension <- function(value, seqs, name) {
  if (grid::is.unit(value)) {
    multiplier <- switch(grid::unitType(value), mm = 1, cm = 10,
      inches = 25.4, points = 25.4 / 72.27, bigpts = 25.4 / 72)
    return(list(native = stats::setNames(rep(0, length(seqs)), seqs),
      mm = stats::setNames(rep(as.numeric(value) * multiplier, length(seqs)), seqs)))
  }
  native <- process_sequence_param(value, seqs, name, NA_real_)
  if (any(!is.finite(native))) ggchord_stop("Track ", name, " must be specified for each displayed sequence")
  list(native = native, mm = stats::setNames(rep(0, length(seqs)), seqs))
}

ggchord_track_frame <- function(pos, id, layout) {
  ref <- layout$sequence_reference
  if (isTRUE(layout$collinear)) return(list(base = cbind(
    ggchord_collinear_x(pos, id, ref), rep(unname(ref$y[id]), length(pos))),
    normal = cbind(rep(0, length(pos)), rep(1, length(pos))), mapper_sign = 1))
  fraction <- pos / unname(ref$lens[id])
  if (ref$orientation[id] != 1) fraction <- 1 - fraction
  angles <- unname(ref$starts[id]) + fraction * unname(ref$ends[id] - ref$starts[id])
  reference <- ref$refs[[id]]
  # Continuous interpolation avoids quantizing short quantitative intervals
  # to the nearest precomputed backbone sample.
  path <- reference$path
  before <- pmax(1L, seq_len(nrow(path)) - 1L)
  after <- pmin(nrow(path), seq_len(nrow(path)) + 1L)
  dx <- path$x[after] - path$x[before]; dy <- path$y[after] - path$y[before]
  norm <- sqrt(dx^2 + dy^2)
  if (any(norm <= 0)) ggchord_stop("Track backbone has an undefined local normal")
  nx <- -dy / norm; ny <- dx / norm
  interpolate <- function(value) stats::approx(reference$angles, value, xout = angles, rule = 2)$y
  normal <- cbind(interpolate(nx), interpolate(ny))
  normal <- normal / sqrt(rowSums(normal^2))
  base <- cbind(interpolate(path$x), interpolate(path$y)) + normal *
    (unname(ref$radius[id]) - reference$r0)
  sample <- map_to_curve_many(reference$angles, unname(ref$radius[id]), reference)
  outward <- map_to_curve_many(reference$angles, unname(ref$radius[id]) + 1, reference) - sample
  sign <- if (stats::median(rowSums(sample * outward)) < 0) -1 else 1
  list(base = ggchord_rotate_points(base, layout$rotation %||% 0),
    normal = ggchord_rotate_points(normal * sign, layout$rotation %||% 0), mapper_sign = sign)
}

ggchord_track_units <- function(normal, layout) {
  context <- ggchord_render_context()
  units <- context$units_per_mm
  if (length(units) != 2L || any(!is.finite(units)) || any(units <= 0))
    units <- rep((layout$text_units_per_inch %||% .35) / 25.4, 2L)
  1 / sqrt((normal[, 1L] / units[1L])^2 + (normal[, 2L] / units[2L])^2)
}

ggchord_track_point <- function(pos, id, offset_native, offset_mm, layout) {
  frame <- ggchord_track_frame(pos, id, layout)
  offset <- offset_native + offset_mm * ggchord_track_units(frame$normal, layout)
  frame$base + frame$normal * offset
}

ggchord_empty_track_geometry <- function() data.frame(x = numeric(), y = numeric(),
  group = integer(), source_row = integer(), accver = character(), track_id = character(),
  label = character(), .component = character(), stringsAsFactors = FALSE)

ggchord_track_annotation_range <- function(registry, seq, layout) {
  points <- list()
  registry[[".automatic-axis"]] <- list(lines = layout$axis_lines, ticks = layout$axis_ticks)
  for (owner in registry) for (data in owner) {
    if (!is.data.frame(data) || !nrow(data) || !"accver" %in% names(data)) next
    data <- data[!is.na(data$accver) & as.character(data$accver) == seq, , drop = FALSE]
    if (!nrow(data)) next
    for (pair in list(c("x", "y"), c("x0", "y0"), c("x1", "y1")))
      if (all(pair %in% names(data))) points[[length(points) + 1L]] <- as.matrix(data[pair])
    if ("label" %in% names(data)) {
      text <- data[!is.na(data$label) & nzchar(as.character(data$label)), , drop = FALSE]
      if (nrow(text)) {
        text$text <- as.character(text$label)
        text$text_x <- text$text_x %||% text$x %||% text$label_x
        text$text_y <- text$text_y %||% text$y %||% text$label_y
        text$text_angle <- text$text_angle %||% text$label_angle %||% 0
        if (length(text$text_x) && length(text$text_y)) {
          text$size <- text$size %||% 2.5
          boxes <- ggchord_text_boxes(text, units_per_inch = layout$text_units_per_inch)
          points[[length(points) + 1L]] <- rbind(cbind(boxes$xmin, boxes$ymin),
            cbind(boxes$xmax, boxes$ymin), cbind(boxes$xmax, boxes$ymax), cbind(boxes$xmin, boxes$ymax))
        }
      }
    }
  }
  if (!length(points)) return(c(0, 0))
  xy <- do.call(rbind, points); xy <- xy[is.finite(xy[, 1]) & is.finite(xy[, 2]), , drop = FALSE]
  frame <- ggchord_track_frame(seq(0, layout$sequence_reference$lens[[seq]], length.out = 129L), seq, layout)
  if (!nrow(xy)) return(c(0, 0))
  nearest <- vapply(seq_len(nrow(xy)), function(i) which.min(rowSums((frame$base -
    matrix(xy[i, ], nrow(frame$base), 2L, byrow = TRUE))^2)), integer(1))
  offsets <- rowSums((xy - frame$base[nearest, , drop = FALSE]) * frame$normal[nearest, , drop = FALSE])
  range(c(0, offsets))
}

ggchord_track_polygons_overlap <- function(a, b) {
  if (max(a$x) <= min(b$x) || max(b$x) <= min(a$x) ||
      max(a$y) <= min(b$y) || max(b$y) <= min(a$y)) return(FALSE)
  if (ggchord_point_in_polygon(a$x[1L], a$y[1L], b) ||
      ggchord_point_in_polygon(b$x[1L], b$y[1L], a)) return(TRUE)
  aa <- rbind(as.matrix(a), as.matrix(a[1L, ])); bb <- rbind(as.matrix(b), as.matrix(b[1L, ]))
  for (i in seq_len(nrow(a))) {
    candidates <- which(pmax(bb[-1L, 1], bb[-nrow(bb), 1]) >= min(aa[i:(i + 1L), 1]) &
      pmin(bb[-1L, 1], bb[-nrow(bb), 1]) <= max(aa[i:(i + 1L), 1]) &
      pmax(bb[-1L, 2], bb[-nrow(bb), 2]) >= min(aa[i:(i + 1L), 2]) &
      pmin(bb[-1L, 2], bb[-nrow(bb), 2]) <= max(aa[i:(i + 1L), 2]))
    for (j in candidates) if (ggchord_segment_intersects(aa[i, ], aa[i + 1L, ], bb[j, ], bb[j + 1L, ])) return(TRUE)
  }
  FALSE
}

ggchord_layout_tracks <- function(plot, layout, registry, inputs) {
  layers <- Filter(function(l) identical(l$ggchord_type, "track"), plot$layers)
  seqs <- layout$seqs; ref <- layout$sequence_reference
  data <- lapply(layers, function(l) ggchord_resolve_layer_input(l, plot$ggchord$data$seq_data))
  for (i in seq_along(data)) {
    d <- data[[i]]; p <- layers[[i]]$ggchord_params
    if (is.null(d)) ggchord_stop("Track layers require explicit data")
    cols <- if (p$track_kind == "line") "pos" else c("start", "end")
    ggchord_require_columns(d, c("accver", cols), "track layer")
    d$accver <- ggchord_conversion_ids(d$accver, "track layer")
    if (any(!d$accver %in% seqs)) ggchord_stop("Track input contains an unknown displayed accession; focus to window IDs explicitly")
    for (name in cols) if (!is.numeric(d[[name]]) || any(!is.finite(d[[name]])) || any(d[[name]] != floor(d[[name]])))
      ggchord_stop("Track genomic coordinates must be finite integers")
    if (p$track_kind != "tile") {
      ggchord_require_columns(d, "track_value", "quantitative track")
      if (!is.numeric(d$track_value) || any(is.infinite(d$track_value))) ggchord_stop("Track values must be numeric; NA is allowed but Inf is not")
    } else if ("track_value" %in% names(d)) ggchord_stop("Tile tracks use track_fill, not numerical height track_value")
    if (!"group" %in% names(d)) d$group <- rep(1L, nrow(d))
    if (anyNA(d$group)) ggchord_stop("Track series group IDs must be non-missing")
    if (isTRUE(p$closed) && (!isTRUE(layout$circular) || (plot$coordinates$circular_gap %||% 0) != 0))
      ggchord_stop("closed=TRUE requires coord_circular(gap=0)")
    if (p$track_kind == "line") {
      groups <- split(seq_len(nrow(d)), interaction(match(d$accver, unique(d$accver)),
        match(d$group, unique(d$group)), drop = TRUE))
      for (g in groups) if (length(g) > 1L && any(diff(d$pos[g]) <= 0))
        ggchord_stop("Line positions must increase strictly within each accession/series; sort or resolve duplicates explicitly")
    }
    data[[i]] <- d
  }
  ids <- vapply(layers, function(l) l$ggchord_params$track_position$track, character(1))
  unique_ids <- unique(ids)
  declarations <- list(); scales <- list(); trained <- list()
  for (id in unique_ids) {
    members <- which(ids == id)
    pos <- layers[[members[1L]]]$ggchord_params$track_position
    width <- ggchord_track_dimension(pos$width, seqs, "width")
    gap <- ggchord_track_dimension(pos$gap, seqs, "gap")
    for (i in members) {
      other <- layers[[i]]$ggchord_params$track_position
      if (other$side != pos$side || other$track_order != pos$track_order ||
          !isTRUE(all.equal(width, ggchord_track_dimension(other$width, seqs, "width"))) ||
          !isTRUE(all.equal(gap, ggchord_track_dimension(other$gap, seqs, "gap"))))
        ggchord_stop("Conflicting placement declarations for track '", id, "'")
    }
    declarations[[id]] <- list(position = pos, width = width, gap = gap,
      layers = vapply(layers[members], function(l) l$ggchord_layer_id, character(1)))
    numeric <- members[vapply(members, function(i) layers[[i]]$ggchord_params$track_kind != "tile", logical(1))]
    if (!length(numeric)) next
    key <- ggchord_track_aes("value", id)
    sc <- plot$scales$get_scales(key)
    if (is.null(sc)) sc <- scale_track_value_continuous(id)
    sc <- sc$clone()
    values <- unlist(lapply(numeric, function(i) c(data[[i]]$track_value,
      if (layers[[i]]$ggchord_params$track_kind == "bar") layers[[i]]$ggchord_params$baseline)), use.names = FALSE)
    finite <- values[!is.na(values)]
    transformed <- suppressWarnings(sc$transform(finite))
    if (any(!is.finite(transformed))) ggchord_stop("Track '", id, "' has values/baseline incompatible with its scale transform")
    sc$train(transformed)
    has_domain <- length(finite) > 0L || !is.null(sc$limits)
    if (has_domain && (length(sc$get_limits()) != 2L || any(!is.finite(sc$get_limits()))))
      ggchord_stop("Track numerical scale limits must define a finite domain")
    for (i in numeric) if (layers[[i]]$ggchord_params$track_kind == "bar" &&
        (!is.finite(sc$map(sc$transform(layers[[i]]$ggchord_params$baseline))) ||
         sc$map(sc$transform(layers[[i]]$ggchord_params$baseline)) < 0 ||
         sc$map(sc$transform(layers[[i]]$ggchord_params$baseline)) > 1))
      ggchord_stop("Track '", id, "' scale excludes its baseline; change limits or explicitly choose an oob policy")
    scales[[id]] <- sc; trained[[id]] <- has_domain
  }
  # Existing entities keep their fixed placement. New bands start beyond
  # their occupied normal range, and further track IDs accumulate outward.
  bands <- list(); band_polygons <- list(); counters <- list()
  annotation_ranges <- stats::setNames(lapply(seqs, function(seq)
    ggchord_track_annotation_range(registry, seq, layout)), seqs)
  ordered <- order(vapply(declarations, function(x) x$position$track_order, numeric(1)), names(declarations))
  for (id in names(declarations)[ordered]) {
    declaration <- declarations[[id]]; p <- declaration$position
    direction <- if (p$side == "outside") 1 else -1
    for (seq in seqs) {
      counter <- paste(seq, p$side, sep = "\r")
      start <- counters[[counter]]
      if (is.null(start)) {
        occupied <- layout$obstacles
        occupied <- if (is.null(occupied)) data.frame() else occupied[occupied$accver == seq, , drop = FALSE]
        edge <- if (nrow(occupied)) max(if (direction > 0) occupied$normal_max else -occupied$normal_min, 0) else 0
        annotation <- annotation_ranges[[seq]]
        edge <- max(edge, if (direction > 0) annotation[2L] else -annotation[1L])
        start <- list(native = max(.03, edge), mm = 0)
      }
      near_native <- start$native + unname(declaration$gap$native[seq])
      near_mm <- start$mm + unname(declaration$gap$mm[seq])
      width_native <- unname(declaration$width$native[seq]); width_mm <- unname(declaration$width$mm[seq])
      sample <- seq(0, ref$lens[[seq]], length.out = 65L)
      frame <- ggchord_track_frame(sample, seq, layout)
      lower <- ggchord_track_point(sample, seq, direction * near_native, direction * near_mm, layout)
      upper <- ggchord_track_point(sample, seq, direction * (near_native + width_native), direction * (near_mm + width_mm), layout)
      if (!isTRUE(layout$collinear) &&
          (ggchord_front_invalid(lower, frame$base) || ggchord_front_invalid(upper, frame$base)))
        ggchord_stop("Track '", id, "' band folds or self-intersects; reduce its width/gap or change sequence curvature/radius")
      offsets <- c(rowSums((lower - frame$base) * frame$normal), rowSums((upper - frame$base) * frame$normal))
      if (isTRUE(layout$circular) && any(rowSums(upper * frame$base) <= 0))
        ggchord_stop("Track inner band reaches the circular centre; reduce width or move tracks outside")
      row <- data.frame(track_id = id, accver = seq, side = p$side,
        order = p$track_order, direction = direction,
        offset_native = direction * near_native, offset_mm = direction * near_mm,
        width_native = width_native, width_mm = width_mm,
        normal_min = min(offsets), normal_max = max(offsets), mapper_sign = frame$mapper_sign,
        measurement_status = (ggchord_render_context()$status %||% "estimated"),
        domain_status = if (is.null(scales[[id]])) "categorical_or_fixed" else if (isTRUE(trained[[id]])) "trained" else "untrained",
        unit = scales[[id]]$ggchord_track_unit %||% NA_character_, stringsAsFactors = FALSE)
      row$source_layers <- I(list(declaration$layers))
      row$limits <- I(list(if (isTRUE(trained[[id]])) scales[[id]]$get_transformation()$inverse(scales[[id]]$get_limits()) else numeric()))
      bands[[length(bands) + 1L]] <- row
      polygon <- rbind(lower, upper[nrow(upper):1L, , drop = FALSE])
      band_polygons[[length(band_polygons) + 1L]] <- data.frame(x = polygon[, 1], y = polygon[, 2])
      counters[[counter]] <- list(native = near_native + width_native, mm = near_mm + width_mm)
    }
  }
  band_table <- ggchord_rbind_fill(bands)
  if (!isTRUE(layout$collinear) && length(seqs) > 1L) {
    for (i in seq_len(nrow(band_table))) for (other in setdiff(seqs, band_table$accver[i])) {
      path <- ggchord_track_frame(seq(0, ref$lens[[other]], length.out = 65L), other, layout)$base
      polygon <- band_polygons[[i]]
      if (any(vapply(seq_len(nrow(path)), function(j)
          ggchord_point_in_polygon(path[j, 1], path[j, 2], polygon), logical(1))))
        ggchord_stop("Track band reaches another sequence backbone; change sequence layout or reduce track width/gap")
    }
    if (nrow(band_table) > 1L) for (i in seq_len(nrow(band_table) - 1L)) for (j in seq.int(i + 1L, nrow(band_table))) {
      if (band_table$accver[i] == band_table$accver[j]) next
      if (ggchord_track_polygons_overlap(band_polygons[[i]], band_polygons[[j]]))
        ggchord_stop("Tracks from different sequences intersect; change curvature/layout or reduce track widths")
    }
  }
  if (isTRUE(layout$collinear) && nrow(band_table)) {
    # A band may not spill over another backbone or another row's full band.
    bounds <- layout$sequence_layout
    for (i in seq_len(nrow(band_table))) {
      b <- band_table[i, ]; at <- match(b$accver, bounds$accver)
      lo <- bounds$y[at] + b$normal_min; hi <- bounds$y[at] + b$normal_max
      others <- which(bounds$y != bounds$y[at] & bounds$xmax >= bounds$xmin[at] & bounds$xmin <= bounds$xmax[at])
      if (any(bounds$y[others] >= lo & bounds$y[others] <= hi))
        ggchord_stop("Track band reaches another backbone; increase seq_gap or reduce track widths")
    }
    if (nrow(band_table) > 1L) for (i in seq_len(nrow(band_table) - 1L)) for (j in seq.int(i + 1L, nrow(band_table))) {
      a <- band_table[i, ]; b <- band_table[j, ]
      ai <- match(a$accver, bounds$accver); bi <- match(b$accver, bounds$accver)
      if (bounds$y[ai] == bounds$y[bi] || bounds$xmax[ai] < bounds$xmin[bi] || bounds$xmax[bi] < bounds$xmin[ai]) next
      if (min(bounds$y[ai] + a$normal_max, bounds$y[bi] + b$normal_max) >
          max(bounds$y[ai] + a$normal_min, bounds$y[bi] + b$normal_min))
        ggchord_stop("Tracks on different display rows overlap; increase seq_gap or reduce track widths")
    }
  }
  diagnostics <- list()
  for (i in seq_along(layers)) {
    layer <- layers[[i]]; p <- layer$ggchord_params; d <- data[[i]]
    id <- unname(ids[i]); owner <- layer$ggchord_layer_id; geometry <- list(); group <- 0L
    sc <- scales[[id]]
    mapped <- if (p$track_kind == "tile") rep(1, nrow(d)) else if (isTRUE(trained[[id]]))
      sc$map(suppressWarnings(sc$transform(d$track_value))) else rep(NA_real_, nrow(d))
    if (any(!is.na(mapped) & (!is.finite(mapped) | mapped < 0 | mapped > 1)))
      ggchord_stop("Track scale oob policy allows values outside its physical band; use limits with censor or squish")
    reason <- ifelse(is.na(mapped), ifelse(is.na(d$track_value %||% rep(NA_real_, nrow(d))), "missing_value", "outside_scale_limits"), "retained")
    clipped <- rep(FALSE, nrow(d))
    band_for <- function(seq) band_table[band_table$track_id == id & band_table$accver == seq, , drop = FALSE]
    at_value <- function(pos, seq, value) {
      b <- band_for(seq)
      ggchord_track_point(pos, seq,
        b$offset_native + b$direction * value * b$width_native,
        b$offset_mm + b$direction * value * b$width_mm, layout)
    }
    add <- function(xy, source, seq, component) {
      group <<- group + 1L
      geometry[[length(geometry) + 1L]] <<- data.frame(x = xy[, 1L], y = xy[, 2L],
        group = group, source_row = rep_len(source, nrow(xy)), accver = unname(seq),
        track_id = id, label = NA_character_, .component = component, stringsAsFactors = FALSE)
    }
    if (p$track_kind == "line") {
      usable <- !is.na(mapped)
      outside <- d$pos < 1 | d$pos > unname(ref$lens[d$accver])
      if (any(outside) && p$bounds == "error") ggchord_stop("Track sample coordinate outside sequence length")
      usable[outside] <- FALSE; reason[outside] <- "outside_window"; clipped[outside] <- TRUE
      series <- split(seq_len(nrow(d)), interaction(match(d$accver, unique(d$accver)),
        match(d$group, unique(d$group)), drop = TRUE))
      for (rows in series) {
        runs <- split(rows[usable[rows]], cumsum(!usable[rows])[usable[rows]])
        for (run in runs) {
          if (length(run) < 2L) { reason[run] <- "isolated_sample"; next }
          seq <- d$accver[run[1L]]; vertices <- list(); sources <- integer()
          for (j in seq_len(length(run) - 1L)) {
            a <- run[j]; b <- run[j + 1L]
            count <- if (isTRUE(layout$collinear)) 2L else 16L
            pos <- seq(d$pos[a], d$pos[b], length.out = count)
            value <- seq(mapped[a], mapped[b], length.out = count)
            vertices[[length(vertices) + 1L]] <- at_value(pos, seq, value)
            sources <- c(sources, rep(a, count - 1L), b)
          }
          if (isTRUE(p$closed) && length(run) == length(rows) && all(usable[rows])) {
            seam <- seq(d$pos[run[length(run)]], d$pos[run[1L]] + ref$lens[[seq]], length.out = 24L)
            value <- seq(mapped[run[length(run)]], mapped[run[1L]], length.out = 24L)
            vertices[[length(vertices) + 1L]] <- at_value(((seam - 1) %% ref$lens[[seq]]) + 1, seq, value)
            sources <- c(sources, rep(run[length(run)], length(seam)))
          }
          add(do.call(rbind, vertices), sources, seq, "track_line")
        }
      }
    } else for (j in seq_len(nrow(d))) {
      if (is.na(mapped[j])) next
      seq <- d$accver[j]; a <- d$start[j]; z <- d$end[j]; len <- ref$lens[[seq]]
      if ((min(a, z) < 1 || max(a, z) > len) && p$bounds == "error") ggchord_stop("Track interval outside sequence length")
      if (max(a, z) < 1 || min(a, z) > len) { reason[j] <- "outside_window"; next }
      if (min(a, z) < 1 || max(a, z) > len) { a <- max(1, min(a, len)); z <- max(1, min(z, len)); clipped[j] <- TRUE }
      if (a > z && !isTRUE(layout$circular)) ggchord_stop("Origin-crossing track intervals require circular coordinates or explicit window splitting")
      pieces <- if (a > z) list(c(a, len), c(1, z)) else list(c(a, z))
      baseline <- if (p$track_kind == "bar") sc$map(sc$transform(p$baseline)) else 0
      for (part in pieces) {
        edges <- c(part[1L] - 1, part[2L])
        pos <- seq(edges[1L], edges[2L], length.out = if (isTRUE(layout$collinear)) 2L else 32L)
        xy <- rbind(at_value(pos, seq, baseline), at_value(rev(pos), seq, mapped[j]))
        add(xy, j, seq, "track_poly")
      }
    }
    out <- if (length(geometry)) ggchord_rbind_fill(geometry) else ggchord_empty_track_geometry()
    registry[[owner]] <- list(track = ggchord_attach_input_columns(out, d))
    inputs[[owner]] <- list(track = d)
    diagnostics[[owner]] <- data.frame(layer_id = rep(owner, nrow(d)), track_id = rep(id, nrow(d)),
      source_row = seq_len(nrow(d)), accver = d$accver,
      retained = seq_len(nrow(d)) %in% out$source_row,
      reason = reason, coordinate_clipped = clipped,
      baseline = rep(if (p$track_kind == "bar") p$baseline else NA_real_, nrow(d)), stringsAsFactors = FALSE)
  }
  obstacle <- band_table
  obstacle$start <- 1; obstacle$end <- unname(ref$lens[obstacle$accver])
  if (!isTRUE(layout$collinear)) {
    low <- pmin(obstacle$normal_min * obstacle$mapper_sign, obstacle$normal_max * obstacle$mapper_sign)
    high <- pmax(obstacle$normal_min * obstacle$mapper_sign, obstacle$normal_max * obstacle$mapper_sign)
    obstacle$normal_min <- low; obstacle$normal_max <- high
  }
  obstacle$side <- ifelse(obstacle$normal_min >= 0, "positive", ifelse(obstacle$normal_max <= 0, "negative", "both"))
  obstacle$priority <- 1; obstacle$source_layer <- paste0("track:", obstacle$track_id)
  layout$track_registry <- band_table; layout$track_scales <- scales
  layout$track_obstacles <- obstacle[names(ggchord_empty_obstacles())]
  layout$track_diagnostics <- ggchord_rbind_fill(diagnostics)
  list(layout = layout, registry = registry, inputs = inputs)
}

ggchord_track_axes <- function(layout, registry) {
  axes <- list()
  for (id in names(layout$track_scales)) {
    sc <- layout$track_scales[[id]]
    band <- layout$track_registry[layout$track_registry$track_id == id, , drop = FALSE][1L, ]
    if (!isTRUE(sc$ggchord_track_axis) || band$domain_status != "trained") next
    owner <- band$source_layers[[1L]][1L]
    # One explicit numerical reference per shared track domain. Automatic
    # breaks use its endpoints; explicit user breaks/labels remain untouched.
    breaks <- if (ggplot2::is_waiver(sc$breaks)) unique(sc$get_limits()) else sc$get_breaks()
    labels <- sc$get_labels(breaks)
    fractions <- sc$map(breaks)
    keep <- is.finite(fractions) & fractions >= 0 & fractions <= 1
    fractions <- fractions[keep]; labels <- labels[keep]
    seq <- band$accver
    anchor <- if (isTRUE(layout$collinear)) 1 else 1 +
      (match(id, names(layout$track_scales)) - 1) *
        (layout$sequence_reference$lens[[seq]] - 1) / max(4L, length(layout$track_scales))
    point <- function(f) ggchord_track_point(rep(anchor, length(f)), seq,
      band$offset_native + band$direction * f * band$width_native,
      band$offset_mm + band$direction * f * band$width_mm, layout)
    ends <- point(c(0, 1)); ticks <- point(fractions)
    normal <- ggchord_track_frame(anchor, seq, layout)$normal[1L, ]
    tangent <- matrix(c(-normal[2L], normal[1L]), nrow = 1L)
    spacing <- ggchord_track_units(tangent, layout)
    label_points <- ticks + matrix(rep(as.numeric(tangent) * spacing * 2.5,
      length(fractions)), ncol = 2L, byrow = TRUE)
    line <- data.frame(x = ends[, 1], y = ends[, 2], group = -1L,
      source_row = NA_integer_, accver = seq, track_id = id, label = NA_character_,
      .component = "track_axis_line")
    text <- ggchord_empty_track_geometry()
    if (length(labels)) text <- data.frame(x = label_points[, 1], y = label_points[, 2],
      group = -seq_along(fractions) - 1L, source_row = NA_integer_, accver = seq,
      track_id = id, label = as.character(labels), .component = "track_axis_text")
    title <- if (ggplot2::is_waiver(sc$name)) id else as.character(sc$name)
    if (length(title) && !is.na(title) && nzchar(title)) {
      value <- point(.5) + tangent * spacing * 5.5
      unit <- sc$ggchord_track_unit
      if (!is.null(unit) && nzchar(unit)) title <- paste0(title, " (", unit, ")")
      title_row <- data.frame(x = value[1, 1], y = value[1, 2], group = -1000L,
        source_row = NA_integer_, accver = seq, track_id = id, label = title,
        .component = "track_axis_text")
      text <- ggchord_rbind_fill(list(text, title_row))
    }
    out <- ggchord_rbind_fill(list(line, text))
    # Add blank input columns so visual mappings can evaluate on axis rows.
    base <- registry[[owner]]$track
    out <- ggchord_rbind_fill(list(base, out))
    registry[[owner]]$track <- out
    axes[[id]] <- ggchord_rbind_fill(list(line, text))
  }
  list(registry = registry, axes = ggchord_rbind_fill(axes) %||% data.frame())
}
