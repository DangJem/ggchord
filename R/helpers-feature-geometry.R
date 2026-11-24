# Shared interval and shape helpers for geom_feature(preset = 'gene') and geom_feature().

ggchord_feature_classes <- function(type) {
  key <- tolower(trimws(as.character(type)))
  key <- gsub("[^a-z0-9]+", "_", key)
  key <- gsub("^_+|_+$", "", key)
  aliases <- c(
    gene = "gene", cds = "gene", resistance_gene = "resistance_gene",
    rep_origin = "origin", replication_origin = "origin", ori = "origin",
    origin = "origin", promoter = "promoter",
    primer = "primer", primer_bind = "primer",
    protein_bind = "binding_site", binding_site = "binding_site",
    operator = "operator", regulatory = "regulatory",
    regulatory_region = "regulatory", mcs = "mcs",
    misc_feature = "misc_feature", repeat_region = "misc_feature",
    ltr = "misc_feature", misc_rna = "regulatory",
    terminator = "regulatory", rbs = "regulatory",
    polya_signal = "regulatory", poly_a_signal = "regulatory",
    enhancer = "regulatory"
  )
  out <- unname(aliases[key])
  out[is.na(out) | !nzchar(out)] <- "misc_feature"
  out
}

ggchord_feature_width_factor <- function(shape) {
  # Shape changes topology, not radial emphasis. Every feature type occupies
  # the same band width so promoters, primers and markers remain legible beside
  # long arrows and do not imply a weaker semantic hierarchy.
  rep(1, length(shape))
}

ggchord_block_arrow_shapes <- c("arrow", "arrow_flush")
ggchord_circular_arrow_head_ratio <- 1.42

ggchord_circular_arrow_band <- function(circular, feature_role) {
  isTRUE(circular) && isTRUE(feature_role %in% c("feature", "gene"))
}

# Minimum visible arc length for a directional glyph. Short biological
# intervals keep their real anchor and metadata, but their drawing interval is
# symmetrically widened just enough to retain direction. The shape factory then
# contracts its head and thickness to produce a compressed glyph instead of
# making a few-base interval look like a conventional long feature.
ggchord_arrow_minimum_length <- function(width, head_length = .04,
                                          heads = 1L) {
  if (heads > 1L) {
    max(.70 * head_length, .55 * width)
  } else {
    max(.32 * head_length, .35 * width)
  }
}

ggchord_arrow_display_interval <- function(a_start, a_end, r0, width,
                                            head_length = .04,
                                            heads = 1L) {
  span <- a_end - a_start
  direction <- if (span < 0) -1 else 1
  minimum_angle <- ggchord_arrow_minimum_length(
    width, head_length, heads
  ) / max(abs(r0), .1)
  if (abs(span) >= minimum_angle) return(c(a_start, a_end))
  midpoint <- (a_start + a_end) / 2
  midpoint + direction * c(-.5, .5) * minimum_angle
}

ggchord_feature_intervals <- function(start, end, length, strand,
                                      circular = FALSE, continuous = FALSE) {
  if (isTRUE(circular) && start > end) {
    pieces <- if (isTRUE(continuous)) list(c(start, end + length)) else
      list(c(start, length), c(1, end))
  } else {
    pieces <- list(c(min(start, end), max(start, end)))
  }
  if (identical(strand, "-")) {
    pieces <- rev(lapply(pieces, rev))
  }
  lapply(seq_along(pieces), function(i) {
    list(start = pieces[[i]][1L], end = pieces[[i]][2L],
         draw_head = i == length(pieces))
  })
}

ggchord_feature_angle_interval <- function(piece, seq_length, orientation,
                                            angle_start, angle_end) {
  fraction <- function(x) {
    if (orientation == 1) x / seq_length else 1 - x / seq_length
  }
  c(
    angle_start + fraction(piece$start) * (angle_end - angle_start),
    angle_start + fraction(piece$end) * (angle_end - angle_start)
  )
}

ggchord_shape_block <- function(a_start, a_end, r0, width, n = 60L) {
  ang <- seq(a_start, a_end, length.out = n)
  list(list(
    angle = c(ang, rev(ang)),
    radius = c(rep(r0 + width / 2, n), rep(r0 - width / 2, n))
  ))
}

ggchord_shape_wedge <- function(a_start, a_end, r0, width) {
  list(list(
    angle = c(a_start, a_start, a_end),
    radius = c(r0 + width / 2, r0 - width / 2, r0)
  ))
}

ggchord_shape_arrow <- function(a_start, a_end, r0, width,
                                head_length = 0.04,
                                head_width = 1,
                                short_feature = "auto",
                                draw_head = TRUE) {
  if (!isTRUE(draw_head) || head_length <= 0) {
    return(ggchord_shape_block(a_start, a_end, r0, width))
  }
  span <- a_end - a_start
  direction <- if (span < 0) -1 else 1
  available <- abs(span) * max(abs(r0), 0.1)
  if (available < head_length) {
    fallback <- short_feature
    if (identical(fallback, "auto")) {
      # Preserve the nominal radial band.  A short biological interval may
      # compress its longitudinal body/head, but must not become visually
      # thinner than neighbouring features merely because it is short.
      head_width <- min(head_width, 1)
      head_length <- min(available * .36, width * .45)
      fallback <- "arrow"
    }
    if (identical(fallback, "error"))
      ggchord_stop("geom_feature(): the visible interval is shorter than its arrow head")
    if (identical(fallback, "block")) {
      return(ggchord_shape_block(a_start, a_end, r0, width))
    }
    if (identical(fallback, "wedge")) {
      return(ggchord_shape_wedge(a_start, a_end, r0, width))
    }
  }
  head_angle <- min(abs(span) * 0.45, head_length / max(abs(r0), 0.1))
  shoulder <- a_end - direction * head_angle
  body <- seq(a_start, shoulder, length.out = 35L)
  half_body <- width / 2
  half_head <- half_body * head_width
  list(list(
    angle = c(body, shoulder, a_end, shoulder, rev(body)),
    radius = c(
      rep(r0 + half_body, length(body)), r0 + half_head, r0,
      r0 - half_head, rep(r0 - half_body, length(body))
    )
  ))
}

ggchord_shape_bidirectional_arrow <- function(a_start, a_end, r0, width,
                                               head_length = 0.04,
                                               head_width = 1,
                                               short_feature = "auto") {
  span <- a_end - a_start
  direction <- if (span < 0) -1 else 1
  available <- abs(span) * max(abs(r0), 0.1)
  if (available < 2 * head_length) {
    if (identical(short_feature, "block")) {
      return(ggchord_shape_block(a_start, a_end, r0, width))
    }
    if (identical(short_feature, "error"))
      ggchord_stop("geom_feature(): the visible interval is shorter than its two arrow heads")
    if (identical(short_feature, "auto")) {
      # Preserve two visible directions and the nominal radial band.  Only the
      # longitudinal head/body dimensions are compressed.
      head_width <- min(head_width, 1)
      head_length <- min(available * .22, width * .42)
    } else if (identical(short_feature, "wedge")) {
      mid <- (a_start + a_end) / 2
      return(list(list(
        angle = c(a_start, mid, a_end, mid),
        radius = c(r0, r0 + width / 2, r0, r0 - width / 2)
      )))
    }
  }
  head_angle <- min(abs(span) * .225,
    head_length / max(abs(r0), .1))
  start_shoulder <- a_start + direction * head_angle
  end_shoulder <- a_end - direction * head_angle
  body <- seq(start_shoulder, end_shoulder, length.out = 35L)
  half_body <- width / 2
  half_head <- half_body * head_width
  list(list(
    angle = c(
      a_start, start_shoulder, body, end_shoulder,
      a_end, end_shoulder, rev(body), start_shoulder
    ),
    radius = c(
      r0, r0 + half_head, rep(r0 + half_body, length(body)),
      r0 + half_head, r0, r0 - half_head,
      rep(r0 - half_body, length(body)), r0 - half_head
    )
  ))
}

ggchord_shape_chevron <- function(a_start, a_end, r0, width) {
  span <- a_end - a_start
  shoulder <- seq(a_start, a_start + 0.68 * span, length.out = 30L)
  list(list(
    angle = c(shoulder, a_end, rev(shoulder), a_start + 0.28 * span),
    radius = c(rep(r0 + width / 2, length(shoulder)), r0,
               rep(r0 - width / 2, length(shoulder)), r0)
  ))
}

# Translate a midpoint symbol through the track projection. The head is
# placed in the plane of the graphic so it stays round on a curved reference.
ggchord_interval_anchor_parts <- function(profile, start, end, center, width,
    backbone, ref, along_metric, normal_metric) {
  shape <- profile$name
  mid <- (start + end) / 2
  span <- abs(end - start)
  ratio <- along_metric / normal_metric
  radius <- if (is.null(ref)) 1 else max(abs(center), .1)
  project <- if (is.null(ref)) cbind else function(x, y) map_to_curve_many(x, y, ref)
  part <- function(x, y, role, paint = "full") {
    list(xy = project(x, y), part = role, paint = paint)
  }
  if (shape == "tick") {
    height <- (profile$body_width %||% (width / normal_metric)) * normal_metric
    return(list(part(rep(mid, 2L), center + c(-.5, .5) * height, "body", "path")))
  }
  if (shape == "diamond") {
    height <- (profile$head_width %||% (width / normal_metric)) * normal_metric
    half <- height * ratio * .34 / radius
    if (span > 0) half <- min(span / 2, max(half, span * .22))
    return(list(part(mid + c(-half, 0, half, 0),
      center + c(0, .5, 0, -.5) * height, "head")))
  }
  head_radius <- if (!is.null(profile$head_width)) profile$head_width * normal_metric / 2 else
    width * if (shape == "lollipop") .58 else .5
  xy <- as.numeric(project(mid, center))
  theta <- seq(0, 2 * pi, length.out = 48L)
  head <- list(xy = cbind(xy[1L] + head_radius * (if (is.null(ref)) ratio else 1) * cos(theta),
    xy[2L] + head_radius * sin(theta)), part = "head", paint = "full")
  if (shape == "point") return(list(head))
  half <- if (is.null(profile$body_width)) width * ratio * .10 / radius else
    profile$body_width * along_metric / (2 * radius)
  if (span > 0) half <- min(half, span * if (is.null(profile$body_width)) .08 else .5)
  stem <- part(mid + c(-half, half, half, -half),
    c(backbone, backbone, center, center), "body")
  list(stem, head)
}

ggchord_feature_geometry <- function(shape, a_start, a_end, r0, width,
                                     sequence_radius, ref,
                                     arrow_head_length = 0.04,
                                     arrow_head_width = NULL,
                                     short_feature = "auto",
                                     draw_head = TRUE,
                                     draw_start_head = FALSE,
                                     bidirectional = FALSE,
                                     fixed_body_width = FALSE,
                                     along_width = width, profile = NULL,
                                     nominal_width = width, along_metric = 1,
                                     normal_metric = 1, draw_tail = TRUE) {
  if (isTRUE(profile$components) && identical(shape, profile$name))
    return(ggchord_profile_geometry(profile, a_start, a_end, r0, width,
      nominal_width, along_metric, normal_metric, draw_head, draw_start_head,
      draw_tail, bidirectional, fixed_body_width))
  if (identical(shape, "arrow") && is.null(arrow_head_width)) {
    if (fixed_body_width) {
      arrow_head_width <- ggchord_circular_arrow_head_ratio
    } else {
      # Retain the nominal maximum thickness by narrowing the shaft.
      stem_fraction <- .72
      width <- width * stem_fraction
      arrow_head_width <- 1 / stem_fraction
    }
  }
  arrow_head_width <- arrow_head_width %||% 1
  if (identical(shape, "arrow_triangle")) {
    if (draw_head && draw_start_head && bidirectional) {
      mid <- (a_start + a_end) / 2
      half <- width * arrow_head_width / 2
      return(list(list(angle = c(a_start, mid, a_end, mid),
        radius = c(r0, r0 + half, r0, r0 - half))))
    }
    if (draw_head)
      return(ggchord_shape_wedge(a_start, a_end, r0, width * arrow_head_width))
    if (draw_start_head)
      return(ggchord_shape_wedge(a_end, a_start, r0, width * arrow_head_width))
    return(ggchord_shape_block(a_start, a_end, r0, width))
  }
  if (shape %in% ggchord_midpoint_shapes) {
    profile <- profile %||% ggchord_resolve_feature_profile(shape)
    return(ggchord_interval_anchor_parts(profile, a_start, a_end, r0, width,
      sequence_radius, ref, along_metric, normal_metric))
  }
  if (identical(shape, "capped_line")) {
    # Put both caps inside the displayed interval, including at seams and
    # adjacent features. Leave a little room for the end strokes.
    span <- a_end - a_start
    direction <- sign(span)
    along <- width / max(abs(r0), .1)
    inset <- min(abs(span) * .10, along * .06)
    from <- a_start + direction * inset
    to <- a_end - direction * inset
    cap_length <- direction * min(abs(to - from) / 2,
      along_width * .56 / max(abs(r0), .1))
    line_width <- max(width * .055, .0015)
    body <- ggchord_shape_block(from, to, r0, line_width)
    start_cap <- ggchord_shape_block(
      from, from + cap_length, r0, width * .48
    )
    end_cap <- ggchord_shape_block(
      to - cap_length, to, r0, width * .48
    )
    return(c(body, start_cap, end_cap))
  }
  if (identical(shape, "arrow_arc")) {
    line_width <- max(width * .24, .0015)
    if (!draw_head && !draw_start_head)
      return(ggchord_shape_block(a_start, a_end, r0, line_width))

    head_length <- min(abs(a_end - a_start) * max(abs(r0), .1) * .28,
      width * .55)
    head_width <- width * .62 / line_width
    # One outline joins shaft and head at the shoulders. The shaft never
    # continues beneath the tapered sides of the head.
    if (bidirectional && draw_head && draw_start_head)
      return(ggchord_shape_bidirectional_arrow(a_start, a_end, r0, line_width,
        head_length = head_length, head_width = head_width, short_feature = "shrink"))
    if (draw_start_head && !draw_head)
      return(ggchord_shape_arrow(a_end, a_start, r0, line_width,
        head_length = head_length, head_width = head_width, short_feature = "shrink"))
    return(ggchord_shape_arrow(a_start, a_end, r0, line_width,
      head_length = head_length, head_width = head_width, short_feature = "shrink"))
  }
  if (isTRUE(bidirectional) &&
      isTRUE(draw_head) && isTRUE(draw_start_head)) {
    interval <- if (identical(short_feature, "auto")) {
      ggchord_arrow_display_interval(
        a_start, a_end, r0, width, arrow_head_length, heads = 2L
      )
    } else c(a_start, a_end)
    return(ggchord_shape_bidirectional_arrow(
      interval[1L], interval[2L], r0, width,
      head_length = arrow_head_length,
      head_width = arrow_head_width,
      short_feature = short_feature
    ))
  }
  switch(
    shape,
    block = ggchord_shape_block(a_start, a_end, r0, width),
    chevron = ggchord_shape_chevron(a_start, a_end, r0, width),
    if (isTRUE(draw_start_head) && !isTRUE(draw_head)) {
      interval <- if (identical(short_feature, "auto")) {
        ggchord_arrow_display_interval(
          a_end, a_start, r0, width, arrow_head_length, heads = 1L
        )
      } else c(a_end, a_start)
      ggchord_shape_arrow(
        interval[1L], interval[2L], r0, width,
        head_length = arrow_head_length,
        head_width = arrow_head_width,
        short_feature = short_feature, draw_head = TRUE
      )
    } else {
      interval <- if (identical(short_feature, "auto") && isTRUE(draw_head)) {
        ggchord_arrow_display_interval(
          a_start, a_end, r0, width, arrow_head_length, heads = 1L
        )
      } else c(a_start, a_end)
      ggchord_shape_arrow(
        interval[1L], interval[2L], r0, width,
        head_length = arrow_head_length,
        head_width = arrow_head_width,
        short_feature = short_feature,
        draw_head = draw_head
      )
    }
  )
}
