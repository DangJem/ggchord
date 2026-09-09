# Shared interval and shape helpers for geom_gene() and geom_feature().

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
    terminator = "regulatory", rbs = "regulatory"
  )
  out <- unname(aliases[key])
  out[is.na(out) | !nzchar(out)] <- "misc_feature"
  out
}

ggchord_feature_class_shape <- function(feature_class) {
  values <- c(
    gene = "arrow", resistance_gene = "arrow", origin = "arrow",
    promoter = "promoter_arrow", primer = "primer_arrow",
    binding_site = "block", operator = "block", regulatory = "compact_arrow",
    mcs = "block", misc_feature = "block"
  )
  out <- unname(values[as.character(feature_class)])
  out[is.na(out)] <- "block"
  out
}

ggchord_feature_width_factor <- function(shape) {
  values <- c(
    arrow = 1, compact_arrow = .72, promoter_arrow = .55,
    primer_arrow = .48, marker = .38, block = .72,
    chevron = .68, lollipop = .55
  )
  out <- unname(values[as.character(shape)])
  out[is.na(out)] <- 1
  out
}

ggchord_feature_preferred_lane <- function(feature_class) {
  values <- c(
    gene = 0L, resistance_gene = 0L, origin = 0L,
    mcs = 1L, binding_site = 1L, operator = 1L,
    regulatory = 1L, promoter = 2L, primer = 2L,
    misc_feature = 1L
  )
  out <- unname(values[as.character(feature_class)])
  out[is.na(out)] <- 0L
  as.integer(out)
}

ggchord_feature_intervals <- function(start, end, length, strand,
                                      circular = FALSE) {
  if (isTRUE(circular) && start > end) {
    pieces <- list(c(start, length), c(1, end))
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

ggchord_shape_marker <- function(a_start, a_end, r0, width) {
  mid <- (a_start + a_end) / 2
  half <- min(abs(a_end - a_start) / 2,
    width * .34 / max(abs(r0), .1))
  direction <- if (a_end < a_start) -1 else 1
  half <- direction * max(half, abs(a_end - a_start) * .22)
  list(list(
    angle = c(mid - half, mid, mid + half, mid),
    radius = c(r0, r0 + width / 2, r0, r0 - width / 2)
  ))
}

ggchord_minimum_display_span <- function(a_start, a_end, r0, width,
                                          multiple) {
  span <- a_end - a_start
  direction <- if (span < 0) -1 else 1
  visual_span <- max(abs(span), width * multiple / max(abs(r0), .1))
  mid <- (a_start + a_end) / 2
  c(mid - direction * visual_span / 2, mid + direction * visual_span / 2)
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
      # Preserve a visible rectangular body and shrink only the head. A full
      # wedge makes short plasmid features look like detached triangles.
      if (available < width * .18) {
        return(ggchord_shape_block(a_start, a_end, r0, width))
      }
      head_length <- available * .36
      fallback <- "arrow"
    }
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
    if (identical(short_feature, "block") || identical(short_feature, "auto")) {
      return(ggchord_shape_block(a_start, a_end, r0, width))
    }
    mid <- (a_start + a_end) / 2
    return(list(list(
      angle = c(a_start, mid, a_end, mid),
      radius = c(r0, r0 + width / 2, r0, r0 - width / 2)
    )))
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

ggchord_shape_lollipop <- function(a_start, a_end, r0, width,
                                   sequence_radius, ref) {
  span <- a_end - a_start
  mid <- (a_start + a_end) / 2
  angle_half <- min(abs(span) * 0.08, width * 0.10 / max(abs(r0), 0.1))
  stem <- list(
    angle = c(mid - angle_half, mid + angle_half,
              mid + angle_half, mid - angle_half),
    radius = c(sequence_radius, sequence_radius, r0, r0)
  )
  theta <- seq(0, 2 * pi, length.out = 48L)
  head_radius <- width * 0.58
  center <- as.numeric(map_to_curve_many(mid, r0, ref)[1L, ])
  delta <- max(abs(span) * 1e-4, 1e-7)
  tangent_pts <- map_to_curve_many(c(mid - delta, mid + delta), rep(r0, 2L), ref)
  tangent <- as.numeric(tangent_pts[2L, ] - tangent_pts[1L, ])
  tangent_norm <- sqrt(sum(tangent^2))
  if (!is.finite(tangent_norm) || tangent_norm <= 1e-12) tangent <- c(1, 0)
  else tangent <- tangent / tangent_norm
  normal <- c(-tangent[2L], tangent[1L])
  head <- list(xy = cbind(
    center[1L] + head_radius *
      (cos(theta) * tangent[1L] + sin(theta) * normal[1L]),
    center[2L] + head_radius *
      (cos(theta) * tangent[2L] + sin(theta) * normal[2L])
  ))
  list(stem, head)
}

ggchord_feature_geometry <- function(shape, a_start, a_end, r0, width,
                                     sequence_radius, ref,
                                     arrow_head_length = 0.04,
                                     arrow_head_width = 1,
                                     arrow_head_style = "shouldered",
                                     short_feature = "auto",
                                     draw_head = TRUE,
                                     bidirectional = FALSE) {
  if (shape %in% c("compact_arrow", "promoter_arrow", "primer_arrow")) {
    available <- abs(a_end - a_start) * max(abs(r0), .1)
    if (available < width * .42) shape <- "marker"
  }
  if (identical(shape, "marker")) {
    return(ggchord_shape_marker(a_start, a_end, r0, width))
  }
  if (identical(shape, "compact_arrow")) {
    return(ggchord_shape_arrow(
      a_start, a_end, r0, width,
      head_length = min(arrow_head_length * .78,
        abs(a_end - a_start) * max(abs(r0), .1) * .38),
      head_width = min(arrow_head_width, 1.18),
      short_feature = short_feature, draw_head = draw_head
    ))
  }
  if (identical(shape, "promoter_arrow")) {
    display <- ggchord_minimum_display_span(
      a_start, a_end, r0, width, multiple = 2.6
    )
    return(ggchord_shape_arrow(
      display[1L], display[2L], r0, width,
      head_length = min(arrow_head_length * .72,
        abs(display[2L] - display[1L]) * max(abs(r0), .1) * .34),
      head_width = 1, short_feature = "block", draw_head = draw_head
        ))
  }
  if (identical(shape, "primer_arrow")) {
    display <- ggchord_minimum_display_span(
      a_start, a_end, r0, width, multiple = 3.0
    )
    return(ggchord_shape_arrow(
      display[1L], display[2L], r0, width,
      head_length = min(arrow_head_length * .60,
        abs(display[2L] - display[1L]) * max(abs(r0), .1) * .42),
      head_width = 1.08, short_feature = "block", draw_head = draw_head
    ))
  }
  if (identical(shape, "arrow") && isTRUE(bidirectional) &&
      isTRUE(draw_head)) {
    return(ggchord_shape_bidirectional_arrow(
      a_start, a_end, r0, width,
      head_length = arrow_head_length,
      head_width = if (identical(arrow_head_style, "flush")) 1 else
        arrow_head_width,
      short_feature = short_feature
    ))
  }
  switch(
    shape,
    block = ggchord_shape_block(a_start, a_end, r0, width),
    chevron = ggchord_shape_chevron(a_start, a_end, r0, width),
    lollipop = ggchord_shape_lollipop(
      a_start, a_end, r0, width, sequence_radius, ref
    ),
    if (identical(arrow_head_style, "triangle") && isTRUE(draw_head)) {
      ggchord_shape_wedge(a_start, a_end, r0, width * arrow_head_width)
    } else {
      ggchord_shape_arrow(
        a_start, a_end, r0, width,
        head_length = arrow_head_length,
        head_width = if (identical(arrow_head_style, "flush")) 1 else
          arrow_head_width,
        short_feature = short_feature,
        draw_head = draw_head
      )
    }
  )
}
