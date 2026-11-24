ggchord_anchor_mark <- function(profile, location, pos, boundaries, size) {
  shape <- profile$name
  glyph <- function(pos, normal, component, part, side = NA_character_, width = size) {
    list(pos = pos, normal = normal, component = component, part = part,
      side = side, width = width)
  }
  local <- function(pos, direction = 1, side = NA_character_) {
    switch(shape,
      block = {
        width <- profile$body_width %||% size
        list(glyph(pos, if (is.na(side)) 0 else direction * width / 2,
          "block", "body", side, width))
      },
      point = list(glyph(pos, 0, "point", "head", side, profile$head_width %||% size)),
      diamond = list(glyph(pos, 0, "diamond", "head", side, profile$head_width %||% size)),
      lollipop = {
        tip <- direction * 1.5 * size
        list(glyph(pos, c(0, tip), "path", "body", side),
          glyph(pos, tip, "point", "head", side, profile$head_width %||% size))
      },
      tick = list(glyph(pos, c(-.5, .5) * (profile$body_width %||% size), "path", "body")),
      cut = list(glyph(pos, sort(c(0, direction * (profile$body_width %||% size) / 2)),
        "path", "body", side)))
  }
  if (shape == "cut" && all(is.na(boundaries))) {
    return(list(glyph(pos, c(-.5, .5) * (profile$body_width %||% size), "path", "body")))
  }
  use_boundaries <- shape == "cut" ||
    (location == "cut" && shape %in% c("block", "diamond", "point", "lollipop"))
  if (!use_boundaries) return(local(pos))
  available <- which(!is.na(boundaries))
  unlist(lapply(available, function(i) {
    local(boundaries[i], c(1, -1)[i], c("top", "bottom")[i])
  }), recursive = FALSE)
}

# Resolve a configurable interval into a silhouette and, when needed,
# separately painted parts. Lengths are projected along the reference and
# widths are projected across it before the profile is drawn.
ggchord_profile_geometry <- function(profile, a_start, a_end, r0, width,
    nominal_width, along_metric, normal_metric, draw_head, draw_start_head,
    draw_tail, bidirectional, fixed_body_width) {
  span <- abs(a_end - a_start)
  radius <- max(abs(r0), .1)
  ends <- c(tail = profile$tail_shape, head = profile$head_shape)
  roles <- names(ends)
  roles <- setNames(roles, roles)
  if (!draw_tail) ends["tail"] <- "none"
  if (draw_start_head && profile$head_shape == "triangle" && profile$tail_shape == "none") {
    ends["tail"] <- profile$head_shape
    roles["tail"] <- "head"
  }
  if (!draw_head) ends["head"] <- "none"
  as_block <- function() {
    polygon <- ggchord_shape_block(a_start, a_end, r0, width)[[1L]]
    polygon$part <- "body"
    polygon$paint <- "full"
    list(polygon)
  }
  if (!span || all(ends == "none") && profile$body_shape != "line")
    return(as_block())

  body_width <- if (profile$body_shape == "line" && is.null(profile$body_width)) {
    nominal_width * (if (profile$name == "capped_line") .055 else .24)
  } else width
  if (profile$name == "arrow" && !fixed_body_width &&
      is.null(profile$body_width) && is.null(profile$head_width))
    body_width <- width * .72
  body_width <- max(body_width, .Machine$double.eps)

  shape_length <- function(side) {
    kind <- ends[[side]]
    if (kind == "none") return(0)
    if (profile$body_shape == "none")
      return(if (bidirectional && all(ends == "triangle")) span / 2 else span)
    requested <- profile[[paste0(roles[[side]], "_length")]]
    if (!is.null(requested)) return(requested * along_metric / radius)
    if (kind == "cap") return(nominal_width * .56 / radius)
    if (kind == "notch") return(span * .28)
    if (profile$name == "chevron") return(span * .32)
    if (profile$name == "arrow_arc") return(min(span * .28, nominal_width * .55 / radius))
    2 * along_metric / radius
  }
  shape_width <- function(side) {
    kind <- ends[[side]]
    if (kind %in% c("none", "notch") || profile$name == "arrow_flush") return(body_width)
    requested <- profile[[paste0(roles[[side]], "_width")]]
    if (!is.null(requested)) return(requested * normal_metric)
    if (kind == "cap") return(max(nominal_width * .48, body_width * 1.3))
    if (profile$name == "arrow_arc") return(nominal_width * .62)
    if (profile$name %in% c("arrow_triangle", "chevron")) return(width)
    body_width * ggchord_circular_arrow_head_ratio
  }
  lengths <- vapply(names(ends), shape_length, numeric(1))
  widths <- vapply(names(ends), shape_width, numeric(1))
  total <- sum(lengths)
  if (profile$body_shape != "none" && total >= span && total > 0) {
    if (profile$short_feature == "error")
      ggchord_stop("geom_feature(): visible interval is too short for the configured ends")
    if (profile$short_feature == "block") return(as_block())
    lengths <- lengths * (.9 * span / total)
  }
  if (profile$name %in% c("arrow", "arrow_flush")) {
    lengths <- pmin(lengths, span * if (bidirectional) .225 else .45)
  }
  limits <- c(0, lengths[["tail"]], span - lengths[["head"]], span)
  if (profile$body_shape == "none") limits[3L] <- limits[2L]
  wrap <- function(xy, part = "body", paint = "full") {
    list(angle = a_start + sign(a_end - a_start) * xy[, 1L],
      radius = r0 + xy[, 2L], part = part, paint = paint)
  }
  node <- function(i, kind, size, part) {
    a <- limits[i]
    b <- limits[i + 1L]
    if (a >= b) return(NULL)
    x <- seq(a, b, length.out = 30L)
    upper <- rep(size / 2, length(x))
    if (kind == "triangle")
      upper <- upper * if (i == 1L) seq(0, 1, length.out = length(x)) else
        seq(1, 0, length.out = length(x))
    list(upper = cbind(x, upper), lower = cbind(x, -upper), part = part,
      notch = if (kind == "notch") c(if (i == 1L) b else a, 0))
  }
  nodes <- Filter(Negate(is.null), list(
    node(1L, ends[["tail"]], widths[["tail"]], roles[["tail"]]),
    node(2L, "body", body_width, "body"),
    node(3L, ends[["head"]], widths[["head"]], roles[["head"]])))
  reverse_rows <- function(x) x[rev(seq_len(nrow(x))), , drop = FALSE]
  first <- nodes[[1L]]
  last <- nodes[[length(nodes)]]
  top <- do.call(rbind, lapply(nodes, `[[`, "upper"))
  bottom <- do.call(rbind, lapply(nodes, `[[`, "lower"))
  outline <- rbind(top, last$notch, reverse_rows(bottom), first$notch)

  styled <- any(grepl("_(fill|colour)$", names(profile)))
  if (!styled) return(list(wrap(outline)))
  fill <- lapply(nodes, function(n) {
    at_start <- identical(n, first)
    at_end <- identical(n, last)
    wrap(rbind(n$upper, if (at_end) n$notch, reverse_rows(n$lower),
      if (at_start) n$notch), n$part, "fill")
  })
  if (!any(grepl("_colour$", names(profile))))
    return(c(fill, list(wrap(outline, paint = "outline"))))

  # When outlines differ by part, trace only the exterior edges. Each
  # shoulder connects to the neighbouring part without drawing a seam.
  edges <- lapply(seq_along(nodes), function(i) {
    n <- nodes[[i]]
    lapply(c("upper", "lower"), function(edge) {
      points <- n[[edge]]
      if (i > 1L) points <- rbind(tail(nodes[[i - 1L]][[edge]], 1L), points)
      wrap(points, n$part, "path")
    })
  })
  first_edge <- rbind(head(first$lower, 1L), first$notch, head(first$upper, 1L))
  last_edge <- rbind(tail(last$upper, 1L), last$notch, tail(last$lower, 1L))
  c(fill, unlist(edges, recursive = FALSE),
    list(wrap(first_edge, first$part, "path"), wrap(last_edge, last$part, "path")))
}
