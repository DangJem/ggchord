# Boxes describe a centre and extent in fractions of their parent key.
ggchord_legend_slot <- function(glyph, box, scale = 1, clip = "off") {
  grobTree(glyph, vp = viewport(x = box[1], y = box[2], width = box[3], height = box[4],
    clip = clip, gp = gpar(lex = scale, cex = scale)))
}

ggchord_legend_slots <- function(n, size) {
  columns <- seq_len(n)
  rows <- ceiling(n / columns)
  score <- pmin(size[1] / columns / 1.9, size[2] / rows / 1.65)
  at <- which.max(score)
  ncol <- columns[at]
  nrow <- rows[at]
  index <- seq_len(n) - 1L
  row <- index %/% ncol
  columns_in_row <- pmin(ncol, n - row * ncol)
  cbind(x = .5 + (index %% ncol - (columns_in_row - 1) / 2) / ncol,
    y = 1 - (row + .5) / nrow, width = .90 / ncol, height = .90 / nrow)
}

ggchord_legend_plan <- function(parts, size) {
  if (!length(parts)) return(list(groups = list(), crowded = FALSE))
  ids <- vapply(seq_along(parts), function(i) {
    part <- parts[[i]]
    if (part$builtin && length(part$source) == 1L) part$source else -i
  }, integer(1))
  groups <- split(parts, match(ids, unique(ids)))
  outer <- if (length(groups) == 1L) matrix(c(.5, .5, 1, 1), nrow = 1L) else
    ggchord_legend_slots(length(groups), size)
  layouts <- lapply(seq_along(groups), function(i) {
    group <- groups[[i]]
    box <- outer[i, ]
    room <- size * box[3:4]
    n <- length(group)
    roles <- vapply(group, `[[`, character(1), "role")
    # The backbone can anchor several relations in one corridor. Dense
    # annotation and Site combinations need more room than that template offers.
    annotations <- sum(roles %in% c("gene", "feature", "primer"))
    annotation_capacity <- if ("site" %in% roles) 1L else 2L
    capacity <- if ("seq" %in% roles && annotations <= annotation_capacity) 5L else 4L
    compact <- n >= 2L && n <= capacity && !anyDuplicated(roles) &&
      all(room >= c(3.3, 2.65) + (n - 2L) * c(.62, .45))
    slots <- if (compact || n == 1L) NULL else ggchord_legend_slots(n, room)
    single_room <- if (is.null(slots)) room else room * slots[1, 3:4]
    crowded <- any(single_room < c(1.45, 1.25)) ||
      any(vapply(group, function(part) part$builtin && .52 * ggchord_legend_mark_mm(part) >
        min(single_room), logical(1)))
    list(parts = group, box = box, room = room, compact = compact,
      slots = slots, crowded = crowded)
  })
  list(groups = layouts, crowded = any(vapply(layouts, `[[`, logical(1), "crowded")))
}

ggchord_legend_mark_mm <- function(part) {
  data <- part$style
  diameter <- if (inherits(part$glyph, "points")) (data$size %||% 0) + (data$stroke %||% 0) else 0
  arrow <- part$params$arrow
  arrow_mm <- if (inherits(arrow, "arrow")) convertWidth(arrow$length[1L], "mm", TRUE) else 0
  max(0, diameter, data$linewidth, arrow_mm, na.rm = TRUE)
}

ggchord_legend_scale_arrow <- function(glyph, scale) {
  if (inherits(glyph$arrow, "arrow")) glyph$arrow$length <- glyph$arrow$length * scale
  if (inherits(glyph, "gTree")) {
    for (i in seq_along(glyph$children))
      glyph$children[[i]] <- ggchord_legend_scale_arrow(glyph$children[[i]], scale)
  }
  glyph
}

ggchord_legend_single <- function(part, box, room) {
  if (part$builtin) {
    mark <- max(.01, ggchord_legend_mark_mm(part))
    scale <- min(1, room / c(4.1, 2.9), min(room) / (1.15 * mark))
    glyph <- if (part$role == "link")
      ggchord_legend_connection(c(.19, .74), c(.81, .26), part, NULL)$line else part$glyph
    return(ggchord_legend_slot(ggchord_legend_scale_arrow(glyph, scale), box, scale, clip = "on"))
  }
  glyph <- part$draw(part$data, part$params, room)
  extent <- 10 * c(attr(glyph, "width") %||% 0, attr(glyph, "height") %||% 0)
  out <- ggchord_legend_slot(glyph, box, clip = "on")
  attr(out, "ggchord_key_crowded") <- any(extent > room)
  out
}

ggchord_legend_render <- function(layout) {
  groups <- lapply(layout$groups, function(layout) {
    box <- layout$box
    room <- layout$room
    parts <- layout$parts
    if (layout$compact)
      return(ggchord_legend_slot(ggchord_contact_key(parts, room), box))
    if (is.null(layout$slots)) return(ggchord_legend_single(parts[[1L]], box, room))
    glyphs <- lapply(seq_along(parts), function(i) {
      slot <- layout$slots[i, ]
      ggchord_legend_single(parts[[i]], slot, room * slot[3:4])
    })
    ggchord_legend_slot(do.call(grobTree, glyphs), box)
  })
  do.call(grobTree, groups)
}

# Use the outline of the key symbol itself to choose a diagram endpoint.
ggchord_legend_feature_anchor <- function(part, box, side, room) {
  glyph <- part$glyph
  paths <- if (inherits(glyph, "gTree")) glyph$children else list(glyph)
  polygons <- Filter(function(path) inherits(path, "polygon"), paths)
  if (!length(polygons)) {
    if (side != "top") return(box[1:2] + c(0, -.31 * box[4]))
    mark <- Filter(function(path) inherits(path, "circle"), paths)[[1L]]
    top <- (as.numeric(mark$y) - .5) * box[4] +
      as.numeric(mark$r) * min(room * box[3:4]) / room[2]
    return(box[1:2] + c(0, top))
  }
  x <- unlist(lapply(polygons, function(path) as.numeric(path$x)))
  y <- unlist(lapply(polygons, function(path) as.numeric(path$y)))
  at <- switch(side, left = x == min(x), right = x == max(x),
    top = y == max(y) & abs(x - .5) == min(abs(x[y == max(y)] - .5)))
  box[1:2] + (c(mean(x[at]), mean(y[at])) - .5) * box[3:4]
}

ggchord_legend_site <- function(part, anchor, top, scale) {
  glyph <- part$glyph
  if (inherits(glyph, "points")) {
    style <- part$style
    stem <- segmentsGrob(anchor[1], anchor[2], anchor[1], top,
      gp = gpar(col = glyph$gp$col, lwd = (style$linewidth %||% .4) * ggplot2::.pt,
        lty = style$linetype %||% 1))
    mark <- ggchord_legend_slot(glyph, c(anchor[1], top, 1, 1))
    return(ggchord_legend_slot(grobTree(stem, mark), c(.5, .5, 1, 1), scale))
  }
  if (inherits(glyph, "polyline")) {
    x <- anchor[1] + c(-.14, -.14, 0, 0)
    y <- c(top, mean(c(anchor[2], top)), mean(c(anchor[2], top)), anchor[2])
    return(grid::polylineGrob(x, y, gp = glyph$gp, default.units = "npc"))
  }
  segmentsGrob(anchor[1], anchor[2], anchor[1], top, gp = glyph$gp)
}

ggchord_legend_connection <- function(start, end, link, ribbon) {
  t <- seq(0, 1, length.out = 25L)
  route <- if (!is.null(ribbon)) ribbon$route else link$route
  along <- if (identical(route, "straight")) t else t^2 * (3 - 2 * t)
  x <- start[1] + (end[1] - start[1]) * along
  y <- start[2] + (end[2] - start[2]) * t
  band <- line <- NULL
  if (!is.null(ribbon)) {
    half <- .095 - .028 * sin(pi * t)
    tip <- ribbon$params$tip %||% "none"
    if (tip %in% c("from", "both")) half <- half * pmin(1, t / .20)
    if (tip %in% c("to", "both")) half <- half * pmin(1, (1 - t) / .20)
    gp <- ribbon$glyph$gp
    if (inherits(gp$fill, "GridLinearGradient")) {
      gradient <- gp$fill
      gradient$x1 <- unit(if (end[1] >= start[1]) 0 else 1, "npc")
      gradient$x2 <- unit(if (end[1] >= start[1]) 1 else 0, "npc")
      gradient$y1 <- unit(1, "npc")
      gradient$y2 <- unit(0, "npc")
      gp$fill <- gradient
    }
    band <- polygonGrob(unit(c(x - half, rev(x + half)), "npc"),
      unit(c(y, rev(y)), "npc"), gp = gp)
  }
  if (!is.null(link)) line <- grid::polylineGrob(x, y, default.units = "npc",
    arrow = link$glyph$arrow, gp = link$glyph$gp)
  list(band = band, line = line)
}

ggchord_contact_key <- function(parts, room) {
  names(parts) <- vapply(parts, `[[`, character(1), "role")
  annotations <- intersect(c("gene", "feature", "primer"), names(parts))
  seq_part <- parts$seq
  site <- parts$site
  relationship <- any(c("link", "ribbon") %in% names(parts))
  n <- length(annotations)
  seq_y <- if (relationship) if (n) .51 else .69 else .36
  seq_box <- c(.5, seq_y, .94, .23)
  seq_gap <- if (identical(seq_part$params$seq_key_style, "double")) .15 * seq_box[4] else 0
  scale <- min(1, room / c(6.5, 5.2))
  symbols <- list()
  anchor <- c(.52, .80)
  if (!is.null(seq_part)) {
    symbols$seq <- ggchord_legend_slot(ggchord_legend_scale_arrow(seq_part$glyph, scale), seq_box)
    anchor <- c(.52, seq_y - seq_gap)
  }
  if (n) {
    y <- if (!is.null(seq_part)) {
      clearance <- switch(seq_part$params$seq_key_style %||% "", band = .06, double = .03, 0)
      (if (relationship) .79 else .63) + clearance
    } else if (relationship) {
      if (is.null(site)) .72 else .56
    } else if (!is.null(site)) .48 else .52
    width <- if (n == 1L) .59 else (if (is.null(site)) .88 else .64) / n
    center <- if (is.null(seq_part) && relationship) {
      if (identical(as.character(parts[[annotations[1L]]]$style$strand), "-")) .63 else .37
    } else if (!is.null(seq_part) && !is.null(site)) .38 else .5
    # A block can frame a gene when both are present. For other feature
    # shapes, adjacent symbols leave both outlines and directions visible.
    nested <- identical(annotations, c("gene", "feature")) &&
      identical(parts$feature$params$arrow_key$shape, "block") &&
      !any(as.character(parts$feature$style$shape) %in% ggchord_feature_shapes)
    boxes <- lapply(seq_along(annotations), function(i)
      c(if (nested) center else center + width * (i - (n + 1) / 2), y,
        if (nested) c(.43, .74)[i] else width, if (nested) c(.22, .46)[i] else .34))
    names(boxes) <- annotations
    for (role in if (nested) rev(annotations) else annotations)
      symbols[[role]] <- ggchord_legend_slot(parts[[role]]$glyph, boxes[[role]])
    if (is.null(seq_part)) {
      role <- if (nested) "feature" else annotations[1L]
      side <- if (identical(as.character(parts[[role]]$style$strand), "-")) "left" else "right"
      anchor <- ggchord_legend_feature_anchor(parts[[role]], boxes[[role]], side, room)
    }
  }
  if (!is.null(site)) {
    foot <- if (!is.null(seq_part)) c(if (n) .79 else anchor[1], seq_y + seq_gap) else
      if (n) ggchord_legend_feature_anchor(parts[[annotations[n]]], boxes[[n]], "top", room) else
        c(.37, .63)
    top <- min(.86, foot[2] + .30)
    symbols$site <- ggchord_legend_site(site, foot, top, min(.65, scale))
    if (is.null(seq_part) && !n) anchor <- foot
  }
  end <- c(if (anchor[1] < .34) .17 else .83, .14)
  relations <- if (relationship) ggchord_legend_connection(anchor, end, parts$link, parts$ribbon)
  if (!is.null(relations$line)) relations$line <- ggchord_legend_scale_arrow(relations$line, scale)
  # The skeleton and feature outlines remain in the foreground at contact.
  glyphs <- Filter(Negate(is.null), c(relations, symbols))
  ggchord_legend_slot(do.call(grobTree, glyphs), c(.5, .5, 1, 1), scale)
}
