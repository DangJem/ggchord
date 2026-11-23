# The V3 artwork uses a 65 by 52 coordinate system with y increasing downwards.
ggchord_legend_artwork_mm <- c(6.5, 5.2)

ggchord_legend_x <- function(x) unit(x / 65, "npc")
ggchord_legend_y <- function(y) unit(1 - y / 52, "npc")

ggchord_legend_slot <- function(glyph, box, scale = 1, clip = "inherit") {
  grobTree(glyph, vp = viewport(x = box[1], y = box[2], width = box[3], height = box[4],
    clip = clip, gp = gpar(lex = scale, cex = scale)))
}

ggchord_legend_slots <- function(n, size) {
  columns <- seq_len(n)
  rows <- ceiling(n / columns)
  at <- which.max(pmin(size[1] / columns / 1.9, size[2] / rows / 1.65))
  ncol <- columns[at]
  nrow <- rows[at]
  index <- seq_len(n) - 1L
  row <- index %/% ncol
  columns_in_row <- pmin(ncol, n - row * ncol)
  cbind(x = .5 + (index %% ncol - (columns_in_row - 1) / 2) / ncol,
    y = 1 - (row + .5) / nrow, width = .90 / ncol, height = .90 / nrow)
}

ggchord_legend_groups <- function(parts) {
  sources <- lapply(parts, function(part) if (part$builtin) part$source)
  categories <- vapply(parts, function(part) part$category %||% "", character(1))
  # Keep a genuinely shared scale together when its layers map different
  # columns. It can join another scale only if its category agrees
  # across all the participating glyphs.
  by_source <- split(rep(categories, lengths(sources)), unlist(sources, use.names = FALSE))
  compatible <- vapply(by_source, function(keys) nzchar(keys[1L]) &&
    all(keys == keys[1L]), logical(1))
  ids <- vapply(seq_along(parts), function(i) {
    source <- sources[[i]]
    if (nzchar(categories[i]) && all(compatible[as.character(source)]))
      return(paste0("category:", categories[i]))
    if (length(source) == 1L) as.character(source) else paste0("part:", i)
  }, character(1))
  # A fixed glyph admitted by show.legend has no scale of its own. Place
  # it in the template only when the key offers one unambiguous mapped group.
  mapped <- unique(ids[lengths(sources) > 0L])
  if (length(mapped) == 1L && !startsWith(mapped, "part:")) {
    fixed <- vapply(parts, function(part) isTRUE(part$fixed_key), logical(1))
    ids[fixed] <- mapped
  }
  split(parts, match(ids, unique(ids)))
}

ggchord_legend_plan <- function(parts, size) {
  if (!length(parts)) return(list(groups = list(), crowded = FALSE))
  groups <- ggchord_legend_groups(parts)
  outer <- if (length(groups) == 1L) matrix(c(.5, .5, 1, 1), nrow = 1L) else
    ggchord_legend_slots(length(groups), size)
  layouts <- lapply(seq_along(groups), function(i) {
    group <- groups[[i]]
    box <- outer[i, ]
    room <- size * box[3:4]
    n <- length(group)
    roles <- vapply(group, `[[`, character(1), "role")
    marks <- vapply(group, ggchord_legend_mark_mm, numeric(1))
    # Check the space used by the artwork and unusually large style variants.
    # A full set of seven distinct roles has a compact template in the atlas.
    needed <- c(3.35, 2.65) + (n - 2L) * c(.51, .43)
    max_mark <- if (n >= 5L) 3.6 else 4.2
    compact <- n >= 2L && !anyDuplicated(roles) && all(room >= needed) && max(marks) <= max_mark
    stagger <- n == 2L && all(roles %in% c("gene", "feature", "primer")) &&
      length(unique(roles)) == 1L && all(room >= c(4.9, 4.1))
    slots <- if (compact || n == 1L) NULL else if (stagger) {
      rbind(c(.35, .65, .69, .46), c(.62, .31, .69, .46))
    } else ggchord_legend_slots(n, room)
    single_room <- if (is.null(slots)) room else room * slots[1L, 3:4]
    crowded <- any(single_room < c(1.45, 1.25)) || max(marks) * .52 > min(single_room)
    list(parts = group, box = box, room = room, compact = compact, separated = length(groups) > 1L,
      slots = slots, crowded = crowded)
  })
  list(groups = layouts, crowded = any(vapply(layouts, `[[`, logical(1), "crowded")))
}

ggchord_legend_mark_mm <- function(part) {
  if (!part$builtin) return(0)
  style <- part$style
  marker <- if (part$role == "site" && (part$site_marker || inherits(part$glyph, "points")))
    (style$size %||% 2) + (style$stroke %||% .4) else 0
  arrow <- part$params$arrow
  arrow_mm <- if (inherits(arrow, "arrow")) convertWidth(arrow$length[1L], "mm", TRUE) else 0
  max(0, marker, style$linewidth, arrow_mm, na.rm = TRUE)
}

ggchord_legend_single <- function(part, box, room, separated) {
  if (part$builtin)
    return(ggchord_legend_slot(ggchord_contact_key(list(part), room, separated), box, clip = "on"))
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
    if (is.null(layout$slots))
      return(ggchord_legend_single(parts[[1L]], box, room, layout$separated))
    glyphs <- lapply(seq_along(parts), function(i) {
      slot <- layout$slots[i, ]
      ggchord_legend_single(parts[[i]], slot, room * slot[3:4], separated = TRUE)
    })
    ggchord_legend_slot(do.call(grobTree, glyphs), box)
  })
  do.call(grobTree, groups)
}

# Match the reference stroke at the geom default and retain the relative
# differences between trained or fixed line widths.
ggchord_legend_pen <- function(part, gp, width, default) {
  if (is.null(gp$lwd)) return(gp)
  linewidth <- part$style$linewidth %||% default
  gp$lwd <- if (is.na(linewidth)) 0 else width * (96 / 25.4 / 10) * linewidth / default
  gp
}

ggchord_legend_translucent <- function(colour, opacity) {
  if (is.null(colour)) return(colour)
  if (inherits(colour, "GridLinearGradient")) {
    colour$colours <- ggchord_legend_translucent(colour$colours, opacity)
    return(colour)
  }
  if (!is.character(colour)) return(colour)
  scales::alpha(colour, grDevices::col2rgb(colour, alpha = TRUE)[4L, ] / 255 * opacity)
}

ggchord_legend_polygon <- function(xy, gp) {
  polygonGrob(ggchord_legend_x(xy[, 1L]), ggchord_legend_y(xy[, 2L]), gp = gp)
}

ggchord_legend_curve <- function(points, straight = FALSE) {
  t <- seq(0, 1, length.out = 25L)
  if (straight) return(cbind(
    x = (1 - t) * points[1L, 1L] + t * points[4L, 1L],
    y = (1 - t) * points[1L, 2L] + t * points[4L, 2L]))
  cbind((1 - t)^3, 3 * (1 - t)^2 * t, 3 * (1 - t) * t^2, t^3) %*% points
}

ggchord_legend_seq_artist <- function(part, y, scale) {
  glyph <- part$glyph
  if (inherits(glyph, "rect")) {
    gp <- ggchord_legend_pen(part, glyph$gp, .85, .9)
    grob <- grid::roundrectGrob(ggchord_legend_x(32.5), ggchord_legend_y(y),
      ggchord_legend_x(49), unit(3.6 / 52, "npc"), r = unit(.08, "mm") * scale, gp = gp)
    return(list(grob = grob, upper = y, lower = y))
  }
  lines <- if (inherits(glyph, "gTree")) glyph$children else list(glyph)
  double <- length(lines) == 2L
  centres <- y + if (double) c(2, -2) else 0
  grobs <- lapply(seq_along(lines), function(i) {
    line <- lines[[i]]
    gp <- ggchord_legend_pen(part, line$gp, if (double) 1.8 else 2.6, .9)
    gp$lineend <- part$params$lineend %||% "round"
    arrow <- line$arrow
    if (!is.null(arrow)) arrow$length <- arrow$length * scale
    segmentsGrob(ggchord_legend_x(8), ggchord_legend_y(centres[i]),
      ggchord_legend_x(57), ggchord_legend_y(centres[i]), gp = gp, arrow = arrow)
  })
  list(grob = do.call(grobTree, grobs), upper = min(centres), lower = max(centres))
}

# Annotation boxes are left, top, width and height. Fit the final outlines,
# keeping the head and body proportions of non-default feature shapes.
ggchord_legend_annotation <- function(part, box, opacity, scale) {
  glyph <- part$glyph
  left <- box[1L]
  top <- box[2L]
  width <- box[3L]
  height <- box[4L]
  strand <- as.character(part$style$strand %||% "+")
  reverse <- identical(strand, "-")
  edge <- c(left + if (reverse) 0 else width, top + height / 2)
  site_x <- left + width * if (reverse) .22 else .78
  if (identical(part$shape, "block")) {
    gp <- ggchord_legend_pen(part, glyph$gp, .75, .25)
    if (!part$alpha_fixed) gp$fill <- ggchord_legend_translucent(gp$fill, opacity)
    rect <- grid::roundrectGrob(ggchord_legend_x(left + width / 2), ggchord_legend_y(top + height / 2),
      ggchord_legend_x(width), unit(height / 52, "npc"), r = unit(.12, "mm") * scale, gp = gp)
    return(list(grob = rect, edge = edge, site = c(site_x, top)))
  }
  outlines <- if (inherits(glyph, "gTree")) glyph$children else list(glyph)
  polygons <- vapply(outlines, inherits, logical(1), "polygon")
  if (identical(part$shape, "primer_arc") && strand %in% c("+", "-") && all(polygons)) {
    gp <- ggchord_legend_pen(part, outlines[[1L]]$gp, 1.8, .25)
    gp$col <- if (is.null(gp$col) || is.na(gp$col)) gp$fill else gp$col
    gp$lineend <- "round"
    x <- function(t) left + width * if (reverse) (1 - t) else t
    profile <- rbind(c(x(0), top + .80 * height), c(x(.24), top + .20 * height),
      c(x(.52), top + .10 * height), c(x(.89), top + .69 * height))
    curve <- ggchord_legend_curve(profile)
    stroke <- grid::polylineGrob(ggchord_legend_x(curve[, 1L]), ggchord_legend_y(curve[, 2L]), gp = gp)
    gp <- ggchord_legend_pen(part, gp, .3, .25)
    tip <- rbind(c(x(.67), top), c(x(1), top + .78 * height), c(x(.65), top + height))
    return(list(grob = grobTree(stroke, ggchord_legend_polygon(tip, gp)), edge = tip[2L, ],
      site = c(x(.77), top + .24 * height)))
  }
  if (all(polygons)) {
    xs <- range(unlist(lapply(outlines, function(g) as.numeric(g$x))))
    ys <- range(unlist(lapply(outlines, function(g) as.numeric(g$y))))
    result <- lapply(outlines, function(g) {
      x <- left + (as.numeric(g$x) - xs[1L]) / diff(xs) * width
      y <- top + (ys[2L] - as.numeric(g$y)) / diff(ys) * height
      gp <- ggchord_legend_pen(part, g$gp, .65, .25)
      if (!part$alpha_fixed) gp$fill <- ggchord_legend_translucent(gp$fill, opacity)
      ggchord_legend_polygon(cbind(x, y), gp)
    })
    x <- unlist(lapply(result, function(g) as.numeric(g$x) * 65))
    y <- unlist(lapply(result, function(g) (1 - as.numeric(g$y)) * 52))
    end <- abs(x - edge[1L]) < 1e-7
    edge[2L] <- mean(range(y[end]))
    site <- c(x[which.min(y)], min(y))
    if (identical(part$shape, "arrow") && length(result) == 1L) {
      x1 <- c(tail(x, -1L), x[1L])
      y1 <- c(tail(y, -1L), y[1L])
      crosses <- x != x1 & pmin(x, x1) <= site_x & pmax(x, x1) >= site_x
      if (any(crosses)) {
        f <- (site_x - x[crosses]) / (x1[crosses] - x[crosses])
        site <- c(site_x, min(y[crosses] + f * (y1[crosses] - y[crosses])))
      }
    }
    return(list(grob = do.call(grobTree, result), edge = edge, site = site))
  }
  # Lollipops retain a round head even inside a narrow annotation box.
  v <- c((left + width / 2) / 65, 1 - (top + height / 2) / 52, width / 65, height / 52)
  list(grob = ggchord_legend_slot(glyph, v, scale),
    edge = c(left + width / 2, top + .81 * height),
    site = c(left + width / 2, top + .19 * height))
}

ggchord_legend_site_artist <- function(part, foot, marker_y) {
  glyph <- part$glyph
  gp <- glyph$gp
  if (inherits(glyph, "points"))
    gp$lwd <- (part$style$linewidth %||% .4) * ggplot2::.pt
  gp <- ggchord_legend_pen(part, gp, 1.25, .4)
  gp$lineend <- "round"
  if (inherits(glyph, "polyline")) {
    x <- foot[1L] + c(-6, -6, 0, 0)
    y <- c(marker_y, mean(c(marker_y, foot[2L])), mean(c(marker_y, foot[2L])), foot[2L])
    return(grid::polylineGrob(ggchord_legend_x(x), ggchord_legend_y(y), gp = gp))
  }
  stem <- segmentsGrob(ggchord_legend_x(foot[1L]), ggchord_legend_y(foot[2L]),
    ggchord_legend_x(foot[1L]), ggchord_legend_y(marker_y), gp = gp)
  if (!part$site_marker && !inherits(glyph, "points")) return(stem)
  if (inherits(glyph, "points")) {
    point <- ggchord_legend_slot(glyph, c(foot[1L] / 65, 1 - marker_y / 52, 1, 1), .36)
  } else {
    radius <- 3.6 * (part$style$size %||% 2) / 2
    gp$fill <- gp$col
    gp <- ggchord_legend_pen(part, gp, .55, .4)
    xy <- cbind(foot[1L] + c(0, radius, 0, -radius), marker_y + c(-radius, 0, radius, 0))
    point <- ggchord_legend_polygon(xy, gp)
  }
  grobTree(stem, point)
}

ggchord_legend_ribbon_artist <- function(part, left, right, opacity) {
  straight <- identical(part$route, "straight")
  a <- ggchord_legend_curve(left, straight)
  b <- ggchord_legend_curve(right, straight)
  centre <- (a + b) / 2
  t <- seq(0, 1, length.out = nrow(a))
  tip <- part$params$tip %||% "none"
  taper <- rep(1, length(t))
  if (tip %in% c("from", "both")) taper <- pmin(taper, t / .20)
  if (tip %in% c("to", "both")) taper <- pmin(taper, (1 - t) / .20)
  a <- centre + (a - centre) * taper
  b <- centre + (b - centre) * taper
  gp <- ggchord_legend_pen(part, part$glyph$gp, .8, .18)
  if (!part$alpha_fixed) {
    gp$fill <- ggchord_legend_translucent(gp$fill, opacity / .6)
    gp$col <- ggchord_legend_translucent(gp$col, .50 / .6)
  }
  if (inherits(gp$fill, "GridLinearGradient")) {
    reverse <- tail(centre[, 1L], 1L) < centre[1L, 1L]
    gp$fill$x1 <- unit(as.numeric(reverse), "npc")
    gp$fill$x2 <- unit(as.numeric(!reverse), "npc")
    gp$fill$y1 <- unit(1, "npc")
    gp$fill$y2 <- unit(0, "npc")
  }
  list(grob = ggchord_legend_polygon(rbind(a, b[rev(seq_len(nrow(b))), ]), gp), centre = centre)
}

ggchord_legend_link_artist <- function(part, points, scale) {
  gp <- ggchord_legend_pen(part, part$glyph$gp, 1.65, .5)
  gp$lineend <- part$params$lineend %||% "round"
  gp$linejoin <- "round"
  arrow <- part$glyph$arrow
  if (!is.null(arrow)) arrow$length <- arrow$length * scale
  grid::polylineGrob(ggchord_legend_x(points[, 1L]), ggchord_legend_y(points[, 2L]), gp = gp, arrow = arrow)
}

ggchord_legend_annotation_boxes <- function(roles, backbone, relations, nested) {
  gene <- "gene" %in% roles
  feature <- "feature" %in% roles
  if (backbone) {
    boxes <- list(gene = c(22, if (relations) 15.5 else 17.5, 24, 7),
      feature = c(22, if (relations) 10 else 13, 25, 8),
      primer = c(19, if (relations) 9 else 14, 27, 9))
    if (gene || feature) boxes$primer <- c(20, 22, 27, 7)
    if (gene && feature) {
      boxes$feature <- c(17, 8, 28, 14)
      boxes$gene <- c(24, 14.5, 19, 7)
    }
  } else {
    row <- if (relations) 20 else 25
    boxes <- list(gene = c(16, row - 4.5, 30, 9),
      feature = c(16, row - 5, 32, 10),
      primer = c(19, row - 7, 27, 9))
    if (gene || feature) boxes$primer <- c(20, row + 7, 27, 9)
    if (gene && feature) {
      boxes$feature <- c(13, row - 8, 34, 14)
      boxes$gene <- c(23, row - 4.5, 19, 7)
    }
  }
  if (gene && feature && !nested) {
    boxes$feature <- c(10, boxes$feature[2L], 20, 10)
    boxes$gene <- c(33, boxes$gene[2L], 21, 7)
  }
  boxes[roles]
}

ggchord_legend_can_nest <- function(gene, feature) {
  if (is.null(gene) || !gene$shape %in% ggchord_block_arrow_shapes ||
      !identical(feature$shape, "block")) return(FALSE)
  if (!feature$alpha_fixed) return(TRUE)
  fills <- c(gene$glyph$gp$fill, feature$glyph$gp$fill)
  if (length(fills) != 2L || !is.character(fills)) return(FALSE)
  rgba <- grDevices::col2rgb(fills, alpha = TRUE)
  rgba[4L, 2L] < 220L || any(rgba[1:3, 1L] != rgba[1:3, 2L])
}

ggchord_legend_relation_art <- function(parts, seq, anchor, from_annotation, from_site, scale) {
  ribbon <- parts$ribbon
  link <- parts$link
  width <- if (is.null(link)) .23 else .26
  left <- right <- control <- NULL
  if (!is.null(seq)) {
    sy <- seq$lower
    left <- rbind(c(21, sy), c(23, sy + 9), c(37, 43), c(44, 49))
    right <- rbind(c(34, sy), c(35, sy + 9), c(49, 34), c(55, 49))
    x <- if (is.null(parts$site) || from_annotation) 35 else 32
    control <- rbind(c(x, sy), c(x + 7, sy + 10), c(45, 41), c(54, 47))
  } else if (from_annotation) {
    direction <- if (anchor[1L] < 24) -1 else 1
    far <- function(x) if (direction == -1) 65 - x else x
    left <- rbind(anchor - c(direction * 8, 0), anchor + c(-direction * 8, 9),
      c(far(41), 40), c(far(44), 48))
    right <- rbind(anchor, anchor + c(direction * 2, 7), c(far(49), 37), c(far(55), 48))
    control <- rbind(anchor, anchor + c(7 * direction, 10),
      c(anchor[1L] - direction, 41), c(far(54), 47))
  } else if (from_site) {
    left <- rbind(anchor - c(3.5, 0), c(28, 34), c(41, 43), c(44, 48))
    right <- rbind(anchor + c(3.5, 0), c(36, 30), c(49, 39), c(55, 48))
    control <- rbind(anchor, c(37, 35), c(45, 41), c(54, 47))
  } else {
    left <- rbind(c(15, 10), c(19, 22), c(31, 37), c(43, 48))
    right <- rbind(c(31, 10), c(34, 19), c(44, 27), c(55, 48))
    control <- rbind(c(17, 12), c(24, 22), c(40, 38), c(54, 47))
  }
  band <- if (!is.null(ribbon)) ggchord_legend_ribbon_artist(ribbon, left, right, width)
  line <- if (!is.null(link)) {
    path <- band$centre %||% ggchord_legend_curve(control, identical(link$route, "straight"))
    ggchord_legend_link_artist(link, path, scale)
  }
  list(ribbon = band$grob, link = line)
}

ggchord_legend_artwork_view <- function(room, role) {
  # Keep the full template for a shared composition. In a separate slot,
  # fit around the single symbol so its unused template margins take no space.
  frame <- switch(role,
    seq = c(32.5, 35, 56, 24),
    gene = c(31, 25, 36, 18),
    feature = c(32, 25, 38, 18),
    primer = c(32.5, 23, 34, 18),
    site = c(30, 18, 22, 22),
    link = c(35.5, 29.5, 44, 43),
    ribbon = c(35, 29, 46, 44),
    c(32.5, 26, 65, 52))
  ratio <- min(room / (frame[3:4] / 10))
  if (nzchar(role)) ratio <- min(ratio, 1)
  extent <- ggchord_legend_artwork_mm * ratio / pmax(room, .Machine$double.eps)
  offset <- c(.5 - frame[1L] / 65, frame[2L] / 52 - .5) * extent
  # Thin strokes can vanish when the key is displayed at screen resolution.
  # Keep a legible stroke in ordinary slots and keep scaling in tight slots.
  ink <- if (nzchar(role)) max(ratio, min(room) / 3.6) else ratio
  list(box = c(.5 + offset, extent), scale = min(ratio, 1.25), ink = min(ink, 1.25))
}

ggchord_contact_key <- function(parts, room, separated = FALSE) {
  names(parts) <- vapply(parts, `[[`, character(1), "role")
  roles <- intersect(c("feature", "gene", "primer"), names(parts))
  has_annotations <- length(roles) > 0L
  has_seq <- !is.null(parts$seq)
  has_site <- !is.null(parts$site)
  has_relations <- any(c("link", "ribbon") %in% names(parts))
  role <- if (separated) names(parts)[1L] else ""
  view <- ggchord_legend_artwork_view(room, role)
  seq <- if (has_seq) ggchord_legend_seq_artist(parts$seq,
    if (has_annotations) 31 else if (has_relations) 18 else 35, view$scale)

  nested <- ggchord_legend_can_nest(parts$gene, parts$feature)
  boxes <- ggchord_legend_annotation_boxes(roles, has_seq, has_relations, nested)
  annotations <- lapply(roles, function(role) {
    opacity <- if (role == "gene") .94 else if (nested) {
      if (has_seq) .25 else .23
    } else if (has_seq) .78 else .82
    ggchord_legend_annotation(parts[[role]], boxes[[role]], opacity, view$scale)
  })
  names(annotations) <- roles
  anchor_annotation <- annotations$gene %||% annotations$feature %||% annotations$primer
  foot <- if (has_seq) c(if (has_annotations) 52 else 32, seq$upper) else
    anchor_annotation$site %||% c(30, 25)
  site <- if (has_site) ggchord_legend_site_artist(parts$site, foot,
    max(5.5, foot[2L] - if (has_seq) 14 else 12))
  connectors <- if (has_relations) ggchord_legend_relation_art(parts, seq,
    anchor_annotation$edge %||% foot, has_annotations, has_site, view$scale)

  grobs <- c(connectors, list(seq$grob), lapply(annotations, `[[`, "grob"), list(site))
  glyph <- do.call(grobTree, Filter(Negate(is.null), unname(grobs)))
  ggchord_legend_slot(glyph, view$box, view$ink)
}
