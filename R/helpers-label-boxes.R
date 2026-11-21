#' Wrap long gene annotation texts at a given character width
#' @noRd
ggchord_label_wrap_text <- function(text, width = NULL) {
  if (is.null(width) || width <= 0) return(text)
  vapply(text, function(t) {
    if (is.na(t) || !nzchar(t)) return(NA_character_)
    paste(strwrap(t, width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# Count collisions for each visible label against the other labels and any
# fixed annotation boxes. Keeping the row-wise counts in one helper ensures
# adaptive fitting and max_overlaps use exactly the same collision semantics.
ggchord_label_conflict_counts <- function(gl, units_per_inch = 0.35,
                                           box_padding = 0,
                                           repel_boxes = NULL) {
  counts <- integer(nrow(gl))
  active <- !is.na(gl$text) & nzchar(gl$text)
  rows <- which(active)
  if (length(rows) == 0L) return(counts)

  boxes <- ggchord_text_boxes(
    gl[rows, , drop = FALSE],
    units_per_inch = units_per_inch,
    box_padding = box_padding
  )
  if (length(rows) > 1L) {
    for (i in seq_len(nrow(boxes))) {
      other <- setdiff(seq_len(nrow(boxes)), i)
      counts[rows[i]] <- counts[rows[i]] + sum(
        ggchord_oriented_box_overlaps(
          boxes[i, , drop = FALSE], boxes[other, , drop = FALSE]
        )
      )
    }
  }
  if (!is.null(repel_boxes) && nrow(repel_boxes) > 0L) {
    for (i in seq_len(nrow(boxes))) {
      counts[rows[i]] <- counts[rows[i]] + sum(
        ggchord_oriented_box_overlaps(
          boxes[i, , drop = FALSE], repel_boxes
        )
      )
    }
  }
  counts
}

# Adapt labels whose initial boxes are crowded or physically too wide for a
# compact perimeter rail. This is deliberately a text-fitting step rather than
# another force/layout parameter: the selected layout remains deterministic and
# is rerun after the text dimensions change.
ggchord_fit_label_text <- function(gl,
                                   fit = c("wrap", "none", "ellipsis", "auto", "shrink"),
                                   max_lines = 2L,
                                   units_per_inch = 0.35,
                                   box_padding = 0.18,
                                   repel_boxes = NULL) {
  fit <- match.arg(fit)
  if (fit %in% c("none", "shrink") || nrow(gl) == 0L) return(gl)
  crowded <- ggchord_label_conflict_counts(
    gl,
    units_per_inch = units_per_inch,
    box_padding = box_padding,
    repel_boxes = repel_boxes
  ) > 0L
  measured <- ggchord_text_boxes(
    gl,
    units_per_inch = units_per_inch,
    box_padding = 0
  )
  # Roughly one inch is a useful upper bound for a default horizontal callout.
  # Measuring the rendered box (rather than counting characters) also works for
  # wide glyphs and non-Latin labels.
  wide <- measured$bw > 1.05 * units_per_inch
  rows <- which((crowded | wide) & !is.na(gl$text) & nzchar(gl$text))
  if (length(rows) == 0L) return(gl)

  fit_one <- function(value) {
    plain <- gsub("[\r\n]+", " ", value)
    plain <- trimws(gsub("[[:space:]]+", " ", plain))
    n_chars <- nchar(plain, type = "width")
    if (!is.finite(n_chars) || n_chars < 2L) return(plain)
    line_width <- max(8L, ceiling(n_chars / max_lines))

    if (identical(fit, "ellipsis")) {
      if (n_chars <= line_width) return(plain)
      return(paste0(substr(plain, 1L, max(1L, line_width - 1L)), "\u2026"))
    }

    lines <- strwrap(plain, width = line_width)
    if (length(lines) <= max_lines) return(paste(lines, collapse = "\n"))
    if (identical(fit, "wrap")) {
      lines <- c(
        lines[seq_len(max_lines - 1L)],
        paste(lines[max_lines:length(lines)], collapse = " ")
      )
      return(paste(lines, collapse = "\n"))
    }

    # auto: preserve as much text as the requested line count permits, then
    # mark only the final overflow with an ellipsis.
    kept <- lines[seq_len(max_lines)]
    kept[max_lines] <- paste0(
      substr(kept[max_lines], 1L, max(1L, line_width - 1L)), "\u2026"
    )
    paste(kept, collapse = "\n")
  }

  gl$text[rows] <- vapply(gl$text[rows], fit_one, character(1), USE.NAMES = FALSE)
  gl
}

# Fit complete text into a declared feature body. Only explicit shrink uses
# this path; it never shortens the text or lets a minimum become a target size.
ggchord_shrink_label_text <- function(gl, available_width, available_height,
                                      min_size = 1.5, units_per_inch = .35) {
  if (!nrow(gl)) return(gl)
  if (!is.numeric(min_size) || length(min_size) != 1L ||
      !is.finite(min_size) || min_size <= 0) ggchord_stop(
    "label_min_size must be one finite positive text size in millimetres")
  if (any(!is.finite(gl$size)) || any(gl$size < min_size)) ggchord_stop(
    "label_min_size cannot exceed the requested text size")
  if (!"label_text_original" %in% names(gl)) gl$label_text_original <- as.character(gl$text)
  gl$label_size_requested <- gl$size
  gl$label_size_min <- min_size
  gl$label_fit_reason <- "inside_unchanged"
  gl$label_inside_fit <- FALSE
  available_width <- rep_len(available_width, nrow(gl))
  available_height <- rep_len(available_height, nrow(gl))
  fits <- function(row, index) {
    box <- ggchord_text_boxes(row, units_per_inch = units_per_inch,
      box_padding = 0)
    all(is.finite(c(available_width[index], available_height[index]))) &&
      available_width[index] > 0 && available_height[index] > 0 &&
      box$w <= available_width[index] && box$h <= available_height[index]
  }
  for (i in seq_len(nrow(gl))) {
    row <- gl[i, , drop = FALSE]
    if (is.na(row$text) || !nzchar(row$text)) next
    if (fits(row, i)) {
      gl$label_inside_fit[i] <- TRUE
      next
    }
    row$size <- min_size
    if (!fits(row, i)) {
      # External fallback keeps the requested font for readability. Failure to
      # fit at the lower bound is reported rather than shrinking past it.
      gl$label_fit_reason[i] <- "minimum_size_exceeded"
      next
    }
    lower <- min_size; upper <- gl$size[i]
    # Physical glyph metrics can be quantized by the backend. A bounded binary
    # search is sufficient and is remeasured on each chosen device/context.
    for (iteration in seq_len(10L)) {
      middle <- (lower + upper) / 2
      row$size <- middle
      if (fits(row, i)) lower <- middle else upper <- middle
    }
    gl$size[i] <- lower
    gl$label_inside_fit[i] <- TRUE
    gl$label_fit_reason[i] <- "inside_shrunk"
  }
  gl$label_size_final <- gl$size
  gl
}

ggchord_label_fit_metadata <- function(gl) {
  if (!nrow(gl)) return(gl)
  if (!"label_text_original" %in% names(gl)) gl$label_text_original <- gl$text
  if (!"label_size_requested" %in% names(gl)) gl$label_size_requested <- gl$size
  gl$label_size_final <- gl$size
  if (!"label_fit_reason" %in% names(gl)) gl$label_fit_reason <- "unchanged"
  if ("feature_label_mode" %in% names(gl)) {
    outside <- !gl$feature_label_mode %in% "inside"
    gl$label_fit_reason[outside & gl$label_fit_reason %in% c("inside_shrunk", "inside_unchanged")] <-
      "collision_repositioned"
  }
  hidden <- is.na(gl$text) | !nzchar(gl$text)
  gl$label_fit_reason[hidden] <- ifelse(gl$label_fit_reason[hidden] == "minimum_size_exceeded",
    "minimum_size_exceeded_hidden", "layout_hidden")
  gl
}

#' De-overlap gene labels
#'
#' Detects overlapping gene label boxes (estimated from the text size) and
#' pushes the labels apart until they no longer collide. Optionally hides
#' labels that still overlap more than `max_overlaps` other labels
#' (ggrepel-style decluttering).
#' @noRd
ggchord_label_deoverlap <- function(gl, units_per_inch = 0.35, seed = 123,
                                    max_overlaps = Inf) {
  if (nrow(gl) < 2) return(gl)
  measured <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = .025
  )
  w <- measured$bw
  h <- measured$bh
  x <- measured$cx
  y <- measured$cy
  n <- nrow(gl)

  # Resolve axis-aligned box overlaps iteratively, pushing labels apart along
  # the axis of least penetration.
  for (iter in seq_len(100)) {
    moved <- FALSE
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        dx <- x[j] - x[i]
        dy <- y[j] - y[i]
        ox <- (w[i] + w[j]) / 2 - abs(dx)
        oy <- (h[i] + h[j]) / 2 - abs(dy)
        if (ox > 0 && oy > 0) {
          if (ox < oy) {
            sgn <- if (dx >= 0) 1 else -1
            x[i] <- x[i] - sgn * ox / 2
            x[j] <- x[j] + sgn * ox / 2
          } else {
            sgn <- if (dy >= 0) 1 else -1
            y[i] <- y[i] - sgn * oy / 2
            y[j] <- y[j] + sgn * oy / 2
          }
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }

  # Optional decluttering: hide labels that still overlap too many others
  if (is.finite(max_overlaps)) {
    n_over <- numeric(n)
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        if (abs(x[i] - x[j]) < (w[i] + w[j]) / 2 &&
            abs(y[i] - y[j]) < (h[i] + h[j]) / 2) {
          n_over[i] <- n_over[i] + 1
          n_over[j] <- n_over[j] + 1
        }
      }
    }
    hide <- n_over > max_overlaps
    if (any(hide)) gl$text[hide] <- NA
  }

  gl$text_x <- gl$text_x + x - measured$cx
  gl$text_y <- gl$text_y + y - measured$cy
  gl
}

#' Hide fixed labels that collide, preserving deterministic input priority
#'
#' Unlike the automatic repel layouts, fixed labels must not silently acquire
#' leader lines or move far away from their feature. This greedy filter keeps
#' the first label in input order and omits only later labels whose oriented
#' text box intersects one that has already been retained.
#' @noRd
ggchord_label_prune_overlaps <- function(gl, units_per_inch = 0.35,
                                         box_padding = 0.015,
                                         repel_boxes = NULL) {
  if (nrow(gl) < 2) return(gl)
  active <- which(!is.na(gl$text) & nzchar(gl$text))
  if (length(active) < 2) return(gl)

  boxes <- ggchord_text_boxes(
    gl[active, , drop = FALSE], units_per_inch = units_per_inch,
    box_padding = box_padding
  )
  kept <- integer(0)
  for (i in seq_along(active)) {
    conflict <- (length(kept) > 0 && any(ggchord_oriented_box_overlaps(
      boxes[i, , drop = FALSE], boxes[kept, , drop = FALSE]
    ))) || (!is.null(repel_boxes) && nrow(repel_boxes) > 0 &&
      any(ggchord_oriented_box_overlaps(
        boxes[i, , drop = FALSE], repel_boxes
      )))
    if (conflict) {
      gl$text[active[i]] <- NA_character_
    } else {
      kept <- c(kept, i)
    }
  }
  gl
}
#' Estimate the axis-aligned text boxes for a set of text labels.
#'
#' `text_x`/`text_y` are the points selected by `hjust`/`vjust`, not the
#' rendered text centre.  This helper returns both the anchor (`x`/`y`) and the
#' centre (`cx`/`cy`) together with the full axis-aligned width/height of the
#' text box (`bw`/`bh`).  The same projection is used by the repulsion solver,
#' the obstacle boxes and the adaptive coordinate limits so all three agree.
#' @noRd
ggchord_text_boxes <- function(df,
                               x_col = "text_x", y_col = "text_y",
                               text_col = "text", angle_col = "text_angle",
                               size_col = "size", hjust_col = "hjust",
                               vjust_col = "vjust",
                               family_col = "family",
                               fontface_col = "fontface",
                               lineheight_col = "lineheight",
                               units_per_inch = 0.35, box_padding = 0) {
  n <- nrow(df)
  if (!n) return(ggchord_project_text_boxes(numeric(), numeric(), numeric(),
    numeric(), numeric(), numeric(), numeric(), units_per_inch, box_padding))

  x <- df[[x_col]]
  y <- df[[y_col]]
  texts <- df[[text_col]]
  sizes <- if (is.null(df[[size_col]])) rep(2.5, n) else df[[size_col]]
  angles <- (if (is.null(df[[angle_col]])) rep(0, n) else df[[angle_col]]) *
    pi / 180
  hjust <- if (is.null(df[[hjust_col]])) rep(0.5, n) else df[[hjust_col]]
  vjust <- if (is.null(df[[vjust_col]])) rep(0.5, n) else df[[vjust_col]]
  families <- if (is.null(df[[family_col]])) rep("", n) else
    as.character(df[[family_col]])
  fontfaces <- if (is.null(df[[fontface_col]])) rep(1, n) else
    df[[fontface_col]]
  fontfaces[is.na(fontfaces)] <- if (is.numeric(fontfaces)) 1L else "plain"
  lineheights <- if (is.null(df[[lineheight_col]])) rep(1.2, n) else
    as.numeric(df[[lineheight_col]])
  families[is.na(families)] <- ""
  lineheights[!is.finite(lineheights) | lineheights <= 0] <- 1.2

  w <- numeric(n)
  h <- numeric(n)
  valid <- !is.na(texts) & nzchar(texts)
  composite <- if ('.ggchord_composite_label' %in% names(df) &&
      all(c('enzyme_label', 'coordinate_label', 'enzyme_fontface') %in% names(df)))
    which(df$.ggchord_composite_label %in% TRUE & valid) else integer()
  uniform_sites <- if ('.ggchord_composite_label' %in% names(df))
    which(df$.ggchord_composite_label %in% FALSE & valid) else integer()
  if (any(valid)) {
    legacy_device <- names(dev.cur()) %in% c("pdf", "postscript")
    unicode <- any(utf8ToInt(paste(texts[valid], collapse = "")) > 127L)
    close_device <- if (dev.cur() > 1L && !(legacy_device && unicode))
      function() invisible(NULL) else ggchord_measurement_device(text = texts[valid])
    on.exit(close_device())
    # Site labels are measured below using the same fragments as their
    # renderer. A generic flat label measurement would be discarded.
    regular <- valid
    regular[c(composite, uniform_sites)] <- FALSE
    for (i in which(regular)) {
      label <- if (isTRUE(df$.axis_plotmath[i])) parse(text = texts[i]) else texts[i]
      grob <- textGrob(
        label, gp = gpar(
          fontsize = sizes[i] * (72.27 / 25.4),
          fontfamily = families[i], fontface = fontfaces[i],
          lineheight = lineheights[i]
        )
      )
      w[i] <- convertWidth(
        grobWidth(grob), "inches", valueOnly = TRUE
      ) * units_per_inch
      h[i] <- convertHeight(
        grobHeight(grob), "inches", valueOnly = TRUE
      ) * units_per_inch
    }
    close_device()
  }

  if (length(composite)) {
    exact <- ggchord_site_label_metrics(df$enzyme_label[composite], df$coordinate_label[composite],
      sizes[composite], families[composite], df$enzyme_fontface[composite],
      lineheights[composite], units_per_inch)
    w[composite] <- exact$width; h[composite] <- exact$height
    hjust[composite] <- .5
  }
  if (length(uniform_sites)) {
    exact <- ggchord_restriction_text_metrics(texts[uniform_sites], sizes[uniform_sites],
      units_per_inch, parse = FALSE, fontface = fontfaces[uniform_sites],
      family = families[uniform_sites], lineheight = lineheights[uniform_sites])
    w[uniform_sites] <- exact$width; h[uniform_sites] <- exact$height
  }
  ggchord_project_text_boxes(x, y, w, h, angles, hjust, vjust, units_per_inch, box_padding)
}

ggchord_project_text_boxes <- function(x, y, w, h, angles, hjust, vjust,
                                       units_per_inch, box_padding) {
  cos_a <- cos(angles)
  sin_a <- sin(angles)
  cx_off <- (0.5 - hjust) * w * cos_a - (0.5 - vjust) * h * sin_a
  cy_off <- (0.5 - hjust) * w * sin_a + (0.5 - vjust) * h * cos_a
  bw <- abs(cos_a) * w + abs(sin_a) * h + 2 * box_padding * units_per_inch
  bh <- abs(sin_a) * w + abs(cos_a) * h + 2 * box_padding * units_per_inch

  columns <- list(
    x = x, y = y,
    cx = x + cx_off, cy = y + cy_off,
    w = w, h = h,
    ow = w + 2 * box_padding * units_per_inch,
    oh = h + 2 * box_padding * units_per_inch,
    angle = angles,
    bw = bw, bh = bh,
    xmin = x + cx_off - bw / 2, xmax = x + cx_off + bw / 2,
    ymin = y + cy_off - bh / 2, ymax = y + cy_off + bh / 2
  )
  # Glyph vectors usually already have equal lengths and automatic row names.
  # Keep normal data.frame recycling and named row semantics for other inputs.
  if (all(lengths(columns) == length(x)) && all(vapply(columns,
      function(v) is.null(names(v)) && !is.object(v), logical(1))))
    return(list2DF(columns))
  do.call(data.frame, columns)
}

# Lay one or more labels out glyph by glyph on the circle selected by their
# text anchor. This is shared by collision detection, leader clipping and the
# renderer so all three use the same curved baseline and local tangents.
ggchord_arc_text_layout <- function(text, units_per_inch = .35,
                                     box_padding = 0, respect_justification = FALSE, line_metrics = NULL, boxes_only = FALSE) {
  if (!nrow(text)) return(list(glyphs = data.frame(), boxes =
    ggchord_text_boxes(data.frame()), owners = integer()))
  glyphs <- list()
  boxes <- list()
  owners <- integer()
  for (i in seq_len(nrow(text))) {
    if (!is.na(text$text[i]) && grepl("\n", text$text[i], fixed=TRUE)) {
      lines <- strsplit(as.character(text$text[i]), "\n", fixed=TRUE)[[1L]]
      centre <- atan2(text$text_y[i],text$text_x[i])
      radius <- sqrt(text$text_x[i]^2+text$text_y[i]^2)
      tangent <- centre*180/pi+90
      difference <- ((text$text_angle[i]-tangent+180)%%360)-180
      direction <- if(abs(difference)>90) -1 else 1
      sample <- text[i,,drop=FALSE];sample$text <- "Mg"
      sample$text_angle <- 0
      pitch <- ggchord_text_boxes(sample,units_per_inch=units_per_inch)$h *
        (text$lineheight[i] %||% 1.2)
      for (line in seq_along(lines)) {
        row <- text[i,,drop=FALSE];row$text <- lines[line]
        line_radius <- radius-direction*((length(lines)+1)/2-line)*pitch
        row$text_x <- line_radius*cos(centre);row$text_y <- line_radius*sin(centre)
        part <- ggchord_arc_text_layout(row, units_per_inch, box_padding, respect_justification, line_metrics, boxes_only)
        if (nrow(part$boxes)) {
          if (!boxes_only) {
            part$glyphs$.arc_parent_group <- text$group[i] %||% i
            part$glyphs$.arc_line <- line
            glyphs[[length(glyphs)+1L]] <- part$glyphs
          }
          boxes[[length(boxes)+1L]] <- part$boxes
          owners <- c(owners,rep(i,nrow(part$boxes)))
        }
      }
      next
    }
    chars <- strsplit(enc2utf8(as.character(text$text[i])), "", fixed = TRUE)[[1L]]
    if (!length(chars) || is.na(text$text[i]) || !nzchar(text$text[i])) next
    # All glyphs on this line inherit the same type settings and viewport.
    # Keep exact per-glyph sizes for the final rotation and obstacle boxes.
    unique_chars <- unique(chars)
    key <- as.character(text$text[i])
    metrics <- if (!is.null(line_metrics)) line_metrics$lines[[key]] else NULL
    if (is.null(metrics)) {
      measured <- text[rep(i, length(unique_chars)), , drop = FALSE]
      measured$text <- unique_chars
      measured$text_x <- 0
      measured$text_y <- 0
      measured$text_angle <- 0
      measured$hjust <- .5
      measured$vjust <- .5
      metrics <- ggchord_text_boxes(
        measured, units_per_inch = units_per_inch, box_padding = 0
      )
      if (!is.null(line_metrics)) line_metrics$lines[[key]] <- metrics[c("w", "h")]
    }
    index <- match(chars, unique_chars)
    width <- metrics$w[index]
    height <- metrics$h[index]
    advance <- width
    fallback <- max(height, na.rm = TRUE) * .34
    advance[!is.finite(advance) | advance <= 0] <- fallback
    tracking <- max(.0008, median(height, na.rm = TRUE) * .018)
    total <- sum(advance) + tracking * max(0, length(chars) - 1L)
    along <- cumsum(c(0, head(advance + tracking, -1L))) + advance / 2 -
      total / 2
    if (respect_justification) along <- along + (.5 - text$hjust[i]) * total
    radius <- sqrt(text$text_x[i]^2 + text$text_y[i]^2)
    centre_angle <- atan2(text$text_y[i], text$text_x[i])
    increasing_tangent <- (centre_angle * 180 / pi + 90) %% 360
    centre_orientation <- text$text_angle[i] %% 360
    orientation_delta <- abs((centre_orientation - increasing_tangent + 180) %%
      360 - 180)
    direction <- if (orientation_delta > 90) -1 else 1
    glyph_angle <- centre_angle + direction * along / max(radius, .1)
    rotation <- glyph_angle * 180 / pi + if (direction > 0) 90 else -90
    # Decide readability once for the complete string.  Normalising every
    # glyph independently made labels that cross a quadrant boundary reverse
    # direction part-way through a word.
    centre_rotation <- centre_angle * 180 / pi +
      if (direction > 0) 90 else -90
    centre_upright <- ggchord_normalize_text_orientation(
      centre_rotation, .5
    )
    string_flip <- ((centre_upright$angle - centre_rotation + 180) %% 360) - 180
    x <- radius * cos(glyph_angle)
    y <- radius * sin(glyph_angle)
    angle <- rotation + string_flip
    if (!boxes_only) {
      glyph <- text[rep(i, length(chars)), , drop = FALSE]
      glyph$text_x <- x; glyph$text_y <- y; glyph$text_angle <- angle
      glyph$hjust <- .5; glyph$vjust <- .5; glyph$text <- chars
      glyph$.arc_parent_group <- text$group[i] %||% i
      glyph$.arc_glyph_index <- seq_along(chars)
      glyphs[[length(glyphs) + 1L]] <- glyph
    }
    boxes[[length(boxes) + 1L]] <- ggchord_project_text_boxes(
      x, y, width, height, angle * pi / 180, .5, .5, units_per_inch, box_padding)
    owners <- c(owners, rep(i, length(chars)))
  }
  list(
    glyphs = if (length(glyphs)) ggchord_rbind_fill(glyphs) else data.frame(),
    boxes = if (length(boxes) == 1L) boxes[[1L]] else if (length(boxes)) do.call(rbind, boxes) else
      ggchord_text_boxes(data.frame()),
    owners = owners
  )
}

ggchord_boxes_overlap_any <- function(candidate, other) {
  if (!nrow(candidate) || is.null(other) || !nrow(other)) return(FALSE)
  # Discard boxes that cannot intersect any glyph in this candidate before
  # testing each glyph individually. Use the same tolerance in both checks.
  tol <- 1e-7
  if (nrow(candidate) > 1L) {
    bounds <- c(min(candidate$xmin), max(candidate$xmax),
      min(candidate$ymin), max(candidate$ymax))
    if (!anyNA(bounds)) {
      possible <- which(bounds[1L] < other$xmax - tol & bounds[2L] > other$xmin + tol &
        bounds[3L] < other$ymax - tol & bounds[4L] > other$ymin + tol)
      if (!length(possible)) return(FALSE)
      if (length(possible) < nrow(other)) other <- other[possible, , drop = FALSE]
    }
  }
  for (i in seq_len(nrow(candidate))) {
    if (any(ggchord_oriented_box_overlaps(candidate, other, tol = tol, i = i))) return(TRUE)
  }
  FALSE
}

ggchord_polygon_bounds <- function(polygons) {
  edges <- vapply(polygons, function(p) {
    if (!nrow(p)) return(c(Inf, -Inf, Inf, -Inf))
    c(range(p$x), range(p$y))
  }, numeric(4L))
  data.frame(xmin = edges[1L, ], xmax = edges[2L, ],
             ymin = edges[3L, ], ymax = edges[4L, ], row.names = NULL)
}

ggchord_boxes_hit_features <- function(boxes, polygons, source_row = NA_integer_,
                                        allow_own = FALSE, tol = .004, bounds = NULL) {
  if (!length(polygons)) return(FALSE)
  for (i in seq_len(nrow(boxes))) {
    candidates <- polygons
    if (!is.null(bounds)) {
      possible <- boxes$xmax[i] >= bounds$xmin - tol & boxes$xmin[i] <= bounds$xmax + tol &
        boxes$ymax[i] >= bounds$ymin - tol & boxes$ymin[i] <= bounds$ymax + tol
      if (!any(possible)) next
      candidates <- polygons[possible]
    }
    if (ggchord_label_hits_features(
      boxes[i, , drop = FALSE], candidates, source_row, allow_own, tol
    )) return(TRUE)
  }
  FALSE
}

ggchord_boxes_corners <- function(boxes) {
  if (!nrow(boxes)) return(matrix(numeric(), ncol = 2L))
  row <- rep(seq_len(nrow(boxes)), each = 4L)
  angle <- boxes$angle[row]
  along <- c(-1, -1, 1, 1) * (boxes$ow[row] / 2)
  across <- c(-1, 1, 1, -1) * (boxes$oh[row] / 2)
  cbind(boxes$cx[row] + along * cos(angle) - across * sin(angle),
        boxes$cy[row] + along * sin(angle) + across * cos(angle))
}

# Test one oriented text rectangle against zero or more oriented rectangles.
# A cheap axis-aligned prefilter is followed by the separating-axis theorem;
# this avoids treating a diagonal label's large empty corner triangles as
# occupied space.
ggchord_oriented_box_overlaps <- function(candidate, other, tol = 1e-7, i = 1L) {
  if (is.null(other) || nrow(other) == 0) return(logical(0))
  possible <- candidate$xmin[i] < other$xmax - tol &
    candidate$xmax[i] > other$xmin + tol &
    candidate$ymin[i] < other$ymax - tol &
    candidate$ymax[i] > other$ymin + tol
  out <- rep(FALSE, nrow(other))
  j <- which(possible)
  if (!length(j)) return(out)
  angle <- (candidate$angle[i] %||% 0)
  other_angle <- other$angle[j] %||% 0
  ca <- cos(angle); sa <- sin(angle)
  cb <- cos(other_angle); sb <- sin(other_angle)
  x <- other$cx[j] - candidate$cx[i]
  y <- other$cy[j] - candidate$cy[i]
  aw <- (candidate$ow[i] %||% candidate$bw[i]) / 2
  ah <- (candidate$oh[i] %||% candidate$bh[i]) / 2
  bw <- (other$ow[j] %||% other$bw[j]) / 2
  bh <- (other$oh[j] %||% other$bh[j]) / 2
  # Project both rectangles on all four local axes. Keep the same dot
  # products (including the squared unit lengths) and boundary tolerance.
  ac <- ca * ca + sa * sa
  bc <- cb * cb + sb * sb
  parallel <- abs(ca * cb + sa * sb)
  perpendicular <- abs(sa * cb - ca * sb)
  out[j] <- abs(x * ca + y * sa) < aw * ac + bw * parallel + bh * perpendicular - tol &
    abs(-x * sa + y * ca) < ah * ac + bw * perpendicular + bh * parallel - tol &
    abs(x * cb + y * sb) < aw * parallel + ah * perpendicular + bw * bc - tol &
    abs(-x * sb + y * cb) < aw * perpendicular + ah * parallel + bh * bc - tol
  out
}

ggchord_point_in_polygon <- function(x, y, polygon) {
  px <- polygon$x
  py <- polygon$y
  n <- length(px)
  if (n < 3L) return(FALSE)
  previous <- c(n, seq_len(n - 1L))
  crosses <- which((py > y) != (py[previous] > y))
  from <- previous[crosses]
  edge_x <- (px[from] - px[crosses]) * (y - py[crosses]) /
    (py[from] - py[crosses]) + px[crosses]
  sum(x < edge_x) %% 2L == 1L
}

ggchord_box_corners <- function(box) {
  angle <- box$angle[1] %||% 0
  axes <- rbind(c(cos(angle), sin(angle)), c(-sin(angle), cos(angle)))
  signs <- rbind(c(-1, -1), c(-1, 1), c(1, 1), c(1, -1))
  cbind(
    box$cx[1] + signs[, 1] * box$ow[1] / 2 * axes[1, 1] +
      signs[, 2] * box$oh[1] / 2 * axes[2, 1],
    box$cy[1] + signs[, 1] * box$ow[1] / 2 * axes[1, 2] +
      signs[, 2] * box$oh[1] / 2 * axes[2, 2]
  )
}

ggchord_segment_intersects <- function(a, b, c, d, tol = 1e-10) {
  cross <- function(u, v) u[1] * v[2] - u[2] * v[1]
  r <- b - a
  s <- d - c
  denominator <- cross(r, s)
  if (abs(denominator) <= tol) return(FALSE)
  offset <- c - a
  t <- cross(offset, s) / denominator
  u <- cross(offset, r) / denominator
  t >= -tol && t <= 1 + tol && u >= -tol && u <= 1 + tol
}

ggchord_box_overlaps_feature <- function(box, polygon, tol = .004) {
  if (nrow(polygon) < 3L) return(FALSE)
  if (box$xmax[1] < min(polygon$x) - tol ||
      box$xmin[1] > max(polygon$x) + tol ||
      box$ymax[1] < min(polygon$y) - tol ||
      box$ymin[1] > max(polygon$y) + tol) return(FALSE)
  corners <- ggchord_box_corners(box)
  if (any(vapply(seq_len(nrow(corners)), function(i) {
    ggchord_point_in_polygon(corners[i, 1], corners[i, 2], polygon)
  }, logical(1)))) return(TRUE)
  if (ggchord_point_in_polygon(box$cx[1], box$cy[1], polygon)) return(TRUE)
  closed_corners <- rbind(corners, corners[1L, , drop = FALSE])
  # Test one box edge against all polygon edges using the same intersection
  # parameters and tolerance as ggchord_segment_intersects().
  px <- polygon$x
  py <- polygon$y
  following <- c(seq.int(2L, length(px)), 1L)
  edge_x <- px[following] - px
  edge_y <- py[following] - py
  for (i in seq_len(nrow(closed_corners) - 1L)) {
    start <- closed_corners[i, ]
    step <- closed_corners[i + 1L, ] - start
    denominator <- step[1] * edge_y - step[2] * edge_x
    use <- which(abs(denominator) > tol)
    dx <- px[use] - start[1]
    dy <- py[use] - start[2]
    along_box <- (dx * edge_y[use] - dy * edge_x[use]) / denominator[use]
    along_polygon <- (dx * step[2] - dy * step[1]) / denominator[use]
    if (any(along_box >= -tol & along_box <= 1 + tol &
            along_polygon >= -tol & along_polygon <= 1 + tol)) return(TRUE)
  }
  # Dense feature paths make a vertex-in-box check an inexpensive and robust
  # boundary-intersection test for the curved polygons produced by ggchord.
  local_x <- (polygon$x - box$cx[1]) * cos(box$angle[1]) +
    (polygon$y - box$cy[1]) * sin(box$angle[1])
  local_y <- -(polygon$x - box$cx[1]) * sin(box$angle[1]) +
    (polygon$y - box$cy[1]) * cos(box$angle[1])
  any(abs(local_x) <= box$ow[1] / 2 + tol &
      abs(local_y) <= box$oh[1] / 2 + tol)
}

ggchord_label_hits_features <- function(box, polygons, source_row = NA_integer_,
                                         allow_own = FALSE, tol = .004) {
  for (polygon in polygons) {
    if (allow_own && !is.na(source_row) &&
        isTRUE(polygon$source_row[1L] == source_row)) next
    if (ggchord_box_overlaps_feature(box, polygon, tol = tol)) return(TRUE)
  }
  FALSE
}

#' Convert physical text dimensions to the current fixed-aspect plot scale.
#'
#' Text is rendered in millimetres, whereas chord geometry is expressed in
#' data units. A fixed data-units-per-inch constant therefore cannot describe
#' the same label on both a small and a large output device. This helper uses
#' the current device dimensions and the undecorated chord span; importantly,
#' it does not feed already-expanded label limits back into the estimate. The
#' latter used to make leader-line clipping grow with the labels themselves and
#' produced conspicuously large, output-size-dependent gaps.
#' @noRd
ggchord_text_obstacle_boxes <- function(seq_labels_df = NULL,
                                        axis_ticks = NULL,
                                        show_axis = FALSE,
                                        units_per_inch = 0.35,
                                        box_padding = 0.05) {
  out <- list()

  if (!is.null(seq_labels_df) && nrow(seq_labels_df) > 0) {
    out[[length(out) + 1]] <- ggchord_text_boxes(
      seq_labels_df,
      x_col = "text_x", y_col = "text_y", text_col = "label",
      angle_col = "text_angle", size_col = "size",
      hjust_col = "hjust", vjust_col = "vjust",
      units_per_inch = units_per_inch, box_padding = box_padding
    )
  }

  if (isTRUE(show_axis) && !is.null(axis_ticks) && nrow(axis_ticks) > 0) {
    axis_labels <- axis_ticks[!is.na(axis_ticks$label), , drop = FALSE]
    if (nrow(axis_labels) > 0) {
      out[[length(out) + 1]] <- ggchord_text_boxes(
        axis_labels,
        x_col = "label_x", y_col = "label_y", text_col = "label",
        angle_col = "label_angle", size_col = "size",
        hjust_col = "label_hjust", vjust_col = "label_vjust",
        units_per_inch = units_per_inch, box_padding = box_padding
      )
    }
  }

  if (length(out) == 0) {
    return(ggchord_text_boxes(data.frame()))
  }
  do.call(rbind, out)
}


# Return TRUE when two leader-line segments cross in their interiors.
ggchord_segments_cross <- function(ax, ay, bx, by, cx, cy, dx, dy,
                                   tol = 1e-10) {
  orient <- function(px, py, qx, qy, rx, ry) {
    (qx - px) * (ry - py) - (qy - py) * (rx - px)
  }
  o1 <- orient(ax, ay, bx, by, cx, cy)
  o2 <- orient(ax, ay, bx, by, dx, dy)
  o3 <- orient(cx, cy, dx, dy, ax, ay)
  o4 <- orient(cx, cy, dx, dy, bx, by)
  o1 * o2 < -tol & o3 * o4 < -tol
}

# Remove crossings by swapping label positions within each sequence-side lane.
#
# For two crossing segments, exchanging their endpoints strictly shortens the
# total segment length. Repeating that 2-opt operation therefore converges and
# preserves the set of label positions; layouts without crossings are left
# byte-for-byte unchanged.
