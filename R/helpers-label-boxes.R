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
                                   fit = c("wrap", "none", "ellipsis", "auto"),
                                   max_lines = 2L,
                                   units_per_inch = 0.35,
                                   box_padding = 0.18,
                                   repel_boxes = NULL) {
  fit <- match.arg(fit)
  if (identical(fit, "none") || nrow(gl) == 0L) return(gl)
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
                               units_per_inch = 0.35, box_padding = 0) {
  n <- nrow(df)
  empty <- data.frame(
    x = numeric(0), y = numeric(0),
    cx = numeric(0), cy = numeric(0),
    w = numeric(0), h = numeric(0),
    ow = numeric(0), oh = numeric(0), angle = numeric(0),
    bw = numeric(0), bh = numeric(0),
    xmin = numeric(0), xmax = numeric(0),
    ymin = numeric(0), ymax = numeric(0),
    stringsAsFactors = FALSE
  )
  if (n == 0) return(empty)

  x <- df[[x_col]]
  y <- df[[y_col]]
  texts <- df[[text_col]]
  sizes <- if (is.null(df[[size_col]])) rep(2.5, n) else df[[size_col]]
  angles <- (if (is.null(df[[angle_col]])) rep(0, n) else df[[angle_col]]) *
    pi / 180
  hjust <- if (is.null(df[[hjust_col]])) rep(0.5, n) else df[[hjust_col]]
  vjust <- if (is.null(df[[vjust_col]])) rep(0.5, n) else df[[vjust_col]]

  w <- numeric(n)
  h <- numeric(n)
  valid <- !is.na(texts) & nzchar(texts)
  if (any(valid)) {
    close_device <- ggchord_measurement_device()
    on.exit(close_device())
    # ggplot2 text sizes are millimetres and are converted to grid font points
    # with `.pt` (72.27 / 25.4). Base graphics' cex is relative to the
    # device's 12-point default. Omitting this conversion underestimates both
    # dimensions by about 2.845 and lets visibly overlapping labels pass the
    # collision test.
    text_cex <- sizes[valid] * (72.27 / 25.4) / 12
    w[valid] <- suppressWarnings(graphics::strwidth(
      texts[valid], units = "inches", cex = text_cex
    )) * units_per_inch
    n_lines <- vapply(strsplit(texts[valid], "\n"), length, integer(1))
    h[valid] <- suppressWarnings(graphics::strheight(
      texts[valid], units = "inches", cex = text_cex
    )) * n_lines * units_per_inch
  }

  cos_a <- cos(angles)
  sin_a <- sin(angles)
  cx_off <- (0.5 - hjust) * w * cos_a - (0.5 - vjust) * h * sin_a
  cy_off <- (0.5 - hjust) * w * sin_a + (0.5 - vjust) * h * cos_a
  bw <- abs(cos_a) * w + abs(sin_a) * h + 2 * box_padding * units_per_inch
  bh <- abs(sin_a) * w + abs(cos_a) * h + 2 * box_padding * units_per_inch

  data.frame(
    x = x, y = y,
    cx = x + cx_off, cy = y + cy_off,
    w = w, h = h,
    ow = w + 2 * box_padding * units_per_inch,
    oh = h + 2 * box_padding * units_per_inch,
    angle = angles,
    bw = bw, bh = bh,
    xmin = x + cx_off - bw / 2, xmax = x + cx_off + bw / 2,
    ymin = y + cy_off - bh / 2, ymax = y + cy_off + bh / 2,
    stringsAsFactors = FALSE
  )
}

# Test one oriented text rectangle against zero or more oriented rectangles.
# A cheap axis-aligned prefilter is followed by the separating-axis theorem;
# this avoids treating a diagonal label's large empty corner triangles as
# occupied space.
ggchord_oriented_box_overlaps <- function(candidate, other, tol = 1e-7) {
  if (is.null(other) || nrow(other) == 0) return(logical(0))
  possible <- candidate$xmin < other$xmax - tol &
    candidate$xmax > other$xmin + tol &
    candidate$ymin < other$ymax - tol &
    candidate$ymax > other$ymin + tol
  out <- rep(FALSE, nrow(other))
  rows <- which(possible)
  if (length(rows) == 0) return(out)

  candidate_angle <- candidate$angle[1] %||% 0
  candidate_ow <- candidate$ow[1] %||% candidate$bw[1]
  candidate_oh <- candidate$oh[1] %||% candidate$bh[1]
  candidate_axes <- rbind(
    c(cos(candidate_angle), sin(candidate_angle)),
    c(-sin(candidate_angle), cos(candidate_angle))
  )
  for (j in rows) {
    other_angle <- other$angle[j] %||% 0
    other_ow <- other$ow[j] %||% other$bw[j]
    other_oh <- other$oh[j] %||% other$bh[j]
    other_axes <- rbind(
      c(cos(other_angle), sin(other_angle)),
      c(-sin(other_angle), cos(other_angle))
    )
    axes <- rbind(candidate_axes, other_axes)
    centre_delta <- c(other$cx[j] - candidate$cx[1],
                      other$cy[j] - candidate$cy[1])
    separated <- FALSE
    for (k in seq_len(nrow(axes))) {
      axis <- axes[k, ]
      candidate_extent <-
        candidate_ow / 2 * abs(sum(axis * candidate_axes[1, ])) +
        candidate_oh / 2 * abs(sum(axis * candidate_axes[2, ]))
      other_extent <-
        other_ow / 2 * abs(sum(axis * other_axes[1, ])) +
        other_oh / 2 * abs(sum(axis * other_axes[2, ]))
      if (abs(sum(centre_delta * axis)) >=
          candidate_extent + other_extent - tol) {
        separated <- TRUE
        break
      }
    }
    out[j] <- !separated
  }
  out
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
