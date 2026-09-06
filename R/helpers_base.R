#' Colour-vision-friendly categorical palette
#' @noRd
chord_palette_discrete <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#56B4E9", "#E69F00", "#6F6F6F", "#F0E442"
)

#' Generate a default categorical palette
#'
#' Returns a colour-vision-friendly discrete palette. Larger palettes use a
#' qualitative HCL palette rather than interpolating through muddy midpoints.
#' @param n Number of colors requested
#' @return A character vector of \code{n} colors
#' @noRd
chord_default_palette <- function(n) {
  if (n <= 0) {
    character(0)
  } else if (n <= length(chord_palette_discrete)) {
    chord_palette_discrete[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, palette = "Dark 3")
  }
}

#' Missing value handling operator
#'
#' Used to safely handle NULL values: returns y if x is NULL, otherwise returns x
#'
#' @param x Any R object (may be NULL)
#' @param y Default value to return when x is NULL
#' @return x if x is not NULL, otherwise y
#' @keywords internal
if_null_else <- function (x, y)
{
  if (is.null(x)) y else x
}


#' Calculate plot extremes
#'
#' Extracts x/y coordinate extremes from all plot elements (sequence arcs, ribbons, gene arrows, etc.) for adjusting the plot range
#'
#' @param allRibbon data.frame, ribbon data (with x, y columns), default NULL
#' @param seqArcs List, sequence arc data (each element is a data frame with x, y, accver), default NULL
#' @param axisLines data.frame, axis line data (with x, y, accver columns), default NULL
#' @param axisTicks data.frame, tick mark data (with x0, y0, x1, y1, label_x, label_y columns), default NULL
#' @param gene_arrows data.frame, gene label data (with text_x, text_y columns), default NULL
#' @param gene_polys data.frame, gene arrow polygon data (with x, y columns), default NULL
#' @param show_axis Logical, whether to include extreme value calculation for axis-related elements, default FALSE
#' @return List containing x_min (minimum x), x_max (maximum x), y_min (minimum y), y_max (maximum y)
#' @keywords internal
get_plot_extremes <- function(allRibbon = NULL, seqArcs = NULL,
                                axisLines = NULL, axisTicks = NULL,
                                gene_arrows = NULL, gene_polys = NULL,
                                seq_labels = NULL, show_axis = FALSE) {
  x_min <- Inf
  x_max <- -Inf
  y_min <- Inf
  y_max <- -Inf

  include <- function(x, y) {
    if (length(x) == 0L || length(y) == 0L) return(invisible(NULL))
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
      x_min <<- min(x_min, min(x[ok]))
      x_max <<- max(x_max, max(x[ok]))
      y_min <<- min(y_min, min(y[ok]))
      y_max <<- max(y_max, max(y[ok]))
    }
    invisible(NULL)
  }

  if (!is.null(allRibbon) && nrow(allRibbon) > 0) {
    include(allRibbon$x, allRibbon$y)
  }
  if (!is.null(seqArcs) && length(seqArcs) > 0) {
    for (arc in seqArcs) {
      if (nrow(arc) > 0) include(arc$x, arc$y)
    }
  }
  if (!is.null(gene_arrows) && nrow(gene_arrows) > 0) {
    include(gene_arrows$text_x, gene_arrows$text_y)
  }
  if (show_axis && !is.null(axisLines) && nrow(axisLines) > 0) {
    include(axisLines$x, axisLines$y)
  }
  if (show_axis && !is.null(axisTicks) && nrow(axisTicks) > 0) {
    include(axisTicks$x0, axisTicks$y0)
    include(axisTicks$x1, axisTicks$y1)
    include(axisTicks$label_x, axisTicks$label_y)
  }
  if (!is.null(gene_polys) && nrow(gene_polys) > 0) {
    include(gene_polys$x, gene_polys$y)
  }
  if (!is.null(seq_labels) && nrow(seq_labels) > 0) {
    include(seq_labels$text_x, seq_labels$text_y)
  }

  list(
    x_min = if (is.finite(x_min)) x_min else NA_real_,
    x_max = if (is.finite(x_max)) x_max else NA_real_,
    y_min = if (is.finite(y_min)) y_min else NA_real_,
    y_max = if (is.finite(y_max)) y_max else NA_real_
  )
}

#' Wrap long gene annotation texts at a given character width
#' @keywords internal
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
#' @keywords internal
ggchord_label_deoverlap <- function(gl, units_per_inch = 0.35, seed = 123,
                                    max_overlaps = Inf) {
  if (nrow(gl) < 2) return(gl)

  close_device <- ggchord_measurement_device()
  on.exit(close_device())
  sizes <- gl$size %||% rep(2.5, nrow(gl))
  w <- suppressWarnings(graphics::strwidth(gl$text, units = "inches",
                                           cex = sizes / 12)) * units_per_inch
  n_lines <- vapply(strsplit(gl$text, "\n"), length, integer(1))
  h <- suppressWarnings(graphics::strheight(gl$text, units = "inches",
                                            cex = sizes / 12)) *
    n_lines * units_per_inch

  x <- gl$text_x
  y <- gl$text_y
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

  gl$text_x <- x
  gl$text_y <- y
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
#' @keywords internal
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
#' @keywords internal
ggchord_device_units_per_inch <- function(x, y,
                                          fallback_inches = 6,
                                          margin_inches = 1.25,
                                          device_inches = NULL) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  x_span <- if (length(x) > 1) diff(range(x)) else 0
  y_span <- if (length(y) > 1) diff(range(y)) else 0
  geometry_span <- max(x_span, y_span, 1)

  if (is.null(device_inches)) {
    # Querying dev.size() on the null device opens R's default device. That
    # side effect could leave Rplots.pdf open after layout-only operations and
    # make a following ggsave() close the wrong device. Use the documented
    # fallback until an actual render device exists.
    device_inches <- if (grDevices::dev.cur() == 1L) {
      c(NA_real_, NA_real_)
    } else tryCatch(
      grDevices::dev.size("in"), error = function(e) c(NA_real_, NA_real_)
    )
  }
  device_short_side <- suppressWarnings(min(device_inches, na.rm = TRUE))
  if (!is.finite(device_short_side) || device_short_side <= 0) {
    usable_inches <- fallback_inches
  } else {
    # Keep a modest allowance for titles and legends, but never replace a
    # genuinely small device with the much larger fallback canvas. The old
    # `< 2` fallback made 4 x 3 inch exports behave as if they were six inches
    # wide and consequently underestimated every text box.
    reserve <- min(margin_inches, device_short_side * 0.4)
    usable_inches <- max(device_short_side - reserve, 0.75)
  }
  geometry_span / usable_inches
}

#' Convert text layers into fixed obstacle rectangles for label repulsion.
#' @keywords internal
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
ggchord_uncross_labels <- function(gl, lanes = NULL, max_swaps = NULL,
                                   endpoint_fun = NULL) {
  n <- nrow(gl)
  if (n < 2) return(list(labels = gl, swaps = 0L))
  if (is.null(lanes)) lanes <- gl$accver %||% rep("all", n)
  if (length(lanes) != n) {
    ggchord_stop("lanes must have one value per gene label")
  }

  active <- !is.na(gl$text) & nzchar(gl$text)
  lane_rows <- split(which(active), as.character(lanes[active]), drop = TRUE)
  total_swaps <- 0L

  for (rows in lane_rows) {
    nr <- length(rows)
    if (nr < 2) next
    swap_limit <- max_swaps %||% max(100L, 4L * nr^2)
    changed <- TRUE
    lane_swaps <- 0L

    while (changed && lane_swaps < swap_limit) {
      changed <- FALSE
      endpoints <- if (is.null(endpoint_fun)) {
        data.frame(x = gl$text_x, y = gl$text_y)
      } else {
        endpoint_fun(gl)
      }
      for (ii in seq_len(nr - 1L)) {
        for (jj in (ii + 1L):nr) {
          i <- rows[ii]
          j <- rows[jj]
          crossed <- ggchord_segments_cross(
            gl$anchor_x[i], gl$anchor_y[i], endpoints$x[i], endpoints$y[i],
            gl$anchor_x[j], gl$anchor_y[j], endpoints$x[j], endpoints$y[j]
          )
          if (!crossed) next

          if (is.null(endpoint_fun)) {
            candidate_x_i <- gl$text_x[j]
            candidate_y_i <- gl$text_y[j]
            candidate_x_j <- gl$text_x[i]
            candidate_y_j <- gl$text_y[i]
            candidate_endpoints <- NULL
          } else {
            candidate <- gl
            candidate$text_x[c(i, j)] <- gl$text_x[c(j, i)]
            candidate$text_y[c(i, j)] <- gl$text_y[c(j, i)]
            candidate_endpoints <- endpoint_fun(candidate)
          }
          old_length <- sqrt((endpoints$x[i] - gl$anchor_x[i])^2 +
                               (endpoints$y[i] - gl$anchor_y[i])^2) +
            sqrt((endpoints$x[j] - gl$anchor_x[j])^2 +
                   (endpoints$y[j] - gl$anchor_y[j])^2)
          new_length <- if (is.null(endpoint_fun)) {
            sqrt((candidate_x_i - gl$anchor_x[i])^2 +
                   (candidate_y_i - gl$anchor_y[i])^2) +
              sqrt((candidate_x_j - gl$anchor_x[j])^2 +
                     (candidate_y_j - gl$anchor_y[j])^2)
          } else {
            sqrt((candidate_endpoints$x[i] - gl$anchor_x[i])^2 +
                   (candidate_endpoints$y[i] - gl$anchor_y[i])^2) +
              sqrt((candidate_endpoints$x[j] - gl$anchor_x[j])^2 +
                     (candidate_endpoints$y[j] - gl$anchor_y[j])^2)
          }
          if (new_length >= old_length - 1e-10) next

          if (is.null(endpoint_fun)) {
            gl$text_x[c(i, j)] <- c(candidate_x_i, candidate_x_j)
            gl$text_y[c(i, j)] <- c(candidate_y_i, candidate_y_j)
          } else {
            gl <- candidate
          }
          lane_swaps <- lane_swaps + 1L
          total_swaps <- total_swaps + 1L
          changed <- TRUE
          break
        }
        if (changed) break
      }
    }
  }

  list(labels = gl, swaps = total_swaps)
}

ggchord_elbow_bends <- function(gl, text_widths, text_heights = NULL,
                                directions = NULL) {
  n <- nrow(gl)
  if (n == 0) return(data.frame(x = numeric(0), y = numeric(0)))
  if (is.null(text_heights)) text_heights <- text_widths
  if (is.null(directions)) {
    directions <- ifelse(gl$hjust < 0.5, "right", "left")
  }
  horiz <- abs(gl$text_x - gl$anchor_x)
  vert <- abs(gl$text_y - gl$anchor_y)
  horizontal_stub <- directions %in% c("left", "right")
  span <- ifelse(horizontal_stub, horiz, vert)
  extent <- ifelse(horizontal_stub, text_widths, text_heights)
  stub_len <- pmin(pmax(0.02, 0.3 * span), pmax(0.3 * extent, 0.04))

  bx <- gl$text_x
  by <- gl$text_y
  left <- directions == "left"
  right <- directions == "right"
  bottom <- directions == "bottom"
  top <- directions == "top"
  bx[left] <- pmin(gl$text_x[left] + stub_len[left], gl$anchor_x[left])
  bx[right] <- pmax(gl$text_x[right] - stub_len[right], gl$anchor_x[right])
  by[bottom] <- pmin(gl$text_y[bottom] + stub_len[bottom], gl$anchor_y[bottom])
  by[top] <- pmax(gl$text_y[top] - stub_len[top], gl$anchor_y[top])
  data.frame(x = bx, y = by)
}

ggchord_label_curve_frame <- function(gl, seq_arcs) {
  n <- nrow(gl)
  frame <- data.frame(
    curve_x = numeric(n), curve_y = numeric(n),
    outward_x = numeric(n), outward_y = numeric(n),
    signed_distance = numeric(n),
    curve_index = integer(n)
  )
  if (n == 0) return(frame)

  arc_ids <- vapply(seq_arcs, function(a) {
    if (nrow(a) == 0) "" else as.character(unique(a$accver)[1])
  }, character(1))

  for (sid in unique(gl$accver)) {
    rows <- which(gl$accver == sid)
    arc_pos <- match(sid, arc_ids)
    if (is.na(arc_pos) || nrow(seq_arcs[[arc_pos]]) < 2) next
    arc <- seq_arcs[[arc_pos]]

    for (i in rows) {
      # Locate the label by its fixed gene anchor, not by the repelled text.
      # This prevents a far-moving label from snapping to another part of a
      # highly curved sequence path.
      d2 <- (arc$x - gl$anchor_x[i])^2 + (arc$y - gl$anchor_y[i])^2
      k <- which.min(d2)
      k0 <- max(1L, k - 1L)
      k1 <- min(nrow(arc), k + 1L)
      tx <- arc$x[k1] - arc$x[k0]
      ty <- arc$y[k1] - arc$y[k0]
      tangent_length <- sqrt(tx^2 + ty^2)
      if (!is.finite(tangent_length) || tangent_length < 1e-12) next

      # Reference paths follow increasing genomic angle. Their right normal
      # is the outside track, including concave curves and paths crossing the
      # origin; an origin-dot-product test would flip sides mid-sequence.
      path_direction <- attr(arc, "ggchord_path_direction") %||% 1
      nx <- path_direction * ty / tangent_length
      ny <- -path_direction * tx / tangent_length
      frame$curve_x[i] <- arc$x[k]
      frame$curve_y[i] <- arc$y[k]
      frame$outward_x[i] <- nx
      frame$outward_y[i] <- ny
      frame$curve_index[i] <- k
      frame$signed_distance[i] <-
        (gl$text_x[i] - arc$x[k]) * nx +
        (gl$text_y[i] - arc$y[k]) * ny
    }
  }
  frame
}

ggchord_enforce_label_side <- function(gl, seq_arcs, side = "auto") {
  if (nrow(gl) == 0 || identical(side, "auto")) return(gl)
  frame <- ggchord_label_curve_frame(gl, seq_arcs)
  want_inside <- identical(side, "inside")
  flip <- (frame$signed_distance < 0) != want_inside
  flip[!is.finite(frame$signed_distance)] <- FALSE
  if (!any(flip)) return(gl)

  # Reflect across the actual local tangent of the sequence curve. Unlike a
  # radius-from-origin approximation, this remains correct when seq_radius,
  # seq_curvature, seq_gap, seq_order or global rotation change the geometry.
  correction <- 2 * frame$signed_distance[flip]
  gl$text_x[flip] <- gl$text_x[flip] -
    correction * frame$outward_x[flip]
  gl$text_y[flip] <- gl$text_y[flip] -
    correction * frame$outward_y[flip]
  gl
}

# Pack one-dimensional label anchors while preserving their gene order.
ggchord_pack_label_axis <- function(preferred, order_value,
                                    before, after, gap = 0) {
  n <- length(preferred)
  if (n < 2) return(preferred)
  ord <- order(order_value, seq_len(n))
  position <- preferred[ord]
  before <- before[ord]
  after <- after[ord]

  # A forward pass creates the minimum legal spacing. Translating the whole
  # lane afterwards is the least-squares fit back to the preferred positions
  # and does not change any of those spacings.
  for (i in 2:n) {
    position[i] <- max(
      position[i], position[i - 1L] + after[i - 1L] + before[i] + gap
    )
  }
  position <- position + mean(preferred[ord] - position)
  out <- numeric(n)
  out[ord] <- position
  out
}

# Pack away from a shared lane centre. The label nearest that centre remains
# close to its gene, while congestion is absorbed by labels farther towards
# either end of the sequence. This avoids giving adjacent labels parallel,
# similarly offset leaders that can intersect on a curved arc.
ggchord_pack_label_axis_outward <- function(preferred, order_value,
                                            before, after, centre,
                                            gap = 0) {
  n <- length(preferred)
  if (n < 2) return(preferred)
  ord <- order(order_value, seq_len(n))
  position <- preferred[ord]
  before <- before[ord]
  after <- after[ord]
  pivot <- which.min(abs(order_value[ord] - centre))

  if (pivot > 1) {
    for (i in seq.int(pivot - 1L, 1L)) {
      position[i] <- min(
        position[i], position[i + 1L] - after[i] - before[i + 1L] - gap
      )
    }
  }
  if (pivot < n) {
    for (i in seq.int(pivot + 1L, n)) {
      position[i] <- max(
        position[i], position[i - 1L] + after[i - 1L] + before[i] + gap
      )
    }
  }
  out <- numeric(n)
  out[ord] <- position
  out
}

# Put horizontal labels on compact cardinal lanes around their own sequence.
# Shared left/right columns used by the auto layout.
ggchord_side_label_columns <- function(gl, seq_arcs,
                                        side = "outside",
                                        units_per_inch = 0.35,
                                        box_padding = 0.25,
                                        point_padding = 0.1,
                                        repel_boxes = NULL,
                                        max_iter = 100, directions = NULL) {
  n <- nrow(gl)
  if (n == 0) return(list(labels = gl, lanes = character(0)))

  active <- !is.na(gl$text) & nzchar(gl$text)
  source_frame <- ggchord_label_curve_frame(gl, seq_arcs)
  anchor_labels <- gl
  anchor_labels$text_x <- anchor_labels$anchor_x
  anchor_labels$text_y <- anchor_labels$anchor_y
  anchor_frame <- ggchord_label_curve_frame(anchor_labels, seq_arcs)
  side_sign <- if (identical(side, "outside")) {
    rep(1, n)
  } else if (identical(side, "inside")) {
    rep(-1, n)
  } else {
    ifelse(source_frame$signed_distance < 0, -1, 1)
  }
  side_sign[!is.finite(side_sign)] <- 1
  desired_x <- side_sign * anchor_frame$outward_x
  desired_y <- side_sign * anchor_frame$outward_y
  # Give every sequence-side group one primary cardinal rail. This keeps the
  # labels belonging to one sequence visually close and prevents neighbouring
  # sequences from competing for the same corners. Parallel rows absorb genuine
  # congestion without treating unused sides as a target to distribute toward.
  if (is.null(directions)) directions <- ifelse(desired_x < 0, "left", "right")
  directions[!active] <- NA_character_
  base_lanes <- paste(gl$accver, directions, sep = "\r")

  gl$text_angle[active] <- 0
  gl$hjust[active] <- c(left = 1, right = 0, top = 0.5, bottom = 0.5)[
    directions[active]
  ]
  gl$vjust[active] <- c(left = 0.5, right = 0.5, top = 0, bottom = 1)[
    directions[active]
  ]
  axis_gap <- max(0.02, 0.025 * units_per_inch)
  base_lane_rows <- split(which(active), base_lanes[active], drop = TRUE)
  rail_gap <- point_padding + box_padding * units_per_inch +
    max(0.04, 0.04 * units_per_inch)
  lanes <- base_lanes
  tracks <- rep(NA_integer_, n)

  # All labels assigned to the same side share one compact vertical column.
  # Packing from feature anchors, rather than previous-pass label positions,
  # avoids layout feedback, preserves their global vertical order and prevents
  # two sequence-specific columns from sending leaders across one another.
  vertical_rows_all <- which(active & directions %in% c("left", "right"))
  vertical_groups <- split(
    vertical_rows_all, directions[vertical_rows_all], drop = TRUE
  )
  template_boxes <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  for (rows in vertical_groups) {
    direction <- directions[rows[1]]
    before <- template_boxes$y[rows] - template_boxes$ymin[rows]
    after <- template_boxes$ymax[rows] - template_boxes$y[rows]
    gl$text_y[rows] <- ggchord_pack_label_axis(
      gl$anchor_y[rows], gl$anchor_y[rows], before, after, gap = axis_gap
    )
    gl$text_x[rows] <- if (direction == "left") {
      min(anchor_frame$curve_x[rows]) - rail_gap
    } else {
      max(anchor_frame$curve_x[rows]) + rail_gap
    }
    tracks[rows] <- 1L
    lanes[rows] <- direction
  }

  # A shared side column must also use a shared endpoint order. Curved and
  # rotated sequences can place two feature anchors in an order that differs
  # from their first packed label positions even when those positions do not
  # overlap. Swap only crossing endpoints, then repack the resulting slots for
  # their actual text heights. Repeating this small 2-opt pass produces a
  # deterministic, compact column whose direct leaders do not cross.
  if (length(vertical_rows_all) > 1L) {
    for (pass in seq_len(8L)) {
      candidate <- gl[vertical_rows_all, , drop = FALSE]
      uncrossed <- ggchord_uncross_labels(
        candidate, lanes = directions[vertical_rows_all]
      )
      if (uncrossed$swaps == 0L) break
      candidate <- uncrossed$labels
      candidate_boxes <- ggchord_text_boxes(
        candidate, units_per_inch = units_per_inch,
        box_padding = box_padding
      )
      candidate_groups <- split(
        seq_along(vertical_rows_all),
        directions[vertical_rows_all], drop = TRUE
      )
      for (local_rows in candidate_groups) {
        before <- candidate_boxes$y[local_rows] -
          candidate_boxes$ymin[local_rows]
        after <- candidate_boxes$ymax[local_rows] -
          candidate_boxes$y[local_rows]
        candidate$text_y[local_rows] <- ggchord_pack_label_axis(
          candidate$text_y[local_rows], candidate$text_y[local_rows],
          before, after, gap = axis_gap
        )
      }
      gl$text_x[vertical_rows_all] <- candidate$text_x
      gl$text_y[vertical_rows_all] <- candidate$text_y
    }
  }

  move_lane_outward <- function(rows, amount) {
    if (!is.finite(amount) || amount <= 0) return()
    direction <- directions[rows[1]]
    if (direction == "left") gl$text_x[rows] <<- gl$text_x[rows] - amount
    if (direction == "right") gl$text_x[rows] <<- gl$text_x[rows] + amount
    if (direction == "bottom") gl$text_y[rows] <<- gl$text_y[rows] - amount
    if (direction == "top") gl$text_y[rows] <<- gl$text_y[rows] + amount
  }

  ensure_vertical_side <- function() {
    for (rows in vertical_groups) {
      direction <- directions[rows[1]]
      dx <- gl$text_x[rows] - anchor_frame$curve_x[rows]
      dy <- gl$text_y[rows] - anchor_frame$curve_y[rows]
      signed <- side_sign[rows] *
        (dx * anchor_frame$outward_x[rows] +
           dy * anchor_frame$outward_y[rows])
      projection <- switch(
        direction,
        left = -side_sign[rows] * anchor_frame$outward_x[rows],
        right = side_sign[rows] * anchor_frame$outward_x[rows],
        bottom = -side_sign[rows] * anchor_frame$outward_y[rows],
        top = side_sign[rows] * anchor_frame$outward_y[rows]
      )
      usable <- is.finite(projection) & projection > 0.25
      if (!any(usable)) next
      target <- max(0.015, 0.025 * units_per_inch)
      amount <- max((target - signed[usable]) / projection[usable], 0)
      # The rail coordinates already place it beyond the relevant curve
      # extreme. This cap is only a local-side correction, never a second
      # layout offset.
      move_lane_outward(rows, min(amount, rail_gap))
    }
  }
  ensure_vertical_side()

  # Sequence and axis text are fixed obstacles. Resolve those conflicts by
  # moving the complete rail outwards, which preserves its alignment and the
  # order of every leader instead of pushing individual labels off the rail.
  if (!is.null(repel_boxes) && nrow(repel_boxes) > 0) {
    for (pass in seq_len(8)) {
      moved <- FALSE
      boxes <- ggchord_text_boxes(
        gl, units_per_inch = units_per_inch, box_padding = 0
      )
      for (rows in vertical_groups) {
        overlap_x <- outer(
          boxes$xmin[rows], repel_boxes$xmax,
          function(a, b) a < b - 1e-7
        ) & outer(
          boxes$xmax[rows], repel_boxes$xmin,
          function(a, b) a > b + 1e-7
        )
        overlap_y <- outer(
          boxes$ymin[rows], repel_boxes$ymax,
          function(a, b) a < b - 1e-7
        ) & outer(
          boxes$ymax[rows], repel_boxes$ymin,
          function(a, b) a > b + 1e-7
        )
        hits <- which(overlap_x & overlap_y, arr.ind = TRUE)
        if (nrow(hits) == 0) next
        direction <- directions[rows[1]]
        amount <- switch(
          direction,
          left = max(boxes$xmax[rows[hits[, 1]]] -
                       repel_boxes$xmin[hits[, 2]] + axis_gap),
          right = max(repel_boxes$xmax[hits[, 2]] -
                        boxes$xmin[rows[hits[, 1]]] + axis_gap),
          bottom = max(boxes$ymax[rows[hits[, 1]]] -
                         repel_boxes$ymin[hits[, 2]] + axis_gap),
          top = max(repel_boxes$ymax[hits[, 2]] -
                      boxes$ymin[rows[hits[, 1]]] + axis_gap)
        )
        move_lane_outward(rows, amount)
        moved <- TRUE
      }
      if (!moved) break
    }
  }
  ensure_vertical_side()


  list(
    labels = gl,
    lanes = lanes,
    directions = directions,
    tracks = tracks
  )
}

# Put labels on the nearest collision-free local offset track. Unlike a
# radius-from-origin layout, every candidate position is measured from the
# actual sequence curve and its local outward normal, so straight, strongly
# curved and differently sized sequences use the same algorithm.
ggchord_offset_label_tracks <- function(gl, seq_arcs,
                                        side = "outside",
                                        orientation = c("horizontal", "arc"),
                                        units_per_inch = 0.35,
                                        box_padding = 0.18,
                                        point_padding = 0.08,
                                        repel_boxes = NULL) {
  orientation <- match.arg(orientation)
  n <- nrow(gl)
  if (n == 0) {
    return(list(labels = gl, lanes = character(0),
                directions = character(0), tracks = integer(0),
                draw_segment = logical(0)))
  }

  active <- !is.na(gl$text) & nzchar(gl$text)
  source_frame <- ggchord_label_curve_frame(gl, seq_arcs)
  anchor_labels <- gl
  anchor_labels$text_x <- anchor_labels$anchor_x
  anchor_labels$text_y <- anchor_labels$anchor_y
  frame <- ggchord_label_curve_frame(anchor_labels, seq_arcs)
  side_sign <- if (identical(side, "outside")) {
    rep(1, n)
  } else if (identical(side, "inside")) {
    rep(-1, n)
  } else {
    ifelse(source_frame$signed_distance < 0, -1, 1)
  }
  side_sign[!is.finite(side_sign)] <- 1
  normal_x <- side_sign * frame$outward_x
  normal_y <- side_sign * frame$outward_y
  directions <- ifelse(
    abs(normal_x) >= abs(normal_y),
    ifelse(normal_x < 0, "left", "right"),
    ifelse(normal_y < 0, "bottom", "top")
  )

  if (identical(orientation, "horizontal")) {
    gl$text_angle[active] <- 0
    gl$hjust[active] <- c(left = 1, right = 0, top = 0.5, bottom = 0.5)[
      directions[active]
    ]
    gl$vjust[active] <- c(left = 0.5, right = 0.5, top = 0, bottom = 1)[
      directions[active]
    ]
  } else {
    tangent_x <- -frame$outward_y
    tangent_y <- frame$outward_x
    angle <- (atan2(tangent_y, tangent_x) * 180 / pi + 360) %% 360
    upside_down <- angle > 90 & angle < 270
    angle[upside_down] <- (angle[upside_down] + 180) %% 360
    gl$text_angle[active] <- angle[active]
    gl$hjust[active] <- 0.5
    gl$vjust[active] <- 0.5
  }

  # Measure the anchor-relative box offset once. The innermost edge of every
  # first-track label is then placed at the same visual clearance from its
  # own sequence curve, even when text justification differs by quadrant.
  gl$text_x <- frame$curve_x
  gl$text_y <- frame$curve_y
  templates <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  centre_offset <-
    (templates$cx - templates$x) * normal_x +
    (templates$cy - templates$y) * normal_y
  text_angle <- gl$text_angle * pi / 180
  text_x_axis_x <- cos(text_angle)
  text_x_axis_y <- sin(text_angle)
  text_y_axis_x <- -sin(text_angle)
  text_y_axis_y <- cos(text_angle)
  # Project the oriented rectangle itself, not its axis-aligned bounding box.
  # Projecting the latter makes a long tangent-aligned arc label look long in
  # the normal direction as well and can push it several unnecessary tracks
  # away from the sequence.
  normal_extent <-
    templates$w * abs(text_x_axis_x * normal_x +
                        text_x_axis_y * normal_y) / 2 +
    templates$h * abs(text_y_axis_x * normal_x +
                        text_y_axis_y * normal_y) / 2 +
    box_padding * units_per_inch
  clearance <- point_padding + max(0.035, 0.04 * units_per_inch)
  base_distance <- pmax(
    clearance - centre_offset + normal_extent,
    clearance
  )

  base_lanes <- paste(gl$accver, side_sign, sep = "\r")
  lane_rows <- split(which(active), base_lanes[active], drop = TRUE)
  tracks <- rep(NA_integer_, n)
  placed_boxes <- NULL
  axis_gap <- max(0.015, 0.02 * units_per_inch)

  overlaps_boxes <- function(candidate, other) {
    any(ggchord_oriented_box_overlaps(candidate, other))
  }

  for (rows in lane_rows) {
    rows <- rows[order(frame$curve_index[rows], rows)]
    lane_step <- max(
      2 * normal_extent[rows] + axis_gap,
      0.08 + 0.08 * units_per_inch,
      na.rm = TRUE
    )
    for (i in rows) {
      selected <- FALSE
      # At most one new track per active label is needed when boxes are
      # finite, with a small reserve for fixed sequence/axis obstacles.
      for (track in seq_len(length(rows) + 8L)) {
        distance <- base_distance[i] + (track - 1L) * lane_step
        # Keep the text centre on its feature's local normal. Tangential
        # nudging can reverse two neighbouring labels and consequently make
        # their leaders cross, even when both text boxes remain disjoint.
        candidate_label <- gl[i, , drop = FALSE]
        candidate_label$text_x <- frame$curve_x[i] + normal_x[i] * distance
        candidate_label$text_y <- frame$curve_y[i] + normal_y[i] * distance
        candidate_box <- ggchord_text_boxes(
          candidate_label,
          units_per_inch = units_per_inch,
          box_padding = box_padding
        )
        if (!overlaps_boxes(candidate_box, placed_boxes) &&
            !overlaps_boxes(candidate_box, repel_boxes)) {
          gl$text_x[i] <- candidate_label$text_x
          gl$text_y[i] <- candidate_label$text_y
          tracks[i] <- track
          placed_boxes <- rbind(placed_boxes, candidate_box)
          selected <- TRUE
        }
        if (selected) break
      }
      if (!selected) {
        track <- length(rows) + 9L
        distance <- base_distance[i] + (track - 1L) * lane_step
        gl$text_x[i] <- frame$curve_x[i] + normal_x[i] * distance
        gl$text_y[i] <- frame$curve_y[i] + normal_y[i] * distance
        tracks[i] <- track
        placed_boxes <- rbind(
          placed_boxes,
          ggchord_text_boxes(
            gl[i, , drop = FALSE], units_per_inch = units_per_inch,
            box_padding = box_padding
          )
        )
      }
    }
  }

  lanes <- paste(base_lanes, tracks, sep = "\r")
  draw_segment <- active
  if (identical(orientation, "arc")) {
    draw_segment <- active & !is.na(tracks) & tracks > 1L
  }
  list(
    labels = gl,
    lanes = lanes,
    directions = directions,
    tracks = tracks,
    draw_segment = draw_segment
  )
}

# Hide only labels that remain conflicted after a deterministic layout. The
# default max_overlaps = Inf therefore retains every label, while a finite
# value behaves as a final decluttering threshold without influencing any
# successfully placed label coordinates.
ggchord_hide_conflicted_labels <- function(gl, max_overlaps = Inf,
                                            units_per_inch = 0.35,
                                            repel_boxes = NULL) {
  if (!is.finite(max_overlaps) || nrow(gl) == 0) return(gl)
  counts <- ggchord_label_conflict_counts(
    gl, units_per_inch = units_per_inch, repel_boxes = repel_boxes
  )
  gl$text[counts > max_overlaps] <- NA_character_
  gl
}

ggchord_label_box_conflicts <- function(gl, units_per_inch = 0.35,
                                        box_padding = 0.25,
                                        repel_boxes = NULL,
                                        tol = 1e-7) {
  active <- !is.na(gl$text) & nzchar(gl$text)
  gl <- gl[active, , drop = FALSE]
  n <- nrow(gl)
  if (n == 0) return(FALSE)

  boxes <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  if (n > 1) {
    dx <- abs(outer(boxes$cx, boxes$cx, "-"))
    dy <- abs(outer(boxes$cy, boxes$cy, "-"))
    overlap <- upper.tri(dx) &
      dx < outer(boxes$bw, boxes$bw, "+") / 2 - tol &
      dy < outer(boxes$bh, boxes$bh, "+") / 2 - tol
    if (any(overlap)) return(TRUE)
  }

  if (!is.null(repel_boxes) && nrow(repel_boxes) > 0) {
    overlap_x <- outer(boxes$xmin, repel_boxes$xmax, function(a, b) a < b - tol) &
      outer(boxes$xmax, repel_boxes$xmin, function(a, b) a > b + tol)
    overlap_y <- outer(boxes$ymin, repel_boxes$ymax, function(a, b) a < b - tol) &
      outer(boxes$ymax, repel_boxes$ymin, function(a, b) a > b + tol)
    if (any(overlap_x & overlap_y)) return(TRUE)
  }

  FALSE
}

# Collapse only elbow stubs that take part in a crossing. This includes corner
# crossings between different cardinal rails; the corresponding leader becomes
# straight, while all conflict-free elbows retain their bend and stub lengths.
ggchord_collapse_crossed_elbows <- function(segments, lanes,
                                            max_passes = NULL) {
  if (nrow(segments) < 2 || length(lanes) == 0) return(segments)
  groups <- unique(segments$group)
  max_passes <- max_passes %||% length(groups)

  for (pass in seq_len(max_passes)) {
    bad <- integer(0)
    for (i in seq_len(nrow(segments) - 1L)) {
      gi <- segments$group[i]
      for (j in (i + 1L):nrow(segments)) {
        gj <- segments$group[j]
        if (gi == gj || is.na(lanes[gi]) || is.na(lanes[gj])) next
        if (ggchord_segments_cross(
          segments$x0[i], segments$y0[i], segments$x1[i], segments$y1[i],
          segments$x0[j], segments$y0[j], segments$x1[j], segments$y1[j]
        )) {
          bad <- c(bad, gi, gj)
        }
      }
    }
    bad <- unique(bad)
    if (length(bad) == 0) break

    changed <- FALSE
    for (g in bad) {
      rows <- which(segments$group == g)
      if (length(rows) < 2) next
      first <- rows[1]
      stub <- rows[length(rows)]
      if (segments$x0[stub] == segments$x1[stub] &&
          segments$y0[stub] == segments$y1[stub]) next
      # The final row ends at the text anchor. Move the preceding bend there
      # and collapse the horizontal stub to a zero-length segment.
      segments$x1[first] <- segments$x1[stub]
      segments$y1[first] <- segments$y1[stub]
      segments$x0[stub] <- segments$x1[stub]
      segments$y0[stub] <- segments$y1[stub]
      changed <- TRUE
    }
    if (!changed) break
  }

  segments
}

#' Rebuild straight leader-line segments after a final label de-overlap pass.
#' @keywords internal
ggchord_repel_segments <- function(gl, min_segment_length = 0.5) {
  n <- nrow(gl)
  empty <- data.frame(x0 = numeric(0), y0 = numeric(0),
                      x1 = numeric(0), y1 = numeric(0),
                      group = integer(0), stringsAsFactors = FALSE)
  if (n == 0) return(empty)

  if ("anchor_x" %in% names(gl)) {
    ax <- gl$anchor_x
    ay <- gl$anchor_y
  } else {
    ax <- gl$text_x
    ay <- gl$text_y
  }
  seg_dist <- sqrt((gl$text_x - ax)^2 + (gl$text_y - ay)^2)
  visible <- if ("text" %in% names(gl)) {
    !is.na(gl$text) & nzchar(gl$text)
  } else {
    rep(TRUE, n)
  }
  keep_seg <- seg_dist > min_segment_length & visible
  data.frame(
    x0 = ax[keep_seg], y0 = ay[keep_seg],
    x1 = gl$text_x[keep_seg], y1 = gl$text_y[keep_seg],
    group = which(keep_seg),
    stringsAsFactors = FALSE
  )
}

# Split leader lines at the real oriented rectangles of other labels. Covered
# pieces can be faded, clipped, or shown in full. The target label itself is a
# hard boundary when include_own is TRUE; it is never represented by a faded
# line running through its own text.
ggchord_clip_segments_to_labels <- function(segments, gl,
                                            units_per_inch = 0.35,
                                            padding = 0.01,
                                            include_own = FALSE,
                                            overlap = c("fade", "clip", "show"),
                                            overlap_alpha = 0.18) {
  overlap <- match.arg(overlap)
  if (nrow(segments) == 0) return(segments)
  if (!is.numeric(overlap_alpha) || length(overlap_alpha) != 1L ||
      !is.finite(overlap_alpha) || overlap_alpha < 0 || overlap_alpha > 1) {
    ggchord_stop("overlap_alpha must be one finite number in [0, 1]")
  }
  annotate <- function(x, occluded = FALSE, alpha = 1) {
    x$occluded <- rep(occluded, nrow(x))
    x$alpha <- rep(alpha, nrow(x))
    x
  }
  if (nrow(gl) == 0 || (!isTRUE(include_own) && nrow(gl) < 2)) {
    return(annotate(segments))
  }
  if (identical(overlap, "show") && !isTRUE(include_own)) {
    return(annotate(segments))
  }
  boxes <- ggchord_text_boxes(gl, units_per_inch = units_per_inch)
  pad <- padding * units_per_inch
  angles <- gl$text_angle %||% rep(0, nrow(gl))
  angles[!is.finite(angles)] <- 0
  visible <- !is.na(gl$text) & nzchar(gl$text)
  out <- list()

  inside_interval <- function(origin, delta, lower, upper, tol = 1e-10) {
    if (abs(delta) < tol) {
      if (origin <= lower || origin >= upper) return(NULL)
      return(c(-Inf, Inf))
    }
    sort(c((lower - origin) / delta, (upper - origin) / delta))
  }

  for (s in seq_len(nrow(segments))) {
    dx <- segments$x1[s] - segments$x0[s]
    dy <- segments$y1[s] - segments$y0[s]
    cuts <- data.frame(start = numeric(), end = numeric(), hard = logical())
    others <- if (identical(overlap, "show")) {
      intersect(which(visible), segments$group[s])
    } else {
      which(visible)
    }
    if (!isTRUE(include_own)) {
      others <- setdiff(others, segments$group[s])
    }
    for (i in others) {
      # Intersect in the label's own coordinate system. Using its
      # axis-aligned bounding box over-clips diagonal leaders around rotated
      # arc labels, producing device-size-dependent gaps that are visibly
      # wider than the text itself.
      angle <- angles[i] * pi / 180
      cos_a <- cos(angle)
      sin_a <- sin(angle)
      rel_x <- segments$x0[s] - boxes$cx[i]
      rel_y <- segments$y0[s] - boxes$cy[i]
      local_x0 <- rel_x * cos_a + rel_y * sin_a
      local_y0 <- -rel_x * sin_a + rel_y * cos_a
      local_dx <- dx * cos_a + dy * sin_a
      local_dy <- -dx * sin_a + dy * cos_a
      tx <- inside_interval(local_x0, local_dx,
                            -boxes$w[i] / 2 - pad,
                            boxes$w[i] / 2 + pad)
      ty <- inside_interval(local_y0, local_dy,
                            -boxes$h[i] / 2 - pad,
                            boxes$h[i] / 2 + pad)
      if (is.null(tx) || is.null(ty)) next
      cut_start <- max(0, tx[1], ty[1])
      cut_end <- min(1, tx[2], ty[2])
      if (cut_start >= cut_end - 1e-10) next
      cuts <- rbind(
        cuts,
        data.frame(
          start = cut_start, end = cut_end,
          hard = isTRUE(include_own) && i == segments$group[s]
        )
      )
    }

    breaks <- sort(unique(c(0, 1, cuts$start, cuts$end)))
    if (length(breaks) < 2L) next
    for (p in seq_len(length(breaks) - 1L)) {
      start <- breaks[p]
      end <- breaks[p + 1L]
      if (end - start < 1e-8) next
      midpoint <- (start + end) / 2
      covering <- cuts$start < midpoint & cuts$end > midpoint
      hard <- any(cuts$hard[covering])
      covered <- any(covering & !cuts$hard)
      if (hard || (covered && identical(overlap, "clip"))) next
      faded <- covered && identical(overlap, "fade")
      alpha <- if (faded) overlap_alpha else 1
      if (alpha <= 0) next
      piece <- segments[s, , drop = FALSE]
      piece$x0 <- segments$x0[s] + start * dx
      piece$y0 <- segments$y0[s] + start * dy
      piece$x1 <- segments$x0[s] + end * dx
      piece$y1 <- segments$y0[s] + end * dy
      piece$occluded <- faded
      piece$alpha <- alpha
      out[[length(out) + 1L]] <- piece
    }
  }
  if (length(out) == 0) return(annotate(segments[0, , drop = FALSE]))
  result <- do.call(rbind, out)
  # Draw faint pieces first. Normal pieces then remain crisp at shared
  # endpoints, and the composite geom draws all text after both kinds.
  result[order(!result$occluded, seq_len(nrow(result))), , drop = FALSE]
}

#' Validate the gene leader-line linetype argument
#'
#' `gene_label_segment_linetype` accepts the special value `"auto"` (solid
#' lines, except dashed for labels moved to the other side of their arc) or
#' any valid ggplot2 linetype (character name or numeric dash pattern).
#' @keywords internal
validate_gene_segment_linetype <- function(lt) {
  if (is.null(lt)) return("auto")
  if (identical(lt, "auto")) return("auto")
  if (is.numeric(lt) && length(lt) >= 1 && all(is.finite(lt))) return(lt)
  ok <- c("blank", "solid", "dashed", "dotted",
          "dotdash", "longdash", "twodash")
  if (is.character(lt) && length(lt) >= 1 && all(lt %in% ok)) return(lt)
  ggchord_stop("gene_label_segment_linetype must be 'auto' or a valid ggplot2 ",
       "linetype (e.g. 'solid', 'dashed', 'dotted', or a numeric dash ",
       "pattern)", call. = FALSE)
}

#' Sample the plot content as repulsive points for label repulsion
#'
#' Collects a sparse set of points along the sequence arcs, gene arrows and
#' axes so that repelled gene labels avoid overlapping the plot content.
#' @keywords internal
ggchord_repel_points <- function(seq_arcs, gene_polys, axis_lines, axis_ticks,
                                 show_axis = FALSE) {
  pts <- list()
  if (length(seq_arcs) > 0) {
    for (arc in seq_arcs) {
      if (nrow(arc) == 0) next
      idx <- seq(1, nrow(arc), by = 8)
      pts[[length(pts) + 1]] <- arc[idx, c("x", "y"), drop = FALSE]
    }
  }
  if (nrow(gene_polys) > 0) {
    idx <- seq(1, nrow(gene_polys), by = 5)
    pts[[length(pts) + 1]] <- gene_polys[idx, c("x", "y"), drop = FALSE]
  }
  if (show_axis) {
    if (nrow(axis_lines) > 0) {
      idx <- seq(1, nrow(axis_lines), by = 8)
      pts[[length(pts) + 1]] <- axis_lines[idx, c("x", "y"), drop = FALSE]
    }
    if (nrow(axis_ticks) > 0) {
      pts[[length(pts) + 1]] <- data.frame(
        x = c(axis_ticks$x0, axis_ticks$x1, axis_ticks$label_x),
        y = c(axis_ticks$y0, axis_ticks$y1, axis_ticks$label_y),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(pts) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  do.call(rbind, pts)
}

#' Hide text labels that overlap the plot content or each other
#'
#' Estimates each label box (using the measured text size) and sets
#' \code{label} to NA when the box overlaps the given content points or another
#' label box. The first and last label of each sequence (axis start/end) are
#' always kept.
#' @keywords internal
ggchord_hide_text_overlaps <- function(df, content_pts,
                                       units_per_inch = 0.35) {
  idx <- which(!is.na(df$label))
  if (length(idx) < 1) return(df)
  # keep the axis start/end labels of every sequence visible
  protect <- logical(length(idx))
  if ("accver" %in% names(df) && length(idx) > 1) {
    for (sid in unique(df$accver[idx])) {
      rows <- idx[df$accver[idx] == sid]
      if (length(rows) > 1) {
        protect[match(min(rows), idx)] <- TRUE
        protect[match(max(rows), idx)] <- TRUE
      } else if (length(rows) == 1) {
        protect[match(rows, idx)] <- TRUE
      }
    }
  }
  if (nrow(content_pts) == 0) {
    content_pts <- data.frame(x = numeric(0), y = numeric(0))
  }
  close_device <- ggchord_measurement_device()
  on.exit(close_device())
  sizes <- df$size[idx] %||% rep(3, length(idx))
  w <- suppressWarnings(graphics::strwidth(df$label[idx], units = "inches",
                                           cex = sizes / 12)) * units_per_inch
  h <- suppressWarnings(graphics::strheight(df$label[idx], units = "inches",
                                            cex = sizes / 12)) * units_per_inch

  hide <- logical(length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]
    x <- df$label_x[i]
    y <- df$label_y[i]
    # overlap with content points: label center too close to an element
    if (nrow(content_pts) > 0) {
      d <- sqrt((x - content_pts$x)^2 + (y - content_pts$y)^2)
      if (any(d < min(w[k], h[k]) * 0.5 + 0.05)) hide[k] <- TRUE
    }
    if (!hide[k]) {
      # overlap with another label's box
      for (k2 in seq_along(idx)) {
        if (k2 == k) next
        j <- idx[k2]
        if (abs(x - df$label_x[j]) < (w[k] + w[k2]) / 2 &&
            abs(y - df$label_y[j]) < (h[k] + h[k2]) / 2) {
          hide[k] <- TRUE
          break
        }
      }
    }
  }
  hide <- hide & !protect
  df$label[idx[hide]] <- NA
  df
}
