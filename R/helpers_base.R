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


#' Process panel margin parameters
#'
#' Standardizes input margin parameters into a list containing t (top), r (right), b (bottom), l (left). Supports single-value or list input.
#'
#' @param arg_list Numeric (single value) or list (named/unnamed), margin parameters
#' @return List containing four elements: t, r, b, l (numeric, margin sizes)
#' @keywords internal
process_panel_margin <- function(arg_list) {
  # Initialize result list with default values of 0
  result <- list(t = 0, r = 0, b = 0, l = 0)

  # Check if input is a list
  if (!is.list(arg_list)) {
    # Handle single-value input
    if (is.numeric(arg_list) && length(arg_list) == 1) {
      value <- arg_list
      result <- list(t = value, r = value, b = value, l = value)
      return(result)
    } else {
      warning("Input is not a valid list or single numeric value; default values will be used")
      return(result)
    }
  }

  # Handle empty list
  if (length(arg_list) == 0) {
    return(result)
  }

  # Handle named list
  if (!is.null(names(arg_list)) && all(names(arg_list) != "")) {
    valid_names <- c("t", "r", "b", "l")

    # Iterate over each element of the input list
    for (name in names(arg_list)) {
      if (name %in% valid_names) {
        # Check if value is numeric
        if (is.numeric(arg_list[[name]]) && length(arg_list[[name]]) == 1) {
          result[[name]] <- arg_list[[name]]
        } else {
          warning(paste("Parameter", name, "is not a single numeric value; default 0 will be used"))
        }
      } else {
        warning(paste("Unknown parameter", name, "will be ignored"))
      }
    }
  }
  # Handle unnamed list
  else {
    param_order <- c("t", "r", "b", "l")
    num_args <- length(arg_list)

    # Handle single-element unnamed list
    if (num_args == 1 && is.numeric(arg_list[[1]])) {
      value <- arg_list[[1]]
      result <- list(t = value, r = value, b = value, l = value)
      return(result)
    }

    # Assign values in order
    for (i in 1:min(num_args, length(param_order))) {
      if (is.numeric(arg_list[[i]]) && length(arg_list[[i]]) == 1) {
        result[[param_order[i]]] <- arg_list[[i]]
      } else {
        warning(paste("Parameter at position", i, "is not a single numeric value; default 0 will be used"))
      }
    }
  }

  return(result)
}


#' Calculate plot extremes
#'
#' Extracts x/y coordinate extremes from all plot elements (sequence arcs, ribbons, gene arrows, etc.) for adjusting the plot range
#'
#' @param allRibbon data.frame, ribbon data (with x, y columns), default NULL
#' @param seqArcs List, sequence arc data (each element is a data frame with x, y, seq_id), default NULL
#' @param axisLines data.frame, axis line data (with x, y, seq_id columns), default NULL
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

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
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
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
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
                                          margin_inches = 1.25) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  x_span <- if (length(x) > 1) diff(range(x)) else 0
  y_span <- if (length(y) > 1) diff(range(y)) else 0
  geometry_span <- max(x_span, y_span, 1)

  device_inches <- tryCatch(
    grDevices::dev.size("in"),
    error = function(e) c(NA_real_, NA_real_)
  )
  usable_inches <- suppressWarnings(min(device_inches, na.rm = TRUE)) -
    margin_inches
  if (!is.finite(usable_inches) || usable_inches < 2) {
    usable_inches <- fallback_inches
  }
  geometry_span / usable_inches
}

#' Convert text layers into fixed obstacle rectangles for label repulsion.
#' @keywords internal
ggchord_text_obstacle_boxes <- function(seq_labels_df = NULL,
                                        group_labels = NULL,
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

  if (!is.null(group_labels) && nrow(group_labels) > 0) {
    out[[length(out) + 1]] <- ggchord_text_boxes(
      group_labels,
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
  o1 * o2 < -tol && o3 * o4 < -tol
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
  if (is.null(lanes)) lanes <- gl$seq_id %||% rep("all", n)
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
    if (nrow(a) == 0) "" else as.character(unique(a$seq_id)[1])
  }, character(1))

  for (sid in unique(gl$seq_id)) {
    rows <- which(gl$seq_id == sid)
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

      # Either local normal is valid geometrically. Choose the one pointing
      # away from the chord centre, which works for circular, straight and
      # Bézier sequence paths and is independent of sequence orientation.
      nx <- -ty / tangent_length
      ny <- tx / tangent_length
      if (nx * arc$x[k] + ny * arc$y[k] < 0) {
        nx <- -nx
        ny <- -ny
      }
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

# Put horizontal labels on compact cardinal rails around their own sequence.
# Each label is classified from the actual local normal of the sequence curve,
# so radius, curvature, gap, rotation and sequence orientation are all already
# represented in the rail choice. Left/right rails form vertical columns and
# top/bottom rails form horizontal rows. Within each rail the original gene
# order is retained, making the resulting leaders planar by construction.
ggchord_compact_label_lanes <- function(gl, seq_arcs,
                                        side = "outside",
                                        units_per_inch = 0.35,
                                        box_padding = 0.25,
                                        point_padding = 0.1,
                                        repel_boxes = NULL,
                                        max_iter = 100,
                                        allow_corner = TRUE) {
  n <- nrow(gl)
  if (n == 0) return(list(labels = gl, lanes = character(0)))
  input_gl <- gl

  active <- !is.na(gl$text) & nzchar(gl$text)
  anchor_labels <- gl
  anchor_labels$text_x <- anchor_labels$anchor_x
  anchor_labels$text_y <- anchor_labels$anchor_y
  anchor_frame <- ggchord_label_curve_frame(anchor_labels, seq_arcs)
  side_sign <- if (identical(side, "outside")) {
    rep(1, n)
  } else if (identical(side, "inside")) {
    rep(-1, n)
  } else {
    ifelse(anchor_frame$signed_distance < 0, -1, 1)
  }
  desired_x <- side_sign * anchor_frame$outward_x
  desired_y <- side_sign * anchor_frame$outward_y
  sequence_x_min <- sequence_x_max <- sequence_x_centre <- numeric(0)
  for (arc in seq_arcs) {
    if (is.null(arc) || nrow(arc) == 0 || !any(is.finite(arc$x))) next
    sid <- as.character(arc$seq_id[1])
    xr <- range(arc$x[is.finite(arc$x)])
    sequence_x_min[sid] <- xr[1]
    sequence_x_max[sid] <- xr[2]
    sequence_x_centre[sid] <- mean(xr)
  }
  directions <- rep(NA_character_, n)
  direction_groups <- paste(gl$seq_id, side_sign, sep = "\r")
  for (rows in split(which(active), direction_groups[active], drop = TRUE)) {
    mean_x <- mean(desired_x[rows])
    mean_y <- mean(desired_y[rows])
    directions[rows] <- if (abs(mean_x) >= abs(mean_y)) {
      if (mean_x < 0) "left" else "right"
    } else {
      if (mean_y < 0) "bottom" else "top"
    }
  }
  # A crowded top/bottom sequence can leave an adjacent outer corner unused.
  # Move only the outermost third of its labels onto that side column. For a
  # nearly centred sequence, use the side where its annotated genes cluster.
  # The column is still pushed along the real local outward normals below, so
  # this reduces top/bottom expansion without an origin/radius approximation.
  corner_lane <- rep(FALSE, n)
  corner_threshold <- max(0.18, 0.22 * units_per_inch)
  for (rows in split(which(active), direction_groups[active], drop = TRUE)) {
    if (!isTRUE(allow_corner)) break
    direction <- directions[rows[1]]
    if (!direction %in% c("top", "bottom") || length(rows) < 4L ||
        !all(side_sign[rows] > 0)) next
    sid <- as.character(gl$seq_id[rows[1]])
    centre <- unname(sequence_x_centre[sid])
    if (length(centre) != 1L || !is.finite(centre)) next
    side_bias <- if (abs(centre) >= corner_threshold) {
      centre
    } else {
      mean(gl$anchor_x[rows])
    }
    if (!is.finite(side_bias) || abs(side_bias) < corner_threshold) next
    corner_direction <- if (side_bias < 0) "left" else "right"
    ord <- if (corner_direction == "right") {
      rows[order(gl$anchor_x[rows], decreasing = TRUE)]
    } else {
      rows[order(gl$anchor_x[rows])]
    }
    relief_count <- min(length(ord) - 2L, ceiling(length(ord) / 3))
    relief <- utils::head(ord, relief_count)
    directions[relief] <- corner_direction
    corner_lane[relief] <- TRUE
  }
  directions[!active] <- NA_character_
  base_lanes <- paste(gl$seq_id, directions, sep = "\r")

  gl$text_angle[active] <- 0
  gl$hjust[active] <- c(left = 1, right = 0, top = 0.5, bottom = 0.5)[
    directions[active]
  ]
  gl$vjust[active] <- c(left = 0.5, right = 0.5, top = 0, bottom = 1)[
    directions[active]
  ]
  boxes <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  axis_gap <- max(0.02, 0.025 * units_per_inch)
  # Long horizontal text cannot form one compact row without either a large
  # empty margin or very long leaders. Split only top/bottom lanes into a
  # small number of staggered parallel rows. Their shared endpoint ordering
  # keeps the leaders planar; vertical lanes deliberately remain one column
  # so their text anchors line up exactly.
  rail_index <- rep(1L, n)
  base_lane_rows <- split(which(active), base_lanes[active], drop = TRUE)
  for (rows in base_lane_rows) {
    if (!directions[rows[1]] %in% c("top", "bottom") ||
        length(rows) < 3) next
    ord <- rows[order(gl$anchor_x[rows], rows)]
    occupied <- boxes$bw[ord]
    available <- max(
      diff(range(gl$anchor_x[ord])) + 1.2 * units_per_inch,
      3 * units_per_inch
    )
    rail_count <- min(length(ord), ceiling(
      (sum(occupied) + axis_gap * (length(ord) - 1L)) / available
    ))
    if (rail_count < 2) next
    for (k in seq_along(ord)) {
      rail_index[ord[k]] <- ((k - 1L) %% rail_count) + 1L
    }
  }
  lanes <- paste(base_lanes, rail_index, sep = "\r")
  rail_gap <- point_padding + box_padding * units_per_inch +
    max(0.04, 0.04 * units_per_inch)
  rail_step <- max(0.16, 0.16 * units_per_inch) +
    max(boxes$h[active], 0)

  lane_rows <- split(which(active), lanes[active], drop = TRUE)
  for (rows in lane_rows) {
    direction <- directions[rows[1]]
    # Retain a small amount of seed-dependent variation without allowing it
    # to alter the gene order or the rail chosen for a label.
    if (direction %in% c("left", "right")) {
      preferred <- gl$anchor_y[rows] + 0.12 *
        (gl$text_y[rows] - gl$anchor_y[rows])
      before <- boxes$y[rows] - boxes$ymin[rows]
      after <- boxes$ymax[rows] - boxes$y[rows]
      gl$text_y[rows] <- ggchord_pack_label_axis(
        preferred, gl$anchor_y[rows], before, after, gap = axis_gap
      )
      sid <- as.character(gl$seq_id[rows[1]])
      gl$text_x[rows] <- if (direction == "left") {
        edge <- if (any(corner_lane[rows])) sequence_x_min[sid] else NA_real_
        if (!is.finite(edge)) edge <- min(anchor_frame$curve_x[rows])
        edge - rail_gap
      } else {
        edge <- if (any(corner_lane[rows])) sequence_x_max[sid] else NA_real_
        if (!is.finite(edge)) edge <- max(anchor_frame$curve_x[rows])
        edge + rail_gap
      }
      if (any(corner_lane[rows]) && length(rows) > 1L) {
        ord <- rows[order(gl$anchor_x[rows], rows)]
        step <- max(0.12, 0.40 * units_per_inch)
        offsets <- seq(0, by = step, length.out = length(ord))
        if (direction == "right") {
          gl$text_x[ord] <- gl$text_x[ord] + offsets
        } else {
          gl$text_x[ord] <- gl$text_x[ord] - rev(offsets)
        }
      }
    } else {
      base_rows <- which(base_lanes == base_lanes[rows[1]] & active)
      preferred <- gl$anchor_x[rows] + 0.12 *
        (gl$text_x[rows] - gl$anchor_x[rows])
      before <- boxes$x[rows] - boxes$xmin[rows]
      after <- boxes$xmax[rows] - boxes$x[rows]
      gl$text_x[rows] <- ggchord_pack_label_axis(
        preferred, gl$anchor_x[rows], before, after,
        gap = axis_gap
      )
      gl$text_y[rows] <- if (direction == "bottom") {
        min(anchor_frame$curve_y[base_rows]) - rail_gap -
          (rail_index[rows[1]] - 1L) * rail_step
      } else {
        max(anchor_frame$curve_y[base_rows]) + rail_gap +
          (rail_index[rows[1]] - 1L) * rail_step
      }
    }
  }

  # Parallel rows are staggered, so each row has room for long horizontal
  # text. Alternate row-specific spacing constraints with a weak global
  # monotonic constraint. This keeps endpoints in gene order without forcing
  # different rows to reserve each other's full text widths.
  for (rows in base_lane_rows) {
    if (!directions[rows[1]] %in% c("top", "bottom")) next
    anchor_order <- rows[order(gl$anchor_x[rows], rows)]
    # Rebalance horizontal rails towards the middle of the available top or
    # bottom sector, rather than keeping them under an often one-sided gene
    # cluster. Limit the tangential translation in physical units so nearby
    # free space is used without creating leaders across the whole plot.
    anchor_centre <- mean(gl$anchor_x[rows])
    max_rebalance <- max(0.18, 0.28 * units_per_inch)
    rebalance <- max(-max_rebalance, min(max_rebalance, -anchor_centre))
    target_centre <- mean(gl$text_x[rows]) + rebalance
    for (pass in seq_len(max_iter)) {
      old <- gl$text_x[rows]
      for (index in sort(unique(rail_index[rows]))) {
        same_rail <- rows[rail_index[rows] == index]
        same_rail <- same_rail[order(gl$anchor_x[same_rail], same_rail)]
        if (length(same_rail) < 2) next
        for (k in seq_len(length(same_rail) - 1L)) {
          i <- same_rail[k]
          j <- same_rail[k + 1L]
          required <- (boxes$bw[i] + boxes$bw[j]) / 2 + axis_gap
          shortage <- required - (gl$text_x[j] - gl$text_x[i])
          if (shortage > 0) {
            gl$text_x[i] <- gl$text_x[i] - shortage / 2
            gl$text_x[j] <- gl$text_x[j] + shortage / 2
          }
        }
      }
      if (length(anchor_order) > 1) {
        for (k in seq_len(length(anchor_order) - 1L)) {
          i <- anchor_order[k]
          j <- anchor_order[k + 1L]
          shortage <- axis_gap / 4 - (gl$text_x[j] - gl$text_x[i])
          if (shortage > 0) {
            gl$text_x[i] <- gl$text_x[i] - shortage / 2
            gl$text_x[j] <- gl$text_x[j] + shortage / 2
          }
        }
      }
      if (max(abs(gl$text_x[rows] - old)) < 1e-7) break
    }
    gl$text_x[rows] <- gl$text_x[rows] +
      target_centre - mean(gl$text_x[rows])
  }

  # Neighbouring sequences can share the same top or bottom sector when
  # radii and curvature differ strongly. Preserve endpoint order across those
  # sequences too, alternating it with the full within-row spacing rule.
  horizontal <- which(active & directions %in% c("top", "bottom"))
  if (length(horizontal) > 1) {
    for (pass in seq_len(max_iter)) {
      old <- gl$text_x[horizontal]
      for (rows in lane_rows) {
        if (!directions[rows[1]] %in% c("top", "bottom") ||
            length(rows) < 2) next
        ord <- rows[order(gl$anchor_x[rows], rows)]
        for (k in seq_len(length(ord) - 1L)) {
          i <- ord[k]
          j <- ord[k + 1L]
          required <- (boxes$bw[i] + boxes$bw[j]) / 2 + axis_gap
          shortage <- required - (gl$text_x[j] - gl$text_x[i])
          if (shortage > 0) {
            gl$text_x[i] <- gl$text_x[i] - shortage / 2
            gl$text_x[j] <- gl$text_x[j] + shortage / 2
          }
        }
      }
      for (direction in c("top", "bottom")) {
        rows <- which(active & directions == direction)
        if (length(rows) < 2) next
        ord <- rows[order(gl$anchor_x[rows], rows)]
        for (k in seq_len(length(ord) - 1L)) {
          i <- ord[k]
          j <- ord[k + 1L]
          shortage <- axis_gap / 4 - (gl$text_x[j] - gl$text_x[i])
          if (shortage > 0) {
            gl$text_x[i] <- gl$text_x[i] - shortage / 2
            gl$text_x[j] <- gl$text_x[j] + shortage / 2
          }
        }
      }
      if (max(abs(gl$text_x[horizontal] - old)) < 1e-7) break
    }
  }

  # A curved arc can still make two adjacent diagonal approaches intersect
  # even when their endpoint x order is monotone. Exchange only their nearby
  # parallel row levels while retaining each label's x position. This changes
  # the approach order without assigning either label to a distant slot.
  lane_crossing_count <- function(rows) {
    total <- 0L
    if (length(rows) < 2L) return(total)
    for (ii in seq_len(length(rows) - 1L)) {
      for (jj in seq.int(ii + 1L, length(rows))) {
        i <- rows[ii]
        j <- rows[jj]
        total <- total + ggchord_segments_cross(
          gl$anchor_x[i], gl$anchor_y[i], gl$text_x[i], gl$text_y[i],
          gl$anchor_x[j], gl$anchor_y[j], gl$text_x[j], gl$text_y[j]
        )
      }
    }
    total
  }
  for (pass in seq_len(max_iter)) {
    changed <- FALSE
    for (rows in base_lane_rows) {
      if (!directions[rows[1]] %in% c("top", "bottom") ||
          length(rows) < 2) next
      for (ii in seq_len(length(rows) - 1L)) {
        for (jj in (ii + 1L):length(rows)) {
          i <- rows[ii]
          j <- rows[jj]
          if (rail_index[i] == rail_index[j]) next
          if (ggchord_segments_cross(
            gl$anchor_x[i], gl$anchor_y[i], gl$text_x[i], gl$text_y[i],
            gl$anchor_x[j], gl$anchor_y[j], gl$text_x[j], gl$text_y[j]
          )) {
            before_crossings <- lane_crossing_count(rows)
            old_y <- gl$text_y[c(i, j)]
            gl$text_y[c(i, j)] <- gl$text_y[c(j, i)]
            if (lane_crossing_count(rows) < before_crossings) {
              changed <- TRUE
              break
            }
            gl$text_y[c(i, j)] <- old_y
          }
        }
        if (changed) break
      }
      if (changed) break
    }
    if (!changed) break
  }

  move_lane_outward <- function(rows, amount) {
    if (!is.finite(amount) || amount <= 0) return()
    if (any(corner_lane[rows])) {
      nx <- mean(desired_x[rows])
      ny <- mean(desired_y[rows])
      norm <- sqrt(nx^2 + ny^2)
      if (is.finite(norm) && norm > 1e-10) {
        gl$text_x[rows] <<- gl$text_x[rows] + amount * nx / norm
        gl$text_y[rows] <<- gl$text_y[rows] + amount * ny / norm
        return()
      }
    }
    direction <- directions[rows[1]]
    if (direction == "left") gl$text_x[rows] <<- gl$text_x[rows] - amount
    if (direction == "right") gl$text_x[rows] <<- gl$text_x[rows] + amount
    if (direction == "bottom") gl$text_y[rows] <<- gl$text_y[rows] - amount
    if (direction == "top") gl$text_y[rows] <<- gl$text_y[rows] + amount
  }

  # Keep every rail on the requested side even after the labels have been
  # packed along it. Oblique local normals may require a small outward shift
  # when a long lane extends tangentially beyond its sequence endpoint.
  ensure_requested_side <- function() {
    for (rows in base_lane_rows) {
      direction <- directions[rows[1]]
      dx <- gl$text_x[rows] - anchor_frame$curve_x[rows]
      dy <- gl$text_y[rows] - anchor_frame$curve_y[rows]
      signed <- side_sign[rows] *
        (dx * anchor_frame$outward_x[rows] +
           dy * anchor_frame$outward_y[rows])
      projection <- if (any(corner_lane[rows])) {
        nx <- mean(desired_x[rows])
        ny <- mean(desired_y[rows])
        norm <- sqrt(nx^2 + ny^2)
        if (!is.finite(norm) || norm < 1e-10) {
          rep(1, length(rows))
        } else {
          (nx / norm) * desired_x[rows] +
            (ny / norm) * desired_y[rows]
        }
      } else {
        switch(
          direction,
          left = -side_sign[rows] * anchor_frame$outward_x[rows],
          right = side_sign[rows] * anchor_frame$outward_x[rows],
          bottom = -side_sign[rows] * anchor_frame$outward_y[rows],
          top = side_sign[rows] * anchor_frame$outward_y[rows]
        )
      }
      amount <- max((rail_gap - signed) / pmax(projection, 0.1), 0)
      move_lane_outward(rows, amount)
    }
  }
  ensure_requested_side()

  # Sequence and axis text are fixed obstacles. Resolve those conflicts by
  # moving the complete rail outwards, which preserves its alignment and the
  # order of every leader instead of pushing individual labels off the rail.
  if (!is.null(repel_boxes) && nrow(repel_boxes) > 0) {
    for (pass in seq_len(8)) {
      moved <- FALSE
      boxes <- ggchord_text_boxes(
        gl, units_per_inch = units_per_inch, box_padding = 0
      )
      for (rows in lane_rows) {
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
        base_rows <- which(base_lanes == base_lanes[rows[1]] & active)
        move_lane_outward(base_rows, amount)
        moved <- TRUE
      }
      if (!moved) break
    }
  }
  ensure_requested_side()

  # Resolve the remaining visible-box conflicts between different sequence
  # rails by moving the already outer rail farther out as one rigid unit.
  # This keeps vertical columns and staggered rows aligned.
  for (pass in seq_len(12)) {
    visible_boxes <- ggchord_text_boxes(
      gl, units_per_inch = units_per_inch, box_padding = 0
    )
    dx <- abs(outer(visible_boxes$cx, visible_boxes$cx, "-"))
    dy <- abs(outer(visible_boxes$cy, visible_boxes$cy, "-"))
    hits <- which(
      upper.tri(dx) &
        dx < outer(visible_boxes$bw, visible_boxes$bw, "+") / 2 - 1e-7 &
        dy < outer(visible_boxes$bh, visible_boxes$bh, "+") / 2 - 1e-7,
      arr.ind = TRUE
    )
    if (nrow(hits) == 0) break
    moved <- FALSE
    for (hit in seq_len(nrow(hits))) {
      i <- hits[hit, 1]
      j <- hits[hit, 2]
      if (base_lanes[i] == base_lanes[j]) next
      if (directions[i] != directions[j]) next
      direction <- directions[i]
      lane_i <- which(base_lanes == base_lanes[i] & active)
      lane_j <- which(base_lanes == base_lanes[j] & active)
      amount <- if (direction %in% c("top", "bottom")) {
        (visible_boxes$bh[i] + visible_boxes$bh[j]) / 2 -
          abs(visible_boxes$cy[i] - visible_boxes$cy[j]) + axis_gap
      } else {
        (visible_boxes$bw[i] + visible_boxes$bw[j]) / 2 -
          abs(visible_boxes$cx[i] - visible_boxes$cx[j]) + axis_gap
      }
      outer_i <- switch(
        direction,
        top = mean(gl$text_y[lane_i]) >= mean(gl$text_y[lane_j]),
        bottom = mean(gl$text_y[lane_i]) <= mean(gl$text_y[lane_j]),
        left = mean(gl$text_x[lane_i]) <= mean(gl$text_x[lane_j]),
        right = mean(gl$text_x[lane_i]) >= mean(gl$text_x[lane_j])
      )
      move_lane_outward(if (outer_i) lane_i else lane_j, amount)
      moved <- TRUE
      break
    }
    if (!moved) break
  }
  ensure_requested_side()

  # Corner relief columns can meet a neighbouring horizontal rail at the
  # corner even though each rail is internally collision-free. Resolve only
  # those mixed-direction contacts by moving the complete corner column along
  # its mean local outward normal; individual labels never change order.
  for (pass in seq_len(16)) {
    visible_boxes <- ggchord_text_boxes(
      gl, units_per_inch = units_per_inch, box_padding = 0
    )
    dx <- abs(outer(visible_boxes$cx, visible_boxes$cx, "-"))
    dy <- abs(outer(visible_boxes$cy, visible_boxes$cy, "-"))
    hits <- which(
      upper.tri(dx) &
        dx < outer(visible_boxes$bw, visible_boxes$bw, "+") / 2 - 1e-7 &
        dy < outer(visible_boxes$bh, visible_boxes$bh, "+") / 2 - 1e-7,
      arr.ind = TRUE
    )
    if (nrow(hits) == 0) break
    corner_hit <- which(
      corner_lane[hits[, 1]] | corner_lane[hits[, 2]]
    )[1]
    if (is.na(corner_hit)) break
    i <- hits[corner_hit, 1]
    j <- hits[corner_hit, 2]
    chosen <- if (corner_lane[i]) i else j
    rows <- which(base_lanes == base_lanes[chosen] & active)
    move_lane_outward(rows, max(axis_gap, 0.08 * units_per_inch))
  }
  ensure_requested_side()

  # Resolve residual contacts between perpendicular or otherwise differently
  # directed rails. Move one complete rail outward at a time; this preserves
  # all within-rail ordering and alignment.
  for (pass in seq_len(32)) {
    final_boxes <- ggchord_text_boxes(
      gl, units_per_inch = units_per_inch, box_padding = 0.03
    )
    dx <- abs(outer(final_boxes$cx, final_boxes$cx, "-"))
    dy <- abs(outer(final_boxes$cy, final_boxes$cy, "-"))
    hits <- which(
      upper.tri(dx) &
        dx < outer(final_boxes$bw, final_boxes$bw, "+") / 2 - 1e-7 &
        dy < outer(final_boxes$bh, final_boxes$bh, "+") / 2 - 1e-7,
      arr.ind = TRUE
    )
    if (nrow(hits) == 0) break
    resolvable <- which(
      base_lanes[hits[, 1]] != base_lanes[hits[, 2]]
    )[1]
    if (is.na(resolvable)) break
    i <- hits[resolvable, 1]
    j <- hits[resolvable, 2]
    lane_i <- which(base_lanes == base_lanes[i] & active)
    lane_j <- which(base_lanes == base_lanes[j] & active)
    outward_score <- function(rows) {
      switch(
        directions[rows[1]],
        top = mean(gl$text_y[rows]),
        bottom = -mean(gl$text_y[rows]),
        left = -mean(gl$text_x[rows]),
        right = mean(gl$text_x[rows])
      )
    }
    chosen <- if (corner_lane[i] != corner_lane[j]) {
      if (corner_lane[i]) lane_i else lane_j
    } else if (outward_score(lane_i) >= outward_score(lane_j)) {
      lane_i
    } else {
      lane_j
    }
    move_lane_outward(chosen, max(axis_gap, 0.08 * units_per_inch))
  }
  ensure_requested_side()

  # Obstacle and inter-sequence rail shifts above can change approach angles.
  # Re-run the monotone row-level exchange on the final rail coordinates.
  for (pass in seq_len(max_iter)) {
    changed <- FALSE
    for (rows in base_lane_rows) {
      if (!directions[rows[1]] %in% c("top", "bottom") ||
          length(rows) < 2L) next
      for (ii in seq_len(length(rows) - 1L)) {
        for (jj in seq.int(ii + 1L, length(rows))) {
          i <- rows[ii]
          j <- rows[jj]
          if (rail_index[i] == rail_index[j] ||
              !ggchord_segments_cross(
                gl$anchor_x[i], gl$anchor_y[i], gl$text_x[i], gl$text_y[i],
                gl$anchor_x[j], gl$anchor_y[j], gl$text_x[j], gl$text_y[j]
              )) next
          before_crossings <- lane_crossing_count(rows)
          old_y <- gl$text_y[c(i, j)]
          gl$text_y[c(i, j)] <- gl$text_y[c(j, i)]
          if (lane_crossing_count(rows) < before_crossings) {
            changed <- TRUE
            break
          }
          gl$text_y[c(i, j)] <- old_y
        }
        if (changed) break
      }
      if (changed) break
    }
    if (!changed) break
  }

  # The corner column and its neighbouring horizontal rail have different
  # endpoint axes, so disjoint text boxes alone do not guarantee planar
  # approaches. Alternate final box and crossing checks while moving the
  # complete corner column outward.
  for (pass in seq_len(12)) {
    crossed_corner <- NA_integer_
    rows <- which(active)
    final_boxes <- ggchord_text_boxes(
      gl, units_per_inch = units_per_inch, box_padding = 0.03
    )
    if (length(rows) > 1L) {
      for (ii in seq_len(length(rows) - 1L)) {
        i <- rows[ii]
        for (jj in seq.int(ii + 1L, length(rows))) {
          j <- rows[jj]
          if (!(corner_lane[i] || corner_lane[j])) next
          if (ggchord_oriented_box_overlaps(
            final_boxes[i, , drop = FALSE],
            final_boxes[j, , drop = FALSE]
          )) {
            crossed_corner <- if (corner_lane[i]) i else j
            break
          }
        }
        if (!is.na(crossed_corner)) break
      }
    }
    if (length(rows) > 1L) {
      for (ii in seq_len(length(rows) - 1L)) {
        if (!is.na(crossed_corner)) break
        i <- rows[ii]
        for (jj in seq.int(ii + 1L, length(rows))) {
          j <- rows[jj]
          if (!(corner_lane[i] || corner_lane[j])) next
          if (ggchord_segments_cross(
            gl$anchor_x[i], gl$anchor_y[i], gl$text_x[i], gl$text_y[i],
            gl$anchor_x[j], gl$anchor_y[j], gl$text_x[j], gl$text_y[j]
          )) {
            crossed_corner <- if (corner_lane[i]) i else j
            break
          }
        }
        if (!is.na(crossed_corner)) break
      }
    }
    if (is.na(crossed_corner)) break
    lane <- which(base_lanes == base_lanes[crossed_corner] & active)
    move_lane_outward(lane, max(axis_gap, 0.08 * units_per_inch))
  }
  ensure_requested_side()

  if (isTRUE(allow_corner) && any(corner_lane)) {
    unresolved_corner <- FALSE
    rows <- which(active)
    final_boxes <- ggchord_text_boxes(
      gl, units_per_inch = units_per_inch, box_padding = 0.03
    )
    if (length(rows) > 1L) {
      for (ii in seq_len(length(rows) - 1L)) {
        i <- rows[ii]
        for (jj in seq.int(ii + 1L, length(rows))) {
          j <- rows[jj]
          if (!(corner_lane[i] || corner_lane[j])) next
          if (ggchord_oriented_box_overlaps(
              final_boxes[i, , drop = FALSE],
              final_boxes[j, , drop = FALSE]
            ) || ggchord_segments_cross(
              gl$anchor_x[i], gl$anchor_y[i], gl$text_x[i], gl$text_y[i],
              gl$anchor_x[j], gl$anchor_y[j], gl$text_x[j], gl$text_y[j]
            )) {
            unresolved_corner <- TRUE
            break
          }
        }
        if (unresolved_corner) break
      }
    }
    if (unresolved_corner) {
      return(ggchord_compact_label_lanes(
        input_gl, seq_arcs,
        side = side,
        units_per_inch = units_per_inch,
        box_padding = box_padding,
        point_padding = point_padding,
        repel_boxes = repel_boxes,
        max_iter = max_iter,
        allow_corner = FALSE
      ))
    }
  }

  list(
    labels = gl,
    lanes = lanes,
    directions = directions,
    corner = corner_lane
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
                draw_segment = logical(0), corner = logical(0)))
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

  base_lanes <- paste(gl$seq_id, side_sign, sep = "\r")
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
    draw_segment = draw_segment,
    corner = rep(FALSE, n)
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
  active <- !is.na(gl$text) & nzchar(gl$text)
  rows <- which(active)
  if (length(rows) == 0) return(gl)
  boxes <- ggchord_text_boxes(gl[rows, , drop = FALSE],
                              units_per_inch = units_per_inch)
  counts <- integer(length(rows))
  if (length(rows) > 1) {
    for (i in seq_len(nrow(boxes))) {
      other <- setdiff(seq_len(nrow(boxes)), i)
      counts[i] <- counts[i] + sum(ggchord_oriented_box_overlaps(
        boxes[i, , drop = FALSE], boxes[other, , drop = FALSE]
      ))
    }
  }
  if (!is.null(repel_boxes) && nrow(repel_boxes) > 0) {
    for (i in seq_len(nrow(boxes))) {
      counts[i] <- counts[i] + sum(ggchord_oriented_box_overlaps(
        boxes[i, , drop = FALSE], repel_boxes
      ))
    }
  }
  gl$text[rows[counts > max_overlaps]] <- NA_character_
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

# Collapse only elbow stubs that cause same-lane crossings. The corresponding
# leader becomes straight; all conflict-free elbows retain their original
# bend and stub lengths.
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
        if (gi == gj || is.na(lanes[gi]) || is.na(lanes[gj]) ||
            lanes[gi] != lanes[gj]) next
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

# Remove the portions of leader lines hidden by other labels. Staggered rails
# deliberately reuse horizontal space; clipping at the real text rectangles
# keeps those compact layouts readable without moving an outer-row label far
# away merely because its approach passes behind an inner-row label.
ggchord_clip_segments_to_labels <- function(segments, gl,
                                            units_per_inch = 0.35,
                                            padding = 0.01,
                                            include_own = FALSE) {
  if (nrow(segments) == 0 || nrow(gl) == 0 ||
      (!isTRUE(include_own) && nrow(gl) < 2)) return(segments)
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
    pieces <- matrix(c(0, 1), ncol = 2)
    others <- which(visible)
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

      kept <- list()
      for (p in seq_len(nrow(pieces))) {
        start <- pieces[p, 1]
        end <- pieces[p, 2]
        if (cut_end <= start || cut_start >= end) {
          kept[[length(kept) + 1L]] <- c(start, end)
        } else {
          if (cut_start > start + 1e-10) {
            kept[[length(kept) + 1L]] <- c(start, min(cut_start, end))
          }
          if (cut_end < end - 1e-10) {
            kept[[length(kept) + 1L]] <- c(max(cut_end, start), end)
          }
        }
      }
      if (length(kept) == 0) {
        pieces <- matrix(numeric(0), ncol = 2)
        break
      }
      pieces <- do.call(rbind, kept)
    }

    if (nrow(pieces) == 0) next
    for (p in seq_len(nrow(pieces))) {
      if (pieces[p, 2] - pieces[p, 1] < 1e-8) next
      piece <- segments[s, , drop = FALSE]
      piece$x0 <- segments$x0[s] + pieces[p, 1] * dx
      piece$y0 <- segments$y0[s] + pieces[p, 1] * dy
      piece$x1 <- segments$x0[s] + pieces[p, 2] * dx
      piece$y1 <- segments$y0[s] + pieces[p, 2] * dy
      out[[length(out) + 1L]] <- piece
    }
  }
  if (length(out) == 0) return(segments[0, , drop = FALSE])
  do.call(rbind, out)
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
  if ("seq_id" %in% names(df) && length(idx) > 1) {
    for (sid in unique(df$seq_id[idx])) {
      rows <- idx[df$seq_id[idx] == sid]
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
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
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
