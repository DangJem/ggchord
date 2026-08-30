#' Default categorical palette (Set1)
#'
#' The nine Set1 colors (previously obtained from the RColorBrewer package),
#' hardcoded so the package has no dependency on RColorBrewer. The colors are
#' identical to RColorBrewer's Set1 palette.
#' @keywords internal
chord_palette_set1 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                        "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

#' Generate a default categorical palette
#'
#' Returns the first \code{n} Set1 colors, interpolating with
#' \code{colorRampPalette()} when \code{n} exceeds 9.
#' @param n Number of colors requested
#' @return A character vector of \code{n} colors
#' @keywords internal
chord_default_palette <- function(n) {
  if (n <= 9) {
    chord_palette_set1[seq_len(n)]
  } else {
    colorRampPalette(chord_palette_set1)(n)
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
    bw = bw, bh = bh,
    xmin = x + cx_off - bw / 2, xmax = x + cx_off + bw / 2,
    ymin = y + cy_off - bh / 2, ymax = y + cy_off + bh / 2,
    stringsAsFactors = FALSE
  )
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
    return(data.frame(xmin = numeric(0), ymin = numeric(0),
                      xmax = numeric(0), ymax = numeric(0),
                      stringsAsFactors = FALSE))
  }
  boxes <- do.call(rbind, out)
  data.frame(xmin = boxes$xmin, ymin = boxes$ymin,
             xmax = boxes$xmax, ymax = boxes$ymax,
             stringsAsFactors = FALSE)
}

# Separate axis-aligned label boxes from one another and from fixed boxes.
#
# Collision detection is vectorised because the old nested R loops spent most
# of their time visiting pairs that did not overlap.  Displacements are still
# applied in label order, preserving the deterministic behaviour of the
# original solver while only looping over actual collisions.
ggchord_separate_boxes <- function(x, y, bw, bh, cx_off, cy_off,
                                   repel_boxes = NULL, max_iter = 500,
                                   x_lim = c(-Inf, Inf),
                                   y_lim = c(-Inf, Inf), tol = 1e-7) {
  n <- length(x)
  if (n == 0 || max_iter <= 0) {
    return(list(x = x, y = y, iterations = 0L))
  }

  if (is.null(repel_boxes) || nrow(repel_boxes) == 0) {
    ob_cx <- ob_cy <- ob_hw <- ob_hh <- numeric(0)
  } else {
    ob_cx <- (repel_boxes$xmin + repel_boxes$xmax) / 2
    ob_cy <- (repel_boxes$ymin + repel_boxes$ymax) / 2
    ob_hw <- (repel_boxes$xmax - repel_boxes$xmin) / 2
    ob_hh <- (repel_boxes$ymax - repel_boxes$ymin) / 2
  }

  for (iter in seq_len(max_iter)) {
    old_x <- x
    old_y <- y
    cx <- x + cx_off
    cy <- y + cy_off
    moved <- FALSE

    if (n > 1) {
      # Matrices use labels as rows/columns. Restricting to the upper triangle
      # gives each unordered pair once; ordering restores the former i/j loop.
      dx <- outer(cx, cx, function(a, b) b - a)
      dy <- outer(cy, cy, function(a, b) b - a)
      ox <- outer(bw, bw, "+") / 2 - abs(dx)
      oy <- outer(bh, bh, "+") / 2 - abs(dy)
      hits <- which(upper.tri(dx) & ox > 0 & oy > 0, arr.ind = TRUE)
      if (nrow(hits) > 1) {
        hits <- hits[order(hits[, 1], hits[, 2]), , drop = FALSE]
      }
      if (nrow(hits) > 0) {
        moved <- TRUE
        for (hit in seq_len(nrow(hits))) {
          i <- hits[hit, 1]
          j <- hits[hit, 2]
          if (ox[i, j] <= oy[i, j]) {
            sgn <- if (abs(dx[i, j]) < 1e-8) {
              if (stats::runif(1) < 0.5) -1 else 1
            } else sign(dx[i, j])
            x[i] <- x[i] - sgn * ox[i, j] / 2
            x[j] <- x[j] + sgn * ox[i, j] / 2
          } else {
            sgn <- if (abs(dy[i, j]) < 1e-8) {
              if (stats::runif(1) < 0.5) -1 else 1
            } else sign(dy[i, j])
            y[i] <- y[i] - sgn * oy[i, j] / 2
            y[j] <- y[j] + sgn * oy[i, j] / 2
          }
        }
      }
    }

    if (length(ob_cx) > 0) {
      dx_ob <- outer(cx, ob_cx, "-")
      dy_ob <- outer(cy, ob_cy, "-")
      ox_ob <- outer(bw / 2, ob_hw, "+") - abs(dx_ob)
      oy_ob <- outer(bh / 2, ob_hh, "+") - abs(dy_ob)
      hits <- which(ox_ob > 0 & oy_ob > 0, arr.ind = TRUE)
      if (nrow(hits) > 1) {
        hits <- hits[order(hits[, 1], hits[, 2]), , drop = FALSE]
      }
      if (nrow(hits) > 0) {
        moved <- TRUE
        for (hit in seq_len(nrow(hits))) {
          i <- hits[hit, 1]
          k <- hits[hit, 2]
          if (ox_ob[i, k] <= oy_ob[i, k]) {
            sgn <- if (abs(dx_ob[i, k]) < 1e-8) {
              if (stats::runif(1) < 0.5) -1 else 1
            } else sign(dx_ob[i, k])
            x[i] <- x[i] + sgn * ox_ob[i, k]
          } else {
            sgn <- if (abs(dy_ob[i, k]) < 1e-8) {
              if (stats::runif(1) < 0.5) -1 else 1
            } else sign(dy_ob[i, k])
            y[i] <- y[i] + sgn * oy_ob[i, k]
          }
        }
      }
    }

    x <- pmin(pmax(x, x_lim[1]), x_lim[2])
    y <- pmin(pmax(y, y_lim[1]), y_lim[2])
    displacement <- max(abs(x - old_x), abs(y - old_y))
    if (!moved || !is.finite(displacement) || displacement < tol) break
  }

  list(x = x, y = y, iterations = iter)
}

ggchord_sum_by_index <- function(index, value, n) {
  out <- numeric(n)
  if (length(index) == 0) return(out)
  totals <- rowsum(value, index, reorder = FALSE)
  out[as.integer(rownames(totals))] <- totals[, 1]
  out
}

ggchord_repel_labels <- function(gl, units_per_inch = 0.35,
                                 max_overlaps = Inf, box_padding = 0.25,
                                 point_padding = 0.1, min_segment_length = 0.5,
                                 force = 1, seed = 123, repel_points = NULL,
                                 repel_boxes = NULL) {
  n <- nrow(gl)
  empty_segments <- data.frame(x0 = numeric(0), y0 = numeric(0),
                               x1 = numeric(0), y1 = numeric(0),
                               group = integer(0), stringsAsFactors = FALSE)
  if (n == 0) return(list(labels = gl, segments = empty_segments))

  set.seed(seed)
  boxes <- ggchord_text_boxes(gl,
                              units_per_inch = units_per_inch,
                              box_padding = box_padding)
  w <- boxes$w
  h <- boxes$h
  bw <- boxes$bw
  bh <- boxes$bh
  cx_off <- boxes$cx - boxes$x
  cy_off <- boxes$cy - boxes$y

  # Anchors are the original (fixed) label positions next to the genes and are
  # used as the leader-line origins; labels start at their text positions.
  if ("anchor_x" %in% names(gl)) {
    ax <- gl$anchor_x
    ay <- gl$anchor_y
  } else {
    ax <- gl$text_x
    ay <- gl$text_y
  }

  # Start from a small, seed-dependent jitter around the anchor so that
  # different seeds lead to different (but reproducible) layouts.
  x <- gl$text_x + stats::runif(n, -0.5, 0.5) * 0.08 * max(w)
  y <- gl$text_y + stats::runif(n, -0.5, 0.5) * 0.08 * max(w)

  if (is.null(repel_points) || nrow(repel_points) == 0) {
    repel_points <- data.frame(x = numeric(0), y = numeric(0))
  }
  rpx <- repel_points$x
  rpy <- repel_points$y
  n_repel <- length(rpx)

  if (is.null(repel_boxes) || nrow(repel_boxes) == 0) {
    repel_boxes <- data.frame(xmin = numeric(0), ymin = numeric(0),
                              xmax = numeric(0), ymax = numeric(0),
                              stringsAsFactors = FALSE)
  }
  ob <- repel_boxes
  n_boxes <- nrow(ob)

  # Keep labels inside a generous region around the anchors.
  x_lim <- range(ax) + c(-1, 1) * (max(bw) + 1.5)
  y_lim <- range(ay) + c(-1, 1) * (max(bh) + 1.5)

  stable_iterations <- 0L
  for (iter in seq_len(300)) {
    old_x <- x
    old_y <- y
    fx <- numeric(n)
    fy <- numeric(n)

    # 1) short-range repulsion from every anchor point. The pairwise work is
    # evaluated in vectorised C code instead of two interpreted R loops.
    dx_anchor <- outer(x, ax, "-")
    dy_anchor <- outer(y, ay, "-")
    d_anchor <- sqrt(dx_anchor^2 + dy_anchor^2)
    cutoff <- point_padding + 0.4
    near_anchor <- d_anchor < cutoff
    d_anchor[d_anchor < 1e-4] <- 1e-4
    f_anchor <- force * 0.1 * (1 - d_anchor / cutoff)
    f_anchor[!near_anchor] <- 0
    fx <- fx + rowSums(f_anchor * dx_anchor / d_anchor)
    fy <- fy + rowSums(f_anchor * dy_anchor / d_anchor)

    # 1b) repulsion from sampled plot content (arcs, genes, axes).
    if (n_repel > 0) {
      cutoff_content <- point_padding + 0.3
      for (i in seq_len(n)) {
        dx <- x[i] - rpx
        dy <- y[i] - rpy
        d <- sqrt(dx^2 + dy^2)
        keep <- d < cutoff_content
        if (any(keep)) {
          dd <- d[keep]
          dd[dd < 1e-4] <- 1e-4
          f <- force * 0.08 * (1 - dd / cutoff_content)
          fx[i] <- fx[i] + sum(f * dx[keep] / dd)
          fy[i] <- fy[i] + sum(f * dy[keep] / dd)
        }
      }
    }

    # 1c) repulsion from fixed text obstacle rectangles (sequence labels,
    # group labels and axis labels).  Using the rectangles directly means a
    # long label cannot overlap an axis label even when its centre is far from
    # the axis-label anchor.
    if (n_boxes > 0) {
      cutoff_box <- point_padding + 0.25
      for (i in seq_len(n)) {
        cx <- pmin(pmax(x[i], ob$xmin), ob$xmax)
        cy <- pmin(pmax(y[i], ob$ymin), ob$ymax)
        dx <- x[i] - cx
        dy <- y[i] - cy
        d <- sqrt(dx^2 + dy^2)

        outside <- d >= 1e-4 & d < cutoff_box
        if (any(outside)) {
          dd <- d[outside]
          f <- force * 0.10 * (1 - dd / cutoff_box)
          fx[i] <- fx[i] + sum(f * dx[outside] / dd)
          fy[i] <- fy[i] + sum(f * dy[outside] / dd)
        }

        inside <- d < 1e-4
        if (any(inside)) {
          for (k in which(inside)) {
            left <- x[i] - ob$xmin[k]
            right <- ob$xmax[k] - x[i]
            bottom <- y[i] - ob$ymin[k]
            top <- ob$ymax[k] - y[i]
            m <- min(left, right, bottom, top)
            if (identical(m, left)) {
              fx[i] <- fx[i] - force * 0.12
            } else if (identical(m, right)) {
              fx[i] <- fx[i] + force * 0.12
            } else if (identical(m, bottom)) {
              fy[i] <- fy[i] - force * 0.12
            } else {
              fy[i] <- fy[i] + force * 0.12
            }
          }
        }
      }
    }

    # 2) repulsion between actual label boxes. Detect all colliding pairs at
    # once, then accumulate their forces by label index.
    if (n > 1) {
      cx <- x + cx_off
      cy <- y + cy_off
      dx_pair <- outer(cx, cx, function(a, b) b - a)
      dy_pair <- outer(cy, cy, function(a, b) b - a)
      ox_pair <- outer(bw, bw, "+") / 2 - abs(dx_pair)
      oy_pair <- outer(bh, bh, "+") / 2 - abs(dy_pair)
      hits <- which(upper.tri(dx_pair) & ox_pair > 0 & oy_pair > 0,
                    arr.ind = TRUE)
      if (nrow(hits) > 0) {
        ii <- hits[, 1]
        jj <- hits[, 2]
        horizontal <- ox_pair[hits] <= oy_pair[hits]
        if (any(horizontal)) {
          hi <- ii[horizontal]
          hj <- jj[horizontal]
          hd <- dx_pair[hits[horizontal, , drop = FALSE]]
          hs <- sign(hd)
          tied <- abs(hd) < 1e-8
          hs[tied] <- ifelse(stats::runif(sum(tied)) < 0.5, -1, 1)
          hf <- force * ox_pair[hits[horizontal, , drop = FALSE]] /
            pmax(bw[hi] + bw[hj], 1e-8)
          fx <- fx + ggchord_sum_by_index(hi, -hs * hf, n) +
            ggchord_sum_by_index(hj, hs * hf, n)
        }
        if (any(!horizontal)) {
          vi <- ii[!horizontal]
          vj <- jj[!horizontal]
          vd <- dy_pair[hits[!horizontal, , drop = FALSE]]
          vs <- sign(vd)
          tied <- abs(vd) < 1e-8
          vs[tied] <- ifelse(stats::runif(sum(tied)) < 0.5, -1, 1)
          vf <- force * oy_pair[hits[!horizontal, , drop = FALSE]] /
            pmax(bh[vi] + bh[vj], 1e-8)
          fy <- fy + ggchord_sum_by_index(vi, -vs * vf, n) +
            ggchord_sum_by_index(vj, vs * vf, n)
        }
      }
    }

    # 3) spring back toward the label's own starting position.
    x <- x + (gl$text_x - x) * 0.10 + fx * 0.65
    y <- y + (gl$text_y - y) * 0.10 + fy * 0.65
    x <- pmin(pmax(x, x_lim[1]), x_lim[2])
    y <- pmin(pmax(y, y_lim[1]), y_lim[2])
    displacement <- max(abs(x - old_x), abs(y - old_y))
    if (is.finite(displacement) && displacement < 1e-4) {
      stable_iterations <- stable_iterations + 1L
      if (stable_iterations >= 5L) break
    } else {
      stable_iterations <- 0L
    }
  }

  # Deterministic box-separation pass. This also treats text obstacles as hard
  # rectangles rather than as a cloud of points.
  separated <- ggchord_separate_boxes(
    x, y, bw, bh, cx_off, cy_off, repel_boxes = ob,
    max_iter = 500, x_lim = x_lim, y_lim = y_lim
  )
  x <- separated$x
  y <- separated$y

  # Leader lines: from anchor to the final label position.
  seg_dist <- sqrt((x - ax)^2 + (y - ay)^2)
  keep_seg <- seg_dist > min_segment_length
  segments <- data.frame(
    x0 = ax[keep_seg], y0 = ay[keep_seg],
    x1 = x[keep_seg], y1 = y[keep_seg],
    group = which(keep_seg),
    stringsAsFactors = FALSE
  )

  # Optional decluttering: hide labels that still overlap too many others.
  if (is.finite(max_overlaps)) {
    n_over <- numeric(n)
    cx <- x + cx_off
    cy <- y + cy_off
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        if (abs(cx[i] - cx[j]) < (bw[i] + bw[j]) / 2 &&
            abs(cy[i] - cy[j]) < (bh[i] + bh[j]) / 2) {
          n_over[i] <- n_over[i] + 1
          n_over[j] <- n_over[j] + 1
        }
      }
    }
    hide <- n_over > max_overlaps
    if (any(hide)) {
      gl$text[hide] <- NA
      keep_seg[hide] <- FALSE
      segments <- segments[segments$group %in% which(!hide), , drop = FALSE]
    }
  }

  gl$text_x <- x
  gl$text_y <- y
  list(labels = gl, segments = segments)
}

#' Final deterministic de-overlap pass for repelled gene labels.
#'
#' Runs after the horizontal/justification and arc-side adjustments so the
#' solver uses the exact rendered text boxes.  It also treats sequence, group
#' and axis labels as hard rectangular obstacles.
#' @keywords internal
ggchord_repel_labels_final <- function(gl, units_per_inch = 0.35,
                                       box_padding = 0.25,
                                       repel_boxes = NULL,
                                       max_iter = 500) {
  active <- !is.na(gl$text) & nzchar(gl$text)
  if (!any(active)) return(gl)
  active_rows <- which(active)
  work <- gl[active_rows, , drop = FALSE]

  boxes <- ggchord_text_boxes(work,
                              units_per_inch = units_per_inch,
                              box_padding = box_padding)
  bw <- boxes$bw
  bh <- boxes$bh
  cx_off <- boxes$cx - boxes$x
  cy_off <- boxes$cy - boxes$y

  x <- work$text_x
  y <- work$text_y
  if ("anchor_x" %in% names(work)) {
    ax <- work$anchor_x
    ay <- work$anchor_y
  } else {
    ax <- work$text_x
    ay <- work$text_y
  }
  x_lim <- range(c(ax, x)) + c(-1, 1) * (max(bw) + 1.5)
  y_lim <- range(c(ay, y)) + c(-1, 1) * (max(bh) + 1.5)

  if (is.null(repel_boxes) || nrow(repel_boxes) == 0) {
    repel_boxes <- data.frame(xmin = numeric(0), ymin = numeric(0),
                              xmax = numeric(0), ymax = numeric(0),
                              stringsAsFactors = FALSE)
  }
  ob <- repel_boxes

  separated <- ggchord_separate_boxes(
    x, y, bw, bh, cx_off, cy_off, repel_boxes = ob,
    max_iter = max_iter, x_lim = x_lim, y_lim = y_lim
  )
  x <- separated$x
  y <- separated$y

  gl$text_x[active_rows] <- x
  gl$text_y[active_rows] <- y
  gl
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

ggchord_elbow_bends <- function(gl, text_widths) {
  n <- nrow(gl)
  if (n == 0) return(data.frame(x = numeric(0), y = numeric(0)))
  horiz <- abs(gl$text_x - gl$anchor_x)
  stub_len <- pmin(pmax(0.02, 0.3 * horiz),
                   pmax(0.3 * text_widths, 0.04))
  dir <- ifelse(gl$hjust < 0.5, 1, -1)
  bx <- gl$text_x - dir * stub_len
  bx <- ifelse(gl$hjust < 0.5,
               pmax(bx, gl$anchor_x), pmin(bx, gl$anchor_x))
  data.frame(x = bx, y = gl$text_y)
}

ggchord_label_curve_frame <- function(gl, seq_arcs) {
  n <- nrow(gl)
  frame <- data.frame(
    curve_x = numeric(n), curve_y = numeric(n),
    outward_x = numeric(n), outward_y = numeric(n),
    signed_distance = numeric(n)
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

# Put horizontal labels into compact local bands around their own sequence.
# The free force layout supplies seed-dependent starting positions, but this
# pass constrains both radial and tangential drift around each gene anchor.
# Labels therefore use nearby two-dimensional space instead of escaping into
# one large global ring when sequences have very different radii.
ggchord_compact_label_lanes <- function(gl, seq_arcs,
                                        side = "outside",
                                        units_per_inch = 0.35,
                                        box_padding = 0.25,
                                        point_padding = 0.1,
                                        repel_boxes = NULL,
                                        max_iter = 100) {
  n <- nrow(gl)
  if (n == 0) return(list(labels = gl, lanes = character(0)))

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
  lanes <- paste(gl$seq_id, side_sign, sep = "\r")

  tangent_x <- -anchor_frame$outward_y
  tangent_y <- anchor_frame$outward_x
  direction_x <- side_sign * anchor_frame$outward_x
  gl$hjust[active] <- ifelse(direction_x[active] >= 0, 0, 1)
  gl$vjust[active] <- 0.5
  boxes <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  centre_normal <-
    (boxes$cx - boxes$x) * anchor_frame$outward_x +
    (boxes$cy - boxes$y) * anchor_frame$outward_y
  normal_extent <-
    (boxes$bw * abs(anchor_frame$outward_x) +
       boxes$bh * abs(anchor_frame$outward_y)) / 2
  clearance <- pmax(abs(anchor_frame$signed_distance), point_padding) + 0.04
  minimum_distance <- pmax(
    abs(anchor_frame$signed_distance) + 0.04,
    clearance + normal_extent - side_sign * centre_normal
  )
  band_depth <- max(0.22, 0.28 * units_per_inch)
  tangent_limit <- max(0.55, 0.75 * units_per_inch)

  vx <- gl$text_x - anchor_frame$curve_x
  vy <- gl$text_y - anchor_frame$curve_y
  tangent_offset <- vx * tangent_x + vy * tangent_y
  radial_distance <- side_sign *
    (vx * anchor_frame$outward_x + vy * anchor_frame$outward_y)
  tangent_offset <- pmin(pmax(tangent_offset, -tangent_limit), tangent_limit)
  radial_distance <- pmin(
    pmax(radial_distance, minimum_distance), minimum_distance + band_depth
  )
  gl$text_x[active] <- anchor_frame$curve_x[active] +
    tangent_offset[active] * tangent_x[active] +
    side_sign[active] * radial_distance[active] *
      anchor_frame$outward_x[active]
  gl$text_y[active] <- anchor_frame$curve_y[active] +
    tangent_offset[active] * tangent_y[active] +
    side_sign[active] * radial_distance[active] *
      anchor_frame$outward_y[active]

  # The final approach direction can differ from the radial direction when a
  # label uses tangential space. Justify from the actual gene-to-label leader
  # so elbow stubs always enter the empty side of the text.
  gl$hjust[active] <- ifelse(
    gl$text_x[active] >= gl$anchor_x[active], 0, 1
  )

  update_minimum_distance <- function(labels) {
    current_boxes <- ggchord_text_boxes(
      labels, units_per_inch = units_per_inch,
      box_padding = box_padding
    )
    current_centre_normal <-
      (current_boxes$cx - current_boxes$x) * anchor_frame$outward_x +
      (current_boxes$cy - current_boxes$y) * anchor_frame$outward_y
    current_normal_extent <-
      (current_boxes$bw * abs(anchor_frame$outward_x) +
         current_boxes$bh * abs(anchor_frame$outward_y)) / 2
    pmax(
      abs(anchor_frame$signed_distance) + 0.04,
      clearance + current_normal_extent -
        side_sign * current_centre_normal
    )
  }
  minimum_distance <- update_minimum_distance(gl)

  # Alternate one collision correction with a projection back into the local
  # band. If a lane is unusually dense, grow the band gradually instead of
  # sending one label far away in a single step.
  active_rows <- which(active)
  work_boxes <- ggchord_text_boxes(
    gl[active_rows, , drop = FALSE],
    units_per_inch = units_per_inch,
    box_padding = box_padding
  )
  for (iter in seq_len(max_iter)) {
    old_x <- gl$text_x[active_rows]
    old_y <- gl$text_y[active_rows]
    cx_off <- work_boxes$cx - work_boxes$x
    cy_off <- work_boxes$cy - work_boxes$y
    separated <- ggchord_separate_boxes(
      old_x, old_y, work_boxes$bw, work_boxes$bh,
      cx_off, cy_off, repel_boxes = repel_boxes,
      max_iter = 1
    )
    gl$text_x[active_rows] <- separated$x
    gl$text_y[active_rows] <- separated$y
    new_hjust <- ifelse(
      gl$text_x[active] >= gl$anchor_x[active], 0, 1
    )
    if (any(new_hjust != gl$hjust[active])) {
      gl$hjust[active] <- new_hjust
      minimum_distance <- update_minimum_distance(gl)
      work_boxes <- ggchord_text_boxes(
        gl[active_rows, , drop = FALSE],
        units_per_inch = units_per_inch,
        box_padding = box_padding
      )
    }

    extra <- 0.12 * floor((iter - 1L) / 40L)
    vx <- gl$text_x - anchor_frame$curve_x
    vy <- gl$text_y - anchor_frame$curve_y
    tangent_offset <- vx * tangent_x + vy * tangent_y
    radial_distance <- side_sign *
      (vx * anchor_frame$outward_x + vy * anchor_frame$outward_y)
    tangent_offset <- pmin(
      pmax(tangent_offset, -(tangent_limit + extra)),
      tangent_limit + extra
    )
    radial_distance <- pmin(
      pmax(radial_distance, minimum_distance),
      minimum_distance + band_depth + extra
    )
    gl$text_x[active] <- anchor_frame$curve_x[active] +
      tangent_offset[active] * tangent_x[active] +
      side_sign[active] * radial_distance[active] *
        anchor_frame$outward_x[active]
    gl$text_y[active] <- anchor_frame$curve_y[active] +
      tangent_offset[active] * tangent_y[active] +
      side_sign[active] * radial_distance[active] *
        anchor_frame$outward_y[active]
    new_hjust <- ifelse(
      gl$text_x[active] >= gl$anchor_x[active], 0, 1
    )
    if (any(new_hjust != gl$hjust[active])) {
      gl$hjust[active] <- new_hjust
      minimum_distance <- update_minimum_distance(gl)
      work_boxes <- ggchord_text_boxes(
        gl[active_rows, , drop = FALSE],
        units_per_inch = units_per_inch,
        box_padding = box_padding
      )
    }
    displacement <- max(
      abs(gl$text_x[active_rows] - old_x),
      abs(gl$text_y[active_rows] - old_y)
    )
    if (is.finite(displacement) && displacement < 1e-6) break
  }

  list(labels = gl, lanes = lanes)
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
