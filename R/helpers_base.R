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
    w[valid] <- suppressWarnings(graphics::strwidth(
      texts[valid], units = "inches", cex = sizes[valid] / 12
    )) * units_per_inch
    n_lines <- vapply(strsplit(texts[valid], "\n"), length, integer(1))
    h[valid] <- suppressWarnings(graphics::strheight(
      texts[valid], units = "inches", cex = sizes[valid] / 12
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
  n <- nrow(gl)
  if (n == 0) return(gl)

  boxes <- ggchord_text_boxes(gl,
                              units_per_inch = units_per_inch,
                              box_padding = box_padding)
  bw <- boxes$bw
  bh <- boxes$bh
  cx_off <- boxes$cx - boxes$x
  cy_off <- boxes$cy - boxes$y

  x <- gl$text_x
  y <- gl$text_y
  if ("anchor_x" %in% names(gl)) {
    ax <- gl$anchor_x
    ay <- gl$anchor_y
  } else {
    ax <- gl$text_x
    ay <- gl$text_y
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

  gl$text_x <- x
  gl$text_y <- y
  gl
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
  keep_seg <- seg_dist > min_segment_length
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
