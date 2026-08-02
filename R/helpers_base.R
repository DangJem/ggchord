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
get_plot_extremes <- function(allRibbon=NULL, seqArcs=NULL, axisLines=NULL, axisTicks=NULL, gene_arrows=NULL, gene_polys = NULL, seq_labels = NULL, show_axis = FALSE) {
  # Initialize vectors to store x and y coordinates
  x_coords <- numeric(0)
  y_coords <- numeric(0)

  # 1. Process ribbons (allRibbon)
  if (!is.null(allRibbon) && nrow(allRibbon) > 0) {
    x_coords <- c(x_coords, allRibbon$x)
    y_coords <- c(y_coords, allRibbon$y)
  }
  # 2. Process sequence arcs (seqArcs, convert list to data frame)
  if (!is.null(seqArcs) && length(seqArcs) > 0) {
    seq_df <- do.call(rbind, seqArcs)
    if (nrow(seq_df) > 0) {
      x_coords <- c(x_coords, seq_df$x)
      y_coords <- c(y_coords, seq_df$y)
    }
  }
  # 3. Process gene labels (gene_arrows)
  if (!is.null(gene_arrows) && nrow(gene_arrows) > 0) {
    x_coords <- c(x_coords, gene_arrows$text_x)
    y_coords <- c(y_coords, gene_arrows$text_y)
  }
  # 4. Process axis lines (axisLines)
  if (show_axis && !is.null(axisLines) && nrow(axisLines) > 0) {
    x_coords <- c(x_coords, axisLines$x)
    y_coords <- c(y_coords, axisLines$y)
  }
  # 5. Process tick marks (axisTicks)
  if (show_axis && !is.null(axisTicks) && nrow(axisTicks) > 0) {
    # x coordinates: x0, x1, label_x
    x_coords <- c(x_coords, axisTicks$x0, axisTicks$x1, axisTicks$label_x)
    # y coordinates: y0, y1, label_y
    y_coords <- c(y_coords, axisTicks$y0, axisTicks$y1, axisTicks$label_y)
  }
  # 6. Process gene arrow polygons (gene_polys): contains x, y coordinates (polygon vertices)
  if (!is.null(gene_polys) && nrow(gene_polys) > 0) {
    x_coords <- c(x_coords, gene_polys$x)
    y_coords <- c(y_coords, gene_polys$y)
  }
  # 7. Process sequence labels (seq_labels): contains text_x/text_y
  if (!is.null(seq_labels) && nrow(seq_labels) > 0) {
    x_coords <- c(x_coords, seq_labels$text_x)
    y_coords <- c(y_coords, seq_labels$text_y)
  }

  # Filter missing values (NA)
  x_coords <- x_coords[!is.na(x_coords)]
  y_coords <- y_coords[!is.na(y_coords)]

  # Calculate extremes (return list containing min/max for x and y)
  list(
    x_min = min(x_coords), # Left extreme
    x_max = max(x_coords), # Right extreme
    y_min = min(y_coords), # Bottom extreme
    y_max = max(y_coords) # Top extreme
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
ggchord_repel_labels <- function(gl, units_per_inch = 0.35,
                                 max_overlaps = Inf, box_padding = 0.25,
                                 point_padding = 0.1, min_segment_length = 0.5,
                                 force = 1, seed = 123, repel_points = NULL) {
  n <- nrow(gl)
  empty_segments <- data.frame(x0 = numeric(0), y0 = numeric(0),
                               x1 = numeric(0), y1 = numeric(0),
                               group = integer(0), stringsAsFactors = FALSE)
  if (n == 0) return(list(labels = gl, segments = empty_segments))

  set.seed(seed)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  sizes <- gl$size %||% rep(2.5, n)
  w <- suppressWarnings(graphics::strwidth(gl$text, units = "inches",
                                           cex = sizes / 12)) * units_per_inch
  n_lines <- vapply(strsplit(gl$text, "\n"), length, integer(1))
  h <- suppressWarnings(graphics::strheight(gl$text, units = "inches",
                                            cex = sizes / 12)) *
    n_lines * units_per_inch

  # Anchors are the original (fixed) label positions; labels start there.
  ax <- gl$text_x
  ay <- gl$text_y
  x <- ax
  y <- ay

  # Optional obstacle points (sequence arcs, gene arrows, axes) that the labels
  # should avoid.
  if (is.null(repel_points) || nrow(repel_points) == 0) {
    repel_points <- data.frame(x = numeric(0), y = numeric(0))
  }
  rpx <- repel_points$x
  rpy <- repel_points$y
  n_repel <- length(rpx)

  # Keep labels inside a generous region around the anchors.
  x_lim <- range(ax) + c(-1, 1) * (max(w) + 1.5)
  y_lim <- range(ay) + c(-1, 1) * (max(w) + 1.5)

  for (iter in seq_len(300)) {
    fx <- numeric(n)
    fy <- numeric(n)
    # 1) short-range repulsion from every anchor point (labels sit just off
    #    the genes / arcs, but do not fly away)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        dx <- x[i] - ax[j]
        dy <- y[i] - ay[j]
        d <- sqrt(dx^2 + dy^2)
        cutoff <- point_padding + 0.6
        if (d < cutoff) {
          if (d < 1e-4) d <- 1e-4
          f <- force * 0.25 * (1 - d / cutoff)
          fx[i] <- fx[i] + f * dx / d
          fy[i] <- fy[i] + f * dy / d
        }
      }
    }
    # 1b) repulsion from the plot content (arcs, genes, axes) so labels do
    #     not cover other elements
    if (n_repel > 0) {
      cutoff_content <- point_padding + 0.8
      for (i in seq_len(n)) {
        dx <- x[i] - rpx
        dy <- y[i] - rpy
        d <- sqrt(dx^2 + dy^2)
        keep <- d < cutoff_content
        if (any(keep)) {
          dd <- d[keep]
          dd[dd < 1e-4] <- 1e-4
          f <- force * 0.15 * (1 - dd / cutoff_content)
          fx[i] <- fx[i] + sum(f * dx[keep] / dd)
          fy[i] <- fy[i] + sum(f * dy[keep] / dd)
        }
      }
    }
    # 2) repulsion between labels (so their boxes do not overlap)
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        dx <- x[j] - x[i]
        dy <- y[j] - y[i]
        d <- sqrt(dx^2 + dy^2)
        if (d < 1e-4) d <- 1e-4
        cutoff <- ((w[i] + w[j]) / 2) * (1 + box_padding)
        if (d < cutoff) {
          f <- force * ((cutoff - d) / cutoff)
          fx[i] <- fx[i] - f * dx / d
          fy[i] <- fy[i] - f * dy / d
          fx[j] <- fx[j] + f * dx / d
          fy[j] <- fy[j] + f * dy / d
        }
      }
    }
    # 3) spring back toward the own anchor (keeps the label associated)
    x <- x + (ax - x) * 0.12 + fx * 0.5
    y <- y + (ay - y) * 0.12 + fy * 0.5
    # 4) clamp inside the region
    x <- pmin(pmax(x, x_lim[1]), x_lim[2])
    y <- pmin(pmax(y, y_lim[1]), y_lim[2])
  }

  # Leader lines: from anchor to the final label position
  seg_dist <- sqrt((x - ax)^2 + (y - ay)^2)
  keep_seg <- seg_dist > min_segment_length
  segments <- data.frame(
    x0 = ax[keep_seg], y0 = ay[keep_seg],
    x1 = x[keep_seg], y1 = y[keep_seg],
    group = which(keep_seg),
    stringsAsFactors = FALSE
  )

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
      idx <- seq(1, nrow(arc), by = 25)
      pts[[length(pts) + 1]] <- arc[idx, c("x", "y"), drop = FALSE]
    }
  }
  if (nrow(gene_polys) > 0) {
    idx <- seq(1, nrow(gene_polys), by = 15)
    pts[[length(pts) + 1]] <- gene_polys[idx, c("x", "y"), drop = FALSE]
  }
  if (show_axis) {
    if (nrow(axis_lines) > 0) {
      idx <- seq(1, nrow(axis_lines), by = 25)
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
