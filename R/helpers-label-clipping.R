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
#' @noRd
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
#' @noRd
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
#' @noRd
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
