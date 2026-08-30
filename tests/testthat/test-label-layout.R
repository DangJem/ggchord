label_boxes_overlap <- function(layout) {
  gl <- layout$gene_labels
  gl <- gl[!is.na(gl$text) & nzchar(gl$text), , drop = FALSE]
  if (nrow(gl) < 2) return(FALSE)
  boxes <- ggchord:::ggchord_text_boxes(
    gl, units_per_inch = layout$gene_label_clip_units
  )
  for (i in seq_len(nrow(boxes) - 1L)) {
    if (any(ggchord:::ggchord_oriented_box_overlaps(
      boxes[i, , drop = FALSE],
      boxes[(i + 1L):nrow(boxes), , drop = FALSE]
    ))) return(TRUE)
  }
  FALSE
}

same_lane_leaders_cross <- function(layout) {
  seg <- layout$gene_label_segments
  if (nrow(seg) < 2) return(FALSE)
  gl <- layout$gene_labels
  frame <- ggchord:::ggchord_label_curve_frame(gl, layout$seq_arcs)
  lane <- paste(
    gl$seq_id,
    ifelse(frame$signed_distance < 0, "inside", "outside"),
    sep = "\r"
  )

  for (i in seq_len(nrow(seg) - 1L)) {
    for (j in (i + 1L):nrow(seg)) {
      gi <- seg$group[i]
      gj <- seg$group[j]
      if (gi == gj || lane[gi] != lane[gj]) next
      if (ggchord:::ggchord_segments_cross(
        seg$x0[i], seg$y0[i], seg$x1[i], seg$y1[i],
        seg$x0[j], seg$y0[j], seg$x1[j], seg$y1[j]
      )) return(TRUE)
    }
  }
  FALSE
}

any_leaders_cross <- function(layout) {
  seg <- layout$gene_label_segments
  if (nrow(seg) < 2) return(FALSE)
  for (i in seq_len(nrow(seg) - 1L)) {
    for (j in (i + 1L):nrow(seg)) {
      if (seg$group[i] == seg$group[j]) next
      if (ggchord:::ggchord_segments_cross(
        seg$x0[i], seg$y0[i], seg$x1[i], seg$y1[i],
        seg$x0[j], seg$y0[j], seg$x1[j], seg$y1[j]
      )) return(TRUE)
    }
  }
  FALSE
}

leader_hits_other_label <- function(layout) {
  seg <- layout$gene_label_segments
  gl <- layout$gene_labels
  if (nrow(seg) == 0 || nrow(gl) < 2) return(FALSE)
  clip_units <- layout$gene_label_clip_units
  if (is.null(clip_units) || !is.finite(clip_units)) {
    clip_units <- ggchord:::ggchord_device_units_per_inch(
      unlist(lapply(layout$seq_arcs, `[[`, "x"), use.names = FALSE),
      unlist(lapply(layout$seq_arcs, `[[`, "y"), use.names = FALSE)
    )
  }
  boxes <- ggchord:::ggchord_text_boxes(
    gl, units_per_inch = clip_units
  )
  interval <- function(origin, delta, lower, upper, tol = 1e-8) {
    if (abs(delta) < tol) {
      if (origin <= lower + tol || origin >= upper - tol) return(NULL)
      return(c(-Inf, Inf))
    }
    sort(c((lower - origin) / delta, (upper - origin) / delta))
  }
  for (s in seq_len(nrow(seg))) {
    dx <- seg$x1[s] - seg$x0[s]
    dy <- seg$y1[s] - seg$y0[s]
    for (i in setdiff(seq_len(nrow(gl)), seg$group[s])) {
      angle <- boxes$angle[i]
      cos_a <- cos(angle)
      sin_a <- sin(angle)
      rel_x <- seg$x0[s] - boxes$cx[i]
      rel_y <- seg$y0[s] - boxes$cy[i]
      local_x0 <- rel_x * cos_a + rel_y * sin_a
      local_y0 <- -rel_x * sin_a + rel_y * cos_a
      local_dx <- dx * cos_a + dy * sin_a
      local_dy <- -dx * sin_a + dy * cos_a
      tx <- interval(local_x0, local_dx,
                     -boxes$w[i] / 2, boxes$w[i] / 2)
      ty <- interval(local_y0, local_dy,
                     -boxes$h[i] / 2, boxes$h[i] / 2)
      if (is.null(tx) || is.null(ty)) next
      enter <- max(0, tx[1], ty[1])
      exit <- min(1, tx[2], ty[2])
      if (enter < exit - 1e-8) return(TRUE)
    }
  }
  FALSE
}

test_that("repelled labels respect the complete geom_seq geometry", {
  data(seq_data_example)
  data(gene_data_example)
  ids <- seq_data_example$seq_id
  named <- function(x) stats::setNames(x, ids)

  cases <- list(
    default = list(),
    shape = list(
      seq_radius = named(c(0.75, 1.25, 0.95, 1.5)),
      seq_curvature = named(c(0, 0.35, 1.4, -0.6)),
      seq_gap = named(c(0.01, 0.05, 0.025, 0.07))
    ),
    ordered_grouped = list(
      seq_order = rev(ids),
      seq_radius = named(c(1.4, 0.85, 1.2, 0.7)),
      seq_curvature = named(c(1.6, 0, 0.55, 1.1)),
      seq_gap = named(c(0.04, 0.015, 0.06, 0.02)),
      seq_group = named(c("A", "A", "B", "B")),
      seq_group_gap = 0.06,
      seq_group_label_radius = 1.55
    ),
    presentation = list(
      seq_labels = named(paste("Sequence", seq_along(ids))),
      seq_colors = named(c("#3366AA", "#AA6633", "#339966", "#993399")),
      seq_group = named(c("A", "A", "B", "B")),
      seq_group_labels = c(A = "Group A", B = "Group B"),
      seq_group_colors = c(A = "navy", B = "firebrick"),
      linewidth = 1.7,
      show_legend = FALSE,
      legend_position = "left"
    )
  )

  for (case_index in seq_along(cases)) {
    case_name <- names(cases)[case_index]
    side <- "outside"
    for (orientation in c(1, -1)) {
      for (mode in c("aligned", "radial", "arc")) {
      seq_layer <- do.call(
        geom_seq,
        c(cases[[case_index]], list(seq_orientation = orientation))
      )
      p <- ggchord(
        seq_data_example, gene_data = gene_data_example,
        rotation = if (orientation == 1) 23 else 137
      ) + seq_layer + geom_gene() +
        geom_gene_label_repel(
          gene_label_layout = mode,
          gene_label_side = side
        )
      invisible(suppressWarnings(ggplot2::ggplot_build(p)))
      layout <- get_chord_layout()
      context <- paste(case_name, "orientation", orientation,
                       "mode", mode, "side", side)

      expect_false(label_boxes_overlap(layout), info = context)
      expect_false(same_lane_leaders_cross(layout), info = context)
      expect_false(any_leaders_cross(layout), info = context)
      expect_false(leader_hits_other_label(layout), info = context)
      frame <- ggchord:::ggchord_label_curve_frame(
        layout$gene_labels, layout$seq_arcs
      )
      expect_true(all(is.finite(frame$signed_distance)), info = context)
      if (side == "inside") {
        expect_true(all(frame$signed_distance <= 1e-7), info = context)
      } else {
        expect_true(all(frame$signed_distance >= -1e-7), info = context)
      }
      }
    }
  }
})

test_that("hidden repelled labels do not retain leader lines", {
  labels <- data.frame(
    text = c("visible", NA_character_, ""),
    anchor_x = c(0, 0, 0), anchor_y = c(0, 0, 0),
    text_x = c(1, 2, 3), text_y = c(1, 2, 3),
    hjust = c(0, 0, 0)
  )
  seg <- ggchord:::ggchord_repel_segments(labels, min_segment_length = 0)
  expect_equal(seg$group, 1L)
})

test_that("rotated own-label clipping ends at the oriented text edge", {
  label <- data.frame(
    text = "a fairly long label", text_x = 1, text_y = 0,
    text_angle = 45, size = 2.5, hjust = 0.5, vjust = 0.5
  )
  segment <- data.frame(
    x0 = 0, y0 = 0, x1 = 1, y1 = 0, group = 1L,
    stringsAsFactors = FALSE
  )
  units <- 0.5
  padding <- 0.01
  clipped <- ggchord:::ggchord_clip_segments_to_labels(
    segment, label, units_per_inch = units, padding = padding,
    include_own = TRUE
  )
  boxes <- ggchord:::ggchord_text_boxes(label, units_per_inch = units)
  angle <- boxes$angle[1]
  rel_x <- clipped$x1 - boxes$cx[1]
  rel_y <- clipped$y1 - boxes$cy[1]
  local_x <- rel_x * cos(angle) + rel_y * sin(angle)
  local_y <- -rel_x * sin(angle) + rel_y * cos(angle)
  boundary_ratio <- max(
    abs(local_x) / (boxes$w[1] / 2 + padding * units),
    abs(local_y) / (boxes$h[1] / 2 + padding * units)
  )
  expect_equal(nrow(clipped), 1L)
  expect_equal(boundary_ratio, 1, tolerance = 1e-7)
})

test_that("nested sequence radii use compact per-sequence label bands", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(
    seq_data_example, ribbon_data_example, gene_data_example,
    title = "ggchord"
  ) +
    geom_seq(
      seq_radius = c(3.3, 2.5, 1.8, 1.25),
      seq_orientation = c(1, -1, 1, -1)
    ) +
    geom_ribbon(ribbon_alpha = 0.45) +
    geom_gene() +
    geom_gene_label_repel() +
    geom_seq_label() +
    geom_axis()
  invisible(suppressWarnings(ggplot2::ggplot_build(p)))
  layout <- get_chord_layout()
  gl <- layout$gene_labels
  leader_length <- sqrt(
    (gl$text_x - gl$anchor_x)^2 + (gl$text_y - gl$anchor_y)^2
  )
  median_by_sequence <- tapply(leader_length, gl$seq_id, median)

  expect_lt(max(leader_length), 2)
  expect_true(all(median_by_sequence < 1.1))
  expect_false(label_boxes_overlap(layout))
  expect_false(same_lane_leaders_cross(layout))
  expect_false(any_leaders_cross(layout))
  expect_false(leader_hits_other_label(layout))

  blue <- gl$seq_id == "MT118296.1"
  expect_equal(diff(range(gl$text_x[blue])), 0, tolerance = 1e-8)
  expect_true(all(gl$hjust[blue] == 1 & gl$vjust[blue] == 0.5))

  green <- gl$seq_id == "OQ646790.1"
  expect_true(all(gl$hjust[green] == 0.5 & gl$vjust[green] == 1))
  expect_lte(length(unique(round(gl$text_y[green], 8))), 3)
})

test_that("leader clipping follows the physical output size", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  build_at_size <- function(mode, width, height) {
    p <- ggchord(
      seq_data_example, ribbon_data_example, gene_data_example,
      title = "ggchord"
    ) +
      geom_seq(
        seq_radius = c(3.3, 2.5, 1.8, 1.25),
        seq_orientation = c(1, -1, 1, -1)
      ) +
      geom_ribbon(ribbon_alpha = 0.45) +
      geom_gene() +
      geom_gene_label_repel(gene_label_layout = mode) +
      geom_seq_label() +
      geom_axis()
    path <- tempfile(fileext = ".png")
    grDevices::png(path, width = width, height = height,
                   units = "in", res = 72)
    result <- tryCatch({
      invisible(suppressWarnings(ggplot2::ggplot_build(p)))
      layout <- get_chord_layout()
      retained <- with(
        layout$gene_label_segments,
        sum(sqrt((x1 - x0)^2 + (y1 - y0)^2))
      )
      c(units = layout$gene_label_clip_units, retained = retained,
        overlap = label_boxes_overlap(layout),
        cross = any_leaders_cross(layout),
        hit = leader_hits_other_label(layout))
    }, finally = {
      grDevices::dev.off()
      unlink(path)
    })
    result
  }

  for (mode in c("aligned", "radial", "arc")) {
    small <- build_at_size(mode, 6, 6)
    large <- build_at_size(mode, 15, 10)
    expect_gt(unname(small["units"]), unname(large["units"]))
    expect_true(all(is.finite(c(small, large))), info = mode)
    expect_gt(unname(small["retained"]), 0)
    expect_gt(unname(large["retained"]), 0)
    expect_equal(unname(small[c("overlap", "cross", "hit")]),
                 c(0, 0, 0), info = paste(mode, "small"))
    expect_equal(unname(large[c("overlap", "cross", "hit")]),
                 c(0, 0, 0), info = paste(mode, "large"))
  }
})
