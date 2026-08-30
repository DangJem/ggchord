label_boxes_overlap <- function(layout) {
  gl <- layout$gene_labels
  gl <- gl[!is.na(gl$text) & nzchar(gl$text), , drop = FALSE]
  if (nrow(gl) < 2) return(FALSE)
  limits <- ggchord:::ggchord_adaptive_limits(layout)
  boxes <- ggchord:::ggchord_text_boxes(
    gl,
    units_per_inch = max(diff(limits$xlim), diff(limits$ylim), 1) / 6
  )
  dx <- abs(outer(boxes$cx, boxes$cx, "-"))
  dy <- abs(outer(boxes$cy, boxes$cy, "-"))
  any(upper.tri(dx) &
        dx < outer(boxes$bw, boxes$bw, "+") / 2 - 1e-7 &
        dy < outer(boxes$bh, boxes$bh, "+") / 2 - 1e-7)
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
      seq_layer <- do.call(
        geom_seq,
        c(cases[[case_index]], list(seq_orientation = orientation))
      )
      p <- ggchord(
        seq_data_example, gene_data = gene_data_example,
        rotation = if (orientation == 1) 23 else 137
      ) + seq_layer + geom_gene() +
        geom_gene_label_repel(
          seed = 7,
          gene_label_side = side,
          gene_label_orientation = "horizontal",
          gene_label_segment = "elbow"
        )
      invisible(suppressWarnings(ggplot2::ggplot_build(p)))
      layout <- get_chord_layout()
      context <- paste(case_name, "orientation", orientation, "side", side)

      expect_false(label_boxes_overlap(layout), info = context)
      expect_false(same_lane_leaders_cross(layout), info = context)
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
