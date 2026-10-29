# Explicit local validation; no testthat, package check or device matrix.
pkgload::load_all(quiet = TRUE)
local({
  sequences <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
  windows <- data.frame(accver = c("A", "A", "B"),
    start = c(100, 700, 100), end = c(300, 900, 300),
    locus_id = c("A_left", "A_right", "B_left"),
    row_id = c("A", "A", "B"), orientation = c(1, -1, 1))
  view <- segment_ggchord_data(sequences, windows = windows)
  plot <- ggchord(view$seq_data) + geom_seq() +
    scale_seq_position_continuous(sequences = "A_right", coordinates = "source",
      breaks = c(700, 800, 900)) +
    geom_seq_scale(length = 100) +
    coord_collinear(segments = view$segments, segment_gap = 80)
  output <- export_ggchord_layout(plot, include = c("axis", "seq"))
  axes <- output$axis[output$axis$component != "seq_scale", , drop = FALSE]
  bars <- output$axis[output$axis$component == "seq_scale", , drop = FALSE]
  stopifnot(nrow(axes) > 0L, all(axes$accver == "A_right"),
    setequal(axes$label[!is.na(axes$label)], c("700", "800", "900")),
    nrow(bars) > 0L, all(bars$bp_length == 100))
  baseline <- bars[bars$.component == "line", ]
  scale <- output$sequence_layout$units_per_bp[1L]
  stopifnot(abs(diff(range(baseline$x)) / scale - 100) < 1e-7)
  reference <- map_ggchord_position(plot, "A_right", 800, source = TRUE)
  inverse <- map_ggchord_position(plot, "A_right", reference$x,
    inverse = TRUE, source = TRUE)
  stopifnot(abs(inverse$position - 800) < 1e-7)
  # Suppressing axes keeps backbones and an explicit shared scale bar.
  hidden <- ggchord_filter_sequence_axes(get_chord_layout(plot), character())
  stopifnot(nrow(hidden$axis_lines) == 0L, nrow(hidden$axis_ticks) == 0L,
    length(hidden$seq_arcs) == 3L)
  invalid <- tryCatch(ggchord_filter_sequence_axes(hidden, "A"), error = identity)
  stopifnot(inherits(invalid, "error"), grepl("window IDs", conditionMessage(invalid)))
  circular <- ggchord(data.frame(accver = "A", length = 1000)) +
    geom_seq() + coord_circular()
  invalid <- tryCatch(ggchord_seq_scale_geometry(data.frame(.scale_row = 1L),
    geom_seq_scale()$ggchord_params, get_chord_layout(circular)), error = identity)
  stopifnot(inherits(invalid, "error"), grepl("local sequence axes", conditionMessage(invalid)))
})
cat("Selective source/window axes and shared linear bp reference passed.\n")
