# Direct checks of the new public workflows, not a testthat/release run.
pkgload::load_all(quiet = TRUE)
local({
  old <- getwd(); setwd(tempdir()); on.exit(setwd(old), add = TRUE)
  sequences <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
  values <- data.frame(accver = c("A", "B"), start = 100, end = 300, value = c(-2, 8))
  points <- data.frame(accver = "A", pos = seq(100, 600, 100), value = c(1, 2, NA, 4, 3, 2))
  plot <- ggchord(sequences) + geom_seq() +
    geom_track_bar(aes(track_value = value), data = values,
      position = position_track("coverage")) +
    geom_track_line(aes(track_value = value), data = points,
      position = position_track("signal", order = 2)) +
    scale_track_value_continuous("coverage", limits = c(-5, 10), unit = "reads") +
    coord_collinear()
  exported <- export_ggchord_layout(plot, include = c("track", "axis"))
  stopifnot(nrow(exported$track_registry) == 4L,
    all(vapply(exported$track_registry$limits[exported$track_registry$track_id == "coverage"],
      function(x) all(x == c(-5, 10)), logical(1))),
    any(exported$track_diagnostics$reason == "missing_value"),
    any(exported$axis$component == "track_axis"))
  stopifnot(inherits(ggplot2::ggplot_build(plot), "ggplot_built"))
  # Line chunks retain the NA gap; no single group may span both sides.
  line <- exported$track[exported$track$track_id == "signal" &
    exported$track$.component == "track_line", ]
  stopifnot(all(vapply(split(line$source_row, line$group), function(x)
    !(any(x < 3) && any(x > 3)), logical(1))))
  for (coordinate in list(coord_circular(), coord_chord())) {
    ring <- ggchord(sequences[1L, ]) + geom_seq() +
      geom_track_tile(aes(track_fill = category),
        data = transform(values[1L, ], category = "site"),
        position = position_track("kind")) +
      scale_track_fill_manual("kind", values = c(site = "#334455")) +
      geom_track_bar(aes(track_value = value), data = values[1L, ],
        position = position_track("coverage", order = 2)) +
      geom_track_line(aes(track_value = value), data = points,
        position = position_track("signal", order = 3)) + coordinate
    built <- ggplot2::ggplot_build(ring)
    key <- ggchord:::ggchord_track_aes("fill", "kind")
    stopifnot(all(built$data[[2L]][[key]] == "#334455"))
  }
  windows <- data.frame(accver = c("A", "A", "B", "B"),
    start = c(100, 400, 100, 400), end = c(200, 600, 200, 600),
    locus_id = c("A1", "A2", "B1", "B2"))
  view <- segment_ggchord_data(sequences, windows = windows)
  focused <- focus_ggchord_track(points, sequences, windows)
  segmented <- ggchord(view$seq_data) + geom_seq() +
    geom_track_line(aes(track_value = value), data = focused,
      position = position_track("signal")) +
    geom_seq_scale(length = 50) +
    scale_seq_position_continuous(sequences = "A2", coordinates = "source") +
    coord_collinear(segments = view$segments)
  out <- export_ggchord_layout(segmented)
  stopifnot(nrow(out$track) > 0L, all(out$axis_selection$selected == (out$axis_selection$accver == "A2")))
  stopifnot(inherits(ggplot2::ggplot_build(segmented), "ggplot_built"))
  report <- diagnose_ggchord(plot, width = 8, height = 5)
  if (isTRUE(report$measurement$converged)) {
    physical_width <- with(report$track_registry, (normal_max - normal_min) /
      report$measurement$units_per_mm[["y"]])
    stopifnot(all(abs(physical_width - 4) < .15))
  }
  invalid <- tryCatch(ggplot2::ggplot_build(plot +
    scale_track_value_continuous("coverage", limits = c(1, 10))), error = identity)
  stopifnot(inherits(invalid, "error"), grepl("baseline", conditionMessage(invalid)))
  empty <- points[0, ]
  invisible(ggplot2::ggplot_build(ggchord(sequences) + geom_seq() +
    geom_track_line(aes(track_value = value), data = empty,
      position = position_track("empty")) + coord_collinear()))
})
cat("Scoped quantitative tracks, coordinate mappings, gaps and physical bandwidth passed.\n")
