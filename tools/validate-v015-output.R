# Explicit v0.15 output acceptance, run only when requested.
# Rscript tools/validate-v015-output.R [output-directory]
# Teaching inputs are synthetic coordinates and arbitrary numeric values.
pkgload::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else tempfile("ggchord-v015-output-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("png", quietly = TRUE),
  requireNamespace("svglite", quietly = TRUE), requireNamespace("qpdf", quietly = TRUE))

sequences <- data.frame(accver = c("A", "B"), length = c(1000, 800))
bins <- data.frame(accver = rep(c("A", "B"), each = 4L),
  start = rep(c(100, 200, 400, 600), 2L), end = rep(c(150, 250, 450, 650), 2L),
  value = c(-1, .5, 1.5, .3, .4, -1.5, .5, 1))
points <- data.frame(accver = rep(c("A", "B"), each = 6L),
  pos = rep(c(100, 180, 250, 400, 500, 650), 2L),
  value = c(.5, 1, NA, -.5, 0, 1, 1, .5, NA, -1, -.5, .5))
sites <- data.frame(accver = c("A", "B"), start = c(200, 450),
  end = c(300, 550), kind = c("site", "reference"))
features <- data.frame(accver = c("A", "B"), start = c(100, 100), end = c(400, 400),
  strand = c("+", "-"), type = "CDS", label = c("LONG GENE LABEL A", "LONG GENE LABEL B"))

add_tracks <- function(plot, bars, line, tiles) {
  plot + geom_track_bar(aes(track_value = value), data = bars,
    position = position_track("signal", width = unit(4, "mm"))) +
    geom_track_line(aes(track_value = value), data = line, colour = "#273F61",
      position = position_track("signal", width = unit(4, "mm"))) +
    scale_track_value_continuous("signal", limits = c(-2, 2), unit = "arbitrary units") +
    geom_track_tile(aes(track_fill = kind), data = tiles,
      position = position_track("sites", side = "inside", width = unit(2, "mm"))) +
    scale_track_fill_manual("sites", values = c(site = "#468C86", reference = "#D1AB59"))
}
linear <- ggchord(sequences) + geom_seq(seq_gap = .6) +
  geom_feature(data = features, feature_width = .12, position = position_plasmid(.15)) +
  geom_feature_label_repel(data = features, label_fit = "shrink", label_min_size = 1.5,
    position = position_plasmid(.15), size = 4)
linear <- add_tracks(linear, bins, points, sites) +
  scale_seq_position_continuous(sequences = "B") + geom_seq_scale(length = 100) + coord_collinear()
curved <- function(coordinate, curvature = 1) {
  sequence_layer <- if (isTRUE(coordinate$ggchord_circular)) geom_seq(seq_radius = 1.2) else
    geom_seq(seq_radius = 1.2, seq_curvature = curvature, seq_orientation = -1)
  plot <- ggchord(sequences[1L, ]) +
    sequence_layer +
    geom_feature(data = features[1L, ], feature_width = .12, position = position_plasmid(-.2)) +
    geom_feature_label_repel(data = features[1L, ], label_fit = "shrink",
      position = position_plasmid(-.2), size = 4)
  add_tracks(plot, bins[bins$accver == "A", ], points[points$accver == "A", ], sites[1L, ]) + coordinate
}
windows <- data.frame(accver = rep(c("A", "B"), each = 2L),
  start = rep(c(50, 400), 2L), end = rep(c(300, 700), 2L),
  locus_id = c("A1", "A2", "B1", "B2"), row_id = rep(c("A", "B"), each = 2L),
  orientation = c(1, -1, 1, -1))
view <- segment_ggchord_data(sequences, windows = windows)
segmented <- add_tracks(ggchord(view$seq_data) + geom_seq(seq_gap = .6),
  focus_ggchord_track(bins, sequences, windows, type = "interval"),
  focus_ggchord_track(points, sequences, windows),
  focus_ggchord_track(sites, sequences, windows, type = "interval")) +
  scale_seq_position_continuous(sequences = c("A2", "B2"), coordinates = "source") +
  geom_seq_scale(length = 100) + coord_collinear(segments = view$segments, segment_gap = 80)
plots <- list(linear = linear, circular = curved(coord_circular(gap = 10)),
  chord = curved(coord_chord(rotation = 30), .8), segmented = segmented)

previous <- grDevices::dev.cur()
manifest <- list(); warnings <- list()
for (name in names(plots)) {
  plot <- plots[[name]] + labs(title = paste("v0.15", name, "output acceptance"))
  exported <- export_ggchord_layout(plot, include = c("track", "labels", "axis"))
  stopifnot(nrow(exported$track) > 0L, nrow(exported$axis) > 0L,
    length(unique(exported$track_registry$track_id)) == 2L,
    all(vapply(exported$track_registry$limits[exported$track_registry$track_id == "signal"],
      function(x) all(x == c(-2, 2)), logical(1))))
  for (extension in c("png", "pdf", "svg")) {
    file <- file.path(out, paste0(name, ".", extension))
    captured <- character()
    device <- switch(extension, png = function(filename, width, height, ...) {
      grDevices::png(filename, width = width, height = height, units = "in", res = 120,
        type = if (isTRUE(capabilities("aqua"))) "quartz" else getOption("bitmapType"), ...)
    }, pdf = grDevices::pdf, svg = svglite::svglite)
    withCallingHandlers(ggplot2::ggsave(file, plot, width = 10, height = 7,
      device = device, dpi = 120), warning = function(w) {
        captured <<- c(captured, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    stopifnot(file.exists(file), file.info(file)$size > 0,
      identical(grDevices::dev.cur(), previous), is.null(ggchord_render_context()))
    if (extension == "png") stopifnot(identical(dim(png::readPNG(file))[1:2], c(840L, 1200L)))
    if (extension == "pdf") stopifnot(qpdf::pdf_length(file) == 1L)
    if (extension == "svg") stopifnot(any(grepl("<svg", readLines(file, warn = FALSE), fixed = TRUE)))
    manifest[[length(manifest) + 1L]] <- data.frame(plot = name, device = extension,
      file = basename(file), bytes = file.info(file)$size, warning_count = length(captured))
    warnings[[paste(name, extension)]] <- captured
  }
}
write.csv(do.call(rbind, manifest), file.path(out, "manifest.csv"), row.names = FALSE)
saveRDS(warnings, file.path(out, "device-warnings.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
if (any(lengths(warnings))) stop("Output produced device warnings; inspect device-warnings.rds")
cat("PASS: three coordinates, segmented reversals, shared multilayer tracks and PNG/PDF/SVG output.\n")
