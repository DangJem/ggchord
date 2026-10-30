# Explicit PNG/PDF/SVG acceptance for the v0.14.0 linear coordinate workflows.
# Run only during requested output/release acceptance:
# Rscript tools/validate-collinear-output.R [output-directory]
# Existing direct scripts check the coordinate and data contracts in detail.
pkgload::load_all(quiet = TRUE)
source(file.path("examples", "collinear.R"), local = TRUE)
source(file.path("examples", "collinear-workflow.R"), local = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else tempfile("ggchord-collinear-output-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("svglite", quietly = TRUE),
  requireNamespace("png", quietly = TRUE))

# Both feature and ribbon layers must survive independently in a real output.
independent <- ggchord(workflow_sequences, validate = "none") +
  geom_link_ribbon(data = workflow_links[1L, ], fill = "#4477AA", alpha = .2) +
  geom_link_ribbon(data = workflow_links[3L, ], fill = "#CC6677", alpha = .2) +
  geom_seq(seq_orientation = c(A = 1, B = -1, C = 1), seq_gap = .4) +
  geom_feature(data = workflow_features[1L, ], shape = "block",
    feature_fill = "#4477AA", position = position_strand(.12)) +
  geom_feature(data = workflow_features[2L, ], shape = "chevron",
    feature_fill = "#CC6677", position = position_strand(-.12)) +
  geom_seq_label() + coord_collinear(align = "center")
layout <- export_ggchord_layout(independent, include = c("feature", "ribbon"))
stopifnot(length(unique(layout$feature$layer_id)) == 2L,
  length(unique(layout$ribbon$layer_id)) == 2L)

plots <- list(composition = collinear_example, aligned = workflow_aligned,
  source_axis = workflow_source_axis, independent_layers = independent)
# These teaching fixtures are ASCII. Use standard PDF rather than the Quartz
# raster path, which can invert native ggplot2 colour bars on some macOS builds.
pdf_device <- grDevices::pdf
original <- grDevices::dev.cur()
files <- character()
for (name in names(plots)) {
  for (extension in c("png", "pdf", "svg")) {
    file <- file.path(out, paste0(name, ".", extension))
    ggplot2::ggsave(file, plots[[name]], width = 10, height = 7, dpi = 120,
      device = if (extension == "pdf") pdf_device else extension)
    stopifnot(file.exists(file), file.info(file)$size > 0,
      identical(grDevices::dev.cur(), original))
    if (extension == "png") stopifnot(
      identical(dim(png::readPNG(file))[1:2], c(840L, 1200L)))
    files <- c(files, file)
  }
}
write.csv(data.frame(file = basename(files), bytes = file.info(files)$size),
  file.path(out, "manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
cat("Collinear PNG/PDF/SVG outputs and independent layers passed.\n")
