# One explicit device and representative physical/diagnostic paths only.
# This does not invoke testthat, package checks, or an output matrix.
pkgload::load_all(quiet = TRUE)
local({
  # The source package's final NAMESPACE registers this method. Registering in
  # this isolated validation process also supports an in-progress source tree.
  registerS3method("ggplot_gtable", "ggchord_built", ggchord:::ggplot_gtable.ggchord_built,
    envir = asNamespace("ggplot2"))
  plot <- ggchord(data.frame(accver = c("A", "B"), length = c(1000, 800))) +
    geom_seq() + geom_seq_label() + coord_collinear()
  original <- plot$ggchord$ref$layout
  previous <- grDevices::dev.cur()
  report <- diagnose_ggchord(plot, width = 7, height = 4, device = "pdf", colour_vision = TRUE)
  stopifnot(report$measurement$status == "measured",
    all(report$measurement$panel_mm > 0),
    all(report$measurement$panel_mm <= report$measurement$canvas_mm),
    all(is.finite(as.matrix(report$text[c("xmin_mm", "xmax_mm", "ymin_mm", "ymax_mm")]))),
    "font_embedding_and_fallback" %in% report$checks$code,
    report$checks$status[report$checks$code == "font_embedding_and_fallback"] == "not_checked",
    identical(plot$ggchord$ref$layout, original),
    identical(grDevices::dev.cur(), previous), is.null(ggchord:::ggchord_render_context()))
  smaller <- diagnose_ggchord(plot, width = 4, height = 3, device = "pdf")
  stopifnot(smaller$measurement$canvas_mm[["width"]] < report$measurement$canvas_mm[["width"]],
    smaller$measurement$units_per_inch > report$measurement$units_per_inch,
    identical(plot$ggchord$ref$layout, original))

  # clip='off' permits text outside a panel when it still fits the canvas.
  outside <- ggchord(data.frame(accver = "A", length = 1000)) + geom_seq() +
    ggplot2::annotate("text", x = -.05, y = .2, label = "site", size = 1) +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 40, unit = "mm")) +
    coord_collinear(xlim = c(0, 2), ylim = c(-1, 1), fit = "manual", clip = "off")
  checked <- diagnose_ggchord(outside, width = 7, height = 4, device = "pdf")
  stopifnot("small_text" %in% checked$checks$code,
    !"canvas_text_bounds" %in% checked$checks$code)

  # Standard ggsave uses the same registered gtable hook, and does not retain
  # a target device in the next render's context.
  file <- tempfile(fileext = ".pdf")
  on.exit(unlink(file), add = TRUE)
  ggplot2::ggsave(file, plot = plot, device = "pdf", width = 7, height = 4)
  stopifnot(file.exists(file), file.info(file)$size > 0,
    is.null(ggchord:::ggchord_render_context()), identical(grDevices::dev.cur(), previous))
})
cat("Actual-panel context, target-size diagnostics and standard saving passed.\n")
