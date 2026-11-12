# Explicit typography observations, no biological records or release matrix.
pkgload::load_all(quiet = TRUE)
target <- FALSE
native <- TRUE
prefix <- 'v0.27.0-native-text-metrics-explicit-family'
path <- paste0('tools/long-term-audit/', prefix, '.png')
ggchord_diagnostic_device(path, 'png', 6, 4, 110)
grid::grid.newpage()
labels <- c('Long name (123)', 'Long\nname (123)', 'Format / 格式')
observations <- list()
for (i in seq_along(labels)) {
  face <- if (i == 2L) 'bold' else 'plain'
  d <- data.frame(text = labels[i], text_x = .5, text_y = .8-i*.2,
    text_angle = 0, size = 3.2, family = 'Songti SC', fontface = face,
    hjust = .5, vjust = .5, lineheight = 1.2)
  grob <- grid::textGrob(d$text, x = d$text_x, y = d$text_y,
    gp = grid::gpar(fontsize = d$size*72.27/25.4, fontfamily = d$family,
      fontface = d$fontface, lineheight = d$lineheight))
  grid::grid.draw(grob)
  actual <- c(width = grid::convertWidth(grid::grobWidth(grob), 'mm', TRUE),
    height = grid::convertHeight(grid::grobHeight(grob), 'mm', TRUE))
  packed <- if (target) native_candidate(d, units_per_inch = 25.4) else
    ggchord_text_boxes(d, units_per_inch = 25.4)
  reported <- c(width = packed$w, height = packed$h)
  observations[[i]] <- list(text = d$text, family = d$family, face = face,
    actual_mm = actual, layout_mm = reported, difference_mm = reported-actual,
    equal = isTRUE(all.equal(actual, reported)))
}
grDevices::dev.off()
if (native) {
  for (label in labels) {
    before <- grDevices::dev.cur()
    fallback <- ggchord_text_boxes(data.frame(text = label, text_x = 0, text_y = 0,
      text_angle = 0, size = 3.2), units_per_inch = 25.4)
    stopifnot(identical(grDevices::dev.cur(), before),
      all(is.finite(fallback$w)), all(fallback$w > 0),
      all(is.finite(fallback$h)), all(fallback$h > 0))
  }
}
if (target || native) stopifnot(all(vapply(observations, `[[`, logical(1), 'equal')))
jsonlite::write_json(list(status = if (native) 'passed native target-device type controls' else if (target) 'passed target-device prototype type controls' else 'observed, not an assumed pass',
  production_override = target, observations = observations,
  source_sha256 = digest::digest(file = 'R/helpers-label-boxes.R', algo = 'sha256'),
  scope = 'Three explicit type specimens using installed macOS Songti SC; rectangle metrics plus separate visual inspection, not ink outlines, embedding or font portability. Font remains system-installed and is not bundled or redistributed; no additional licence rights are asserted.',
  file = path, png_sha256 = digest::digest(file = path, algo = 'sha256')),
  paste0('tools/long-term-audit/', prefix, '-results.json'), pretty = TRUE, auto_unbox = TRUE)
cat('Native text metric observations recorded\n')
