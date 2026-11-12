# Native composition-furniture checks; no publication/package acceptance.
pkgload::load_all(quiet = TRUE)
library(dplyr)
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
data(seq_data_example); data(ribbon_data_example)
hit <- ribbon_data_example |> slice_head(n = 1)
base <- ggchord(seq_data_example) + geom_seq(show.legend = FALSE)
plain <- ggchord_fusion(base) + fusion_panel(base + coord_collinear(), 'detail') +
  fusion_ribbon(hit, qpanel = 'center', spanel = 'detail', fill = '#548CA0', alpha = .35)
titled <- plain + ggplot2::labs(title = 'Real source intervals', subtitle = 'Independent panels; native plot furniture',
  caption = 'Source: unchanged bundled HSP records', tag = 'A') +
  theme(plot.margin = margin(5, 5, 5, 5, 'mm'), plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 10), plot.caption = element_text(size = 9),
    plot.background = element_rect(fill = '#F3F5F7', colour = NA))
a <- export_ggchord_fusion(plain, 12, 8)
b <- export_ggchord_fusion(titled, 12, 8)
stopifnot(is.null(a$metadata$plot_decoration), b$metadata$plot_decoration$status == 'native_gtable',
  identical(b$metadata$plot_decoration$labels$title, 'Real source intervals'),
  min(b$metadata$content_bounds_mm[c(1, 3)]) >= 5,
  b$metadata$content_bounds_mm[4] < b$metadata$canvas_mm[2],
  identical(a$relationship_data, b$relationship_data), nrow(b$connections) == nrow(a$connections))
# Public source-position mapping independently checks final outer-mm offsets.
close_mapping_device <- ggchord_measurement_device(12, 8)
for (id in c('center', 'detail')) {
  front <- filter(b$fronts, panel_id == id)
  panel <- filter(b$panels, panel_id == id)
  context <- b$diagnostics[[id]]$render_context
  original <- if (id == 'center') base else base + coord_collinear()
  position <- if (id == 'detail') {
    map_ggchord_position(original, unique(front$instance_id), front$position_boundary, source = TRUE)
  } else list(x = front$local_x, y = front$local_y) # Chord has no public inverse-position API.

  expected_x <- panel$panel_xmin+(position$x-context$x_range[1])/diff(context$x_range)*(panel$panel_xmax-panel$panel_xmin)
  expected_y <- panel$panel_ymin+(position$y-context$y_range[1])/diff(context$y_range)*(panel$panel_ymax-panel$panel_ymin)
  stopifnot(all(abs(front$x-expected_x) < .02), all(abs(front$y-expected_y) < .02))
}
close_mapping_device()
# Inspect actual native text in the assembled drawing, not only metadata.
collect_text <- function(g) {
  children <- c(g$grobs, as.list(g$children))
  c(if (inherits(g, 'text')) as.character(g$label), unlist(lapply(children, collect_text), use.names = FALSE))
}
close <- ggchord_measurement_device(12, 8)
grid::grid.newpage(); built <- ggchord_build_fusion(titled)
text <- collect_text(built$grob)
stopifnot(all(c('Real source intervals', 'Independent panels; native plot furniture',
  'Source: unchanged bundled HSP records', 'A') %in% text))
close()
# A later labs/theme update follows native replacement/blank semantics.
updated <- titled + ggplot2::labs(title = 'Updated title', subtitle = NULL, caption = NULL, tag = NULL) +
  theme(plot.title = element_blank())
u <- export_ggchord_fusion(updated, 12, 8)
stopifnot(identical(u$metadata$plot_decoration$labels$title, 'Updated title'))
# No meaningful fake composition axes are created by a complete native theme.
classic <- titled + ggplot2::theme_classic() + ggplot2::theme(plot.margin = margin(4, 4, 4, 4, 'mm'))
c <- export_ggchord_fusion(classic, 12, 8)
stopifnot(c$metadata$plot_decoration$coordinate_axes == 'none; child axes belong to each child plot')
# Parent links resolve the same nested leaf/source after outer title translation.
child <- plain + ggplot2::labs(title = 'Nested source panels') +
  theme(plot.background = element_rect(fill = '#EAF1F5', colour = NA))
parent <- ggchord_fusion(base, center_box = c(.02, .32, .15, .85)) +
  fusion_panel(child, 'nested', box = c(.4, .98, .05, .95)) +
  fusion_ribbon(hit, qpanel = 'center', spanel = 'nested/detail', fill = '#548CA0') +
  ggplot2::labs(title = 'Parent composition', caption = 'Nested titles retain source/instance registration')
n <- export_ggchord_fusion(parent, 16, 9)
stopifnot(nrow(n$connections) == 2, 'nested/detail' %in% n$panels$panel_id,
  n$diagnostics$nested$plot_decoration$labels$title == 'Nested source panels')
# Explicit overlap uses the shared disjoint rendering ledger, including a native
# composition background; one translucent ribbon must not be repainted.
overlap <- titled; overlap$allow_overlap <- TRUE
path <- 'tools/long-term-audit/v0.27.0-fusion-decoration.png'
previous <- grDevices::dev.cur()
ggchord_diagnostic_device(path, 'png', 12, 8, 110)
grid::grid.newpage(); rendered <- ggchord_build_fusion(overlap)
grid::grid.draw(rendered$grob)
grDevices::dev.off()
if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
check_single_paint <- function(rendered, path) {
  d <- rendered$export$paths
  qn <- sum(rendered$export$fronts$endpoint == 'q')
  sn <- sum(rendered$export$fronts$endpoint == 's')
  edge_n <- (nrow(d)-qn-sn)/2
  edge1 <- as.matrix(d[nrow(d)-seq_len(edge_n)+1L, c('x', 'y')])
  edge2 <- as.matrix(d[qn+seq_len(edge_n), c('x', 'y')])
  centers <- (edge1+edge2)/2
  panels <- rendered$export$panels |> filter(coordinate != 'fusion')
  outside <- vapply(seq_len(nrow(centers)), function(i) !any(
    centers[i, 1] >= panels$container_xmin & centers[i, 1] <= panels$container_xmax &
    centers[i, 2] >= panels$container_ymin & centers[i, 2] <= panels$container_ymax), logical(1))
  candidates <- which(outside & rowSums((edge1-edge2)^2) > 4)
  stopifnot(length(candidates) > 0)
  at <- candidates[ceiling(length(candidates)/2)]
  raster <- png::readPNG(path)
  canvas <- rendered$export$metadata$canvas_mm
  col <- round(centers[at, 1]/canvas[1]*dim(raster)[2])
  row <- round((1-centers[at, 2]/canvas[2])*dim(raster)[1])
  observed <- raster[row, col, 1:3]
  # Quartz can colour-manage nominal RGB. Compare one paint on the SAME
  # device implementation rather than locking raw sRGB output values.
  reference <- tempfile(fileext = '.png')
  ggchord_diagnostic_device(reference, 'png', 1, 1, 110)
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = '#F3F5F7', col = NA))
  grid::grid.rect(gp = grid::gpar(fill = scales::alpha('#548CA0', .35), col = NA))
  grDevices::dev.off()
  calibration <- png::readPNG(reference)
  expected <- calibration[ceiling(dim(calibration)[1]/2), ceiling(dim(calibration)[2]/2), 1:3]
  unlink(reference)
  stopifnot(max(abs(observed-expected)) < .02)
}
if (requireNamespace('png', quietly = TRUE)) {
  check_single_paint(rendered, path)
  # Exercise the older nested-connection viewport path too. Parent foregrounds
  # must not disable the ancestor's disjoint clipping rectangles.
  child_only <- plain + ggplot2::labs(title = 'Nested opacity check') +
    theme(plot.background = element_rect(fill = '#F3F5F7', colour = NA))
  nested_only <- ggchord_fusion(base, allow_overlap = TRUE, center_box = c(.02, .28, .2, .8)) +
    fusion_panel(child_only, 'child', box = c(.34, .98, .05, .95))
  nested_path <- tempfile(fileext = '.png')
  ggchord_diagnostic_device(nested_path, 'png', 16, 9, 110)
  grid::grid.newpage(); nested_render <- ggchord_build_fusion(nested_only)
  grid::grid.draw(nested_render$grob)
  grDevices::dev.off()
  check_single_paint(nested_render, nested_path)
  unlink(nested_path)
}
# Fixed decorations and margins cannot be silently reduced to claim a fit.
fail(export_ggchord_fusion(titled + theme(plot.margin = margin(100, 100, 100, 100, 'mm')), 4, 4),
  'insufficient composition space')
jsonlite::write_json(list(status = 'passed', native_title_subtitle_caption_tag = TRUE,
  actual_drawing_text = TRUE, later_updates_preserved = TRUE, no_dummy_axes = TRUE,
  source_records_unchanged = TRUE, nested_source_registration = TRUE,
  transparent_connection_once = requireNamespace('png', quietly = TRUE),
  explicit_margin_failure = TRUE, content_bounds_mm = b$metadata$content_bounds_mm),
  'tools/long-term-audit/v0.27.0-decoration-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Fusion native plot decorations: passed\n')
