# One explicit typography control, not a biological observation or release matrix.
pkgload::load_all(quiet = TRUE)
path <- 'tools/long-term-audit/v0.27.0-site-label-metrics.png'
ggchord_diagnostic_device(path, 'png', 5, 3, 110)
d <- data.frame(x = .5, y = .5, .component = 'label', source_row = 1,
  label = 'Long\nname (123)', enzyme_label = 'Long\nname', coordinate_label = '(123)',
  enzyme_fontface = 'bold', label_order = 'enzyme_position',
  colour = '#202020', alpha = 1, size = 3.2, family = 'sans', fontface = 1,
  angle = 0, hjust = .5, vjust = .5, lineheight = 1.2,
  .ggchord_composite_label = TRUE)
coordinate <- list(transform = function(data, panel) data)
grob <- GeomRestrictionSite$draw_panel(d, NULL, coordinate)
grid::grid.newpage()
grid::pushViewport(grid::viewport(width = grid::unit(100, 'mm'),
  height = grid::unit(50, 'mm'), xscale = c(0, 1), yscale = c(0, 1)))
grid::grid.draw(grob)
texts <- function(grob) {
  if (inherits(grob, 'text')) return(list(grob))
  unlist(lapply(grob$children %||% list(), texts), recursive = FALSE)
}
parts <- texts(grob)
stopifnot(length(parts) == 2L)
bounds <- lapply(parts, function(part) {
  x <- grid::convertX(part$x, 'mm', valueOnly = TRUE)
  y <- grid::convertY(part$y, 'mm', valueOnly = TRUE)
  w <- grid::convertWidth(grid::grobWidth(part), 'mm', valueOnly = TRUE)
  h <- grid::convertHeight(grid::grobHeight(part), 'mm', valueOnly = TRUE)
  c(xmin = x-part$hjust*w, xmax = x+(1-part$hjust)*w,
    ymin = y-part$vjust*h, ymax = y+(1-part$vjust)*h)
})
actual <- c(xmin = min(vapply(bounds, `[[`, 0, 'xmin')),
  xmax = max(vapply(bounds, `[[`, 0, 'xmax')),
  ymin = min(vapply(bounds, `[[`, 0, 'ymin')),
  ymax = max(vapply(bounds, `[[`, 0, 'ymax')))
context <- list(panel_mm = c(width = 100, height = 50), panel_bounds_mm = c(xmin = 0, ymin = 0))
styled <- ggchord_diagnostic_text_styles(d, geom_restriction_site(
  data = data.frame(accver = 'format', enzyme = 'format control', position = 1)))
diagnostic <- ggchord_diagnostic_text(styled, coordinate, NULL, context)
reported <- unlist(diagnostic[c('xmin_mm', 'xmax_mm', 'ymin_mm', 'ymax_mm')], use.names = FALSE)
native <- d
native$text <- native$label; native$text_x <- 50; native$text_y <- 25; native$text_angle <- 0
box <- ggchord_text_boxes(native, units_per_inch = 25.4)
obstacle <- unlist(box[c('xmin', 'xmax', 'ymin', 'ymax')], use.names = FALSE)
print(rbind(rendered = unname(actual), diagnostic = reported, obstacle = obstacle))
stopifnot(isTRUE(all.equal(unname(actual), reported)),
  isTRUE(all.equal(unname(actual), obstacle)))
grid::popViewport(); grDevices::dev.off()
jsonlite::write_json(list(status = 'passed one wrapped composite typography control',
  actual_rendered_fragment_envelope_mm = as.list(actual), diagnostic_envelope_mm = as.list(reported),
  obstacle_envelope_mm = as.list(obstacle), file = path,
  png_sha256 = digest::digest(file = path, algo = 'sha256'),
  scope = 'Logical grid text envelopes on one PNG device; not exact glyph ink, embedded fonts or platform certification'),
  'tools/long-term-audit/v0.27.0-site-label-metrics-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Rendered site fragments, diagnostic bounds and routing bounds agree\n')
