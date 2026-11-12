# Component style and physical packing controls; no package test suite.
pkgload::load_all(quiet = TRUE)
path <- 'tools/long-term-audit/v0.27.0-site-layout-style.png'
ggchord_diagnostic_device(path, 'png', 5, 3, 110)
source <- data.frame(accver = 'format', length = 1000)
sites <- data.frame(accver = 'format', enzyme = 'format control', position = 1)
p <- ggchord(source) + geom_restriction_site(data = sites) +
  theme_ggchord_plasmid() +
  theme(ggchord.restriction.label = element_text(size = 13, family = 'sans', lineheight = 1.4))
layer <- p$layers[[1L]]
style <- ggchord_resolve_theme_component(p, layer, 'text_params', 'ggchord.restriction.label')
stopifnot(isTRUE(all.equal(style$size, 13/ggplot2::.pt)), style$family == 'sans',
  style$lineheight == 1.4, is.null(layer$geom_params$text_params$size))
explicit <- geom_restriction_site(data = sites, label_size = 2.5, fontface = 'italic', family = 'serif')
fixed <- ggchord_resolve_theme_component(p, explicit, 'text_params', 'ggchord.restriction.label')
stopifnot(fixed$size == 2.5, fixed$family == 'serif', fixed$fontface == 'italic')
check <- function(style, composite) {
  label <- data.frame(x = 1, y = 1, source_row = 1, restriction_component = 'label',
    label = 'Long\nname (123)', size = style$size, family = style$family,
    fontface = style$fontface, lineheight = style$lineheight,
    .ggchord_composite_label = composite, enzyme_label = 'Long\nname',
    coordinate_label = '(123)', enzyme_fontface = 'bold', label_order = 'enzyme_position',
    label_side = 'outside', hjust = .5, vjust = .5)
  leader <- label[rep(1L, 2L), ]; leader$restriction_component <- 'leader'
  leader$x <- c(.7, 1); leader$y <- c(.7, 1)
  d <- rbind(label, leader)
  tasks <- ggchord_perimeter_tasks(list(format = list(restriction_site = d)), 25.4)
  stopifnot(nrow(tasks) == 1L, identical(tasks$fontface, as.character(style$fontface)))
  actual <- if (composite) ggchord_site_label_metrics(tasks$enzyme_label,
    tasks$coordinate_label, style$size, style$family, 'bold', style$lineheight) else {
    grob <- grid::textGrob(tasks$text, gp = grid::gpar(fontsize = style$size*72.27/25.4,
      fontfamily = style$family, fontface = style$fontface, lineheight = style$lineheight))
    data.frame(width = grid::convertWidth(grid::grobWidth(grob), 'mm', TRUE),
      height = grid::convertHeight(grid::grobHeight(grob), 'mm', TRUE))
  }
  stopifnot(isTRUE(all.equal(tasks$w, actual$width+.03*25.4)),
    isTRUE(all.equal(tasks$h, actual$height+.03*25.4)))
  list(composite = composite, effective_style = style,
    packing_mm = c(width = tasks$w, height = tasks$h), actual_mm = unlist(actual))
}
a <- check(style, TRUE); b <- check(fixed, FALSE)
grid::grid.newpage(); grid::grid.text('Physical packing and explicit style controls')
grDevices::dev.off()
jsonlite::write_json(list(status = 'passed', implicit = a, explicit_uniform_face = b,
  original_layer_preserved = TRUE,
  scope = 'One PNG composite packing envelope and component precedence; platforms unverified'),
  'tools/long-term-audit/v0.27.0-site-layout-style-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Site packing and component style precedence: passed\n')
