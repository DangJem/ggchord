# Geometry review only. Grey polygons and red measured boxes are diagnostic
# overlays, not source feature colours or new biological measurements.
library(ggplot2)
folder <- 'tools/long-term-audit'
cases <- c('pTRE_Tight_BI', 'pCAMBIA1300', 'pSpCas9_BB_2A_GFP_PX458')
for (id in cases) {
  value <- jsonlite::fromJSON(file.path(folder, paste0('v0.27.0-quality-identity-', id, '.json')))
  contact <- value$feature_contacts[value$feature_contacts$relation == 'unresolved', ]
  stopifnot(nrow(contact) == 1L)
  text <- value$text[value$text$label_id == contact$label_id, ]
  owner <- strsplit(contact$feature_id, '/', fixed = TRUE)[[1L]][1L]
  polygons <- value$feature_polygons
  polygons <- polygons[polygons$source_id == contact$feature_source_id & polygons$layer_id == owner, ]
  bounds_x <- range(c(text$xmin_mm, text$xmax_mm)) + c(-2, 2)
  bounds_y <- range(c(text$ymin_mm, text$ymax_mm)) + c(-3, 3)
  zoom <- 5
  plot <- ggplot() + geom_polygon(data = polygons,
    aes(x_mm, y_mm, group = polygon_id), fill = '#ECECE8', colour = '#596363', linewidth = .25) +
    geom_rect(data = text, aes(xmin = xmin_mm, xmax = xmax_mm, ymin = ymin_mm, ymax = ymax_mm),
      fill = NA, colour = '#B24E45', linewidth = .25) +
    geom_text(data = text, aes(centre_x_mm, centre_y_mm, label = label,
      angle = angle_radians*180/pi, size = I(size_mm*zoom)), show.legend = FALSE) +
    coord_fixed(xlim = bounds_x, ylim = bounds_y, expand = FALSE) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  ggsave(file.path(folder, paste0('v0.27.0-axis-feature-contact-', id, '.png')), plot,
    width = diff(bounds_x)*zoom/25.4, height = diff(bounds_y)*zoom/25.4,
    dpi = 150, bg = 'white', limitsize = FALSE)
  cat('Diagnostic envelopes:', id, '\n')
}
