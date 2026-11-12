# Shared grid fragments for site-label rendering and physical measurement.
ggchord_site_label_parts <- function(enzyme, coordinate, size, family = '',
    enzyme_face = 'plain', lineheight = 1.2, colour = NULL) {
  common <- list(fontsize = size*(72.27/25.4), fontfamily = family,
    lineheight = lineheight, col = colour)
  a <- grid::textGrob(enzyme, gp = do.call(grid::gpar, c(common, list(fontface = enzyme_face))))
  b <- grid::textGrob(coordinate, gp = do.call(grid::gpar, c(common, list(fontface = 'plain'))))
  has_a <- !is.na(enzyme) && nzchar(enzyme)
  has_b <- !is.na(coordinate) && nzchar(coordinate)
  gap <- grid::unit(if (has_a && has_b) .65 else 0, 'mm')
  list(enzyme = a, coordinate = b, gap = gap, has_enzyme = has_a, has_coordinate = has_b)
}

ggchord_site_label_metrics <- function(enzyme, coordinate, size, family = '',
    enzyme_face = 'plain', lineheight = 1.2, units_per_inch = 25.4) {
  n <- length(enzyme)
  size <- rep_len(size, n); family <- rep_len(family, n)
  enzyme_face <- rep_len(enzyme_face, n); lineheight <- rep_len(lineheight, n)
  width <- height <- numeric(n)
  close <- if (grDevices::dev.cur() > 1L) function() invisible(NULL) else
    ggchord_measurement_device(text = c(enzyme, coordinate))
  on.exit(close())
  for (i in seq_len(n)) {
    part <- ggchord_site_label_parts(enzyme[i], coordinate[i], size[i], family[i],
      enzyme_face[i], lineheight[i])
    a_w <- if (part$has_enzyme) grid::grobWidth(part$enzyme) else grid::unit(0, 'mm')
    b_w <- if (part$has_coordinate) grid::grobWidth(part$coordinate) else grid::unit(0, 'mm')
    a_h <- if (part$has_enzyme) grid::grobHeight(part$enzyme) else grid::unit(0, 'mm')
    b_h <- if (part$has_coordinate) grid::grobHeight(part$coordinate) else grid::unit(0, 'mm')
    width[i] <- grid::convertWidth(a_w+b_w+part$gap, 'inches', valueOnly = TRUE)*units_per_inch
    height[i] <- max(grid::convertHeight(a_h, 'inches', valueOnly = TRUE),
      grid::convertHeight(b_h, 'inches', valueOnly = TRUE))*units_per_inch
  }
  data.frame(width = width, height = height)
}
