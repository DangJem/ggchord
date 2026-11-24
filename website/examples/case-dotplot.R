library(ggchord)
library(ggplot2)

observed <- flagship_dotplot_example
sequences <- observed$seq_data
overview <- ggchord(sequences) +
  geom_alignment_segment(aes(colour = pident, linewidth = length),
    data = observed$links, alpha = 0.95) +
  scale_linewidth_continuous(
    range = c(0.72, 1.65), breaks = c(1000, 5000, 10000, 20000),
    labels = scales::label_comma(), name = "HSP length (bp)",
    guide = guide_legend(nrow = 1, order = 2,
      override.aes = list(colour = "#465F58"),
      theme = theme(legend.key.width = grid::unit(12, "mm")))
  ) +
  layout_dotplot(x_order = sequences$accver, y_order = sequences$accver,
    x_gap = observed$axis_gap, y_gap = observed$axis_gap) +
  scale_colour_stepsn(
    name = "Identity (%)", colours = pal_pident(), limits = c(0, 100),
    breaks = c(0, 50, 80, 90, 95, 100),
    guide = guide_coloursteps(direction = "horizontal", order = 1,
      even.steps = FALSE, show.limits = TRUE,
      theme = theme(legend.key.width = grid::unit(170, "mm"),
        legend.key.height = grid::unit(3.3, "mm")))
  ) +
  coord_dotplot() +
  theme_ggchord_publication(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
    legend.box = "vertical", legend.box.just = "center",
    legend.box.spacing = grid::unit(4, "mm"), legend.spacing.y = grid::unit(2, "mm"),
    plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_rect(fill = "#FFFCF5", colour = NA),
    plot.title.position = "panel", plot.title = element_text(size = 19, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11.5, colour = "#52645D", margin = margin(b = 10)),
    axis.text = element_text(size = 10.5, colour = "#465F58"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
    axis.text.x.top = element_text(size = 12, face = "bold", margin = margin(b = 6)),
    axis.text.y.left = element_text(size = 12, face = "bold", margin = margin(r = 6)),
    plot.margin = margin(8, 9, 7, 9))

# Put accessions at cell centres and number the original reference endpoints.
first <- c(1, head(cumsum(sequences$length + observed$axis_gap), -1) + 1)
last <- first + sequences$length - 1
centres <- (first + last) / 2
position_breaks <- as.vector(rbind(first, last))
position_labels <- as.vector(rbind(rep("1", length(first)), scales::label_comma()(sequences$length)))
counts <- observed$cell_labels

overview +
  scale_x_continuous("Query position (bp)", breaks = position_breaks, labels = position_labels,
    guide = guide_axis(n.dodge = 2),
    sec.axis = dup_axis(name = NULL, breaks = centres, labels = sequences$accver,
      guide = guide_axis(n.dodge = 1))) +
  scale_y_continuous("Subject sequence", breaks = centres, labels = sequences$accver,
    sec.axis = dup_axis(name = "Subject position (bp)", breaks = position_breaks,
      labels = position_labels, guide = guide_axis(check.overlap = TRUE))) +
  geom_text(aes(x = x, y = y, label = label), data = dplyr::filter(counts, HSP > 0),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.45, colour = "#284A3D", fontface = "bold") +
  geom_text(aes(x = x, y = y, label = label), data = dplyr::filter(counts, HSP == 0),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.1, colour = "#697F75") +
  labs(title = "Alignment dotplot",
    subtitle = "32 HSPs across six observed sequence pairs · One query–subject ordering per pair")
