library(ggchord)
library(ggplot2)

observed <- flagship_dotplot_example
sequences <- observed$seq_data
overview <- ggchord(sequences) +
  geom_alignment_segment(aes(colour = pident, linewidth = length),
    data = observed$links, alpha = 0.9) +
  scale_linewidth_continuous(range = c(0.95, 2.0), name = "HSP length (bp)") +
  layout_dotplot(x_order = sequences$accver, y_order = sequences$accver,
    x_gap = observed$axis_gap, y_gap = observed$axis_gap) +
  scale_colour_link_pident() +
  coord_dotplot() +
  theme_ggchord_publication() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
    plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_rect(fill = "#FFFCF5", colour = NA),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10.5),
    axis.title = element_text(size = 11.5),
    legend.title = element_text(size = 10.5), legend.text = element_text(size = 10),
    axis.text.x.top = element_text(size = 11, face = "bold", margin = margin(b = 5)),
    axis.text.y.left = element_text(size = 11, face = "bold", margin = margin(r = 5)),
    plot.margin = margin(8, 11, 8, 11))

overview +
  geom_text(aes(x = x, y = y, label = label), data = observed$cell_labels,
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2, colour = "#465F58") +
  labs(title = "Six observed sequence pairs · 32 HSPs",
    subtitle = "One query→subject direction per pair; empty cells contain no input records")
