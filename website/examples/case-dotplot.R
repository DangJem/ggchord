library(ggchord)
library(ggplot2)

observed <- flagship_dotplot_example
sequences <- observed$seq_data
overview <- ggchord(sequences) +
  geom_alignment_segment(aes(colour = pident, linewidth = length),
    data = observed$links, alpha = 0.9) +
  scale_linewidth_continuous(range = c(0.8, 1.8), name = "HSP length (bp)") +
  layout_dotplot(x_order = sequences$accver, y_order = sequences$accver,
    x_gap = observed$axis_gap, y_gap = observed$axis_gap) +
  coord_dotplot() +
  theme_ggchord_publication() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
    plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_rect(fill = "#FFFCF5", colour = NA),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(12, 16, 12, 16))

overview +
  geom_text(aes(x = x, y = y, label = label), data = observed$cell_labels,
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.8, colour = "#566E67") +
  labs(title = "Six observed sequence pairs · 32 HSPs",
    subtitle = "One query→subject direction per pair; empty cells contain no input records")
