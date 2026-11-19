library(ggchord)

dotplot_pair_counts <- flagship_dotplot_example$pair_counts
dotplot_cells <- flagship_dotplot_example$cell_labels
flagship_hsp <- flagship_dotplot_example$links
flagship_sequences <- flagship_dotplot_example$seq_data
flagship_paper <- "#FFFCF5"
flagship_theme <- theme_ggchord_publication() +
  theme(
    plot.background = element_rect(fill = flagship_paper, colour = NA),
    panel.background = element_rect(fill = flagship_paper, colour = NA),
    legend.background = element_rect(fill = flagship_paper, colour = NA),
    legend.key = element_rect(fill = flagship_paper, colour = NA),
    legend.position = "right", legend.box = "vertical", legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.margin = margin(12, 16, 12, 16)
  )
flag_dotplot <- ggchord(prepare_ggchord_data(flagship_sequences, flagship_hsp)) +
  geom_alignment_segment(aes(link_colour = pident, link_linewidth = length), alpha = 0.85) +
  scale_link_linewidth_continuous(range = c(0.8, 1.8), name = "HSP length (bp)") +
  coord_dotplot(x_seqs = flagship_sequences$accver, y_seqs = flagship_sequences$accver,
    gap = 4000) +
  ggplot2::geom_text(data = dotplot_cells, aes(x = x, y = y, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.8, colour = "#566E67") +
  labs(title = "Six observed sequence pairs · 32 HSPs",
    subtitle = "One query→subject direction per pair; empty cells contain no input records") +
  flagship_theme +
  theme(
    panel.grid = element_blank(), legend.position = "bottom", legend.box = "horizontal",
    axis.text = element_text(size = 9),
    axis.text.x.top = element_text(size = 10, face = "bold", margin = margin(b = 9)),
    axis.text.y.left = element_text(size = 10, face = "bold", margin = margin(r = 9)),
    axis.text.y.right = element_text(size = 9),
    axis.ticks.length = unit(1.5, "mm"),
    plot.margin = margin(12, 16, 12, 16)
  )
flag_dotplot
