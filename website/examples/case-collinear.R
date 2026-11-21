library(ggchord)
library(ggplot2)
library(dplyr)

comparison <- flagship_comparison_example$collinear
# Stable entity keys let the selected names follow their own CDS arrows.
identify_genes <- function(records) {
  records |> mutate(feature_id = paste(accver, start, end, strand, sep = "/"), label = anno)
}
genes <- identify_genes(comparison$features)
labels <- identify_genes(comparison$labels)
sequence_colours <- c(
  "MT108731.1" = "#E64B35", "MT118296.1" = "#233E4A",
  "OQ646790.1" = "#269B82", "OR222515.1" = "#D59B17"
)
strand_colours <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")

overview <- ggchord(comparison$seq_data) +
  geom_link_ribbon(aes(fill = pident), data = comparison$adjacent_links,
    alpha = 0.50, colour = "#53776A", linewidth = 0.18) +
  geom_link_ribbon(aes(fill = pident), data = comparison$cross_links,
    alpha = 0.33, colour = "#547568", linewidth = 0.28, linetype = "dashed") +
  geom_seq(aes(colour = accver), linewidth = 1.1,
    arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed")) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = 2.2), thickness = 2.2, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_for("genes", scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"), guide = guide_legend(nrow = 1, order = 1))) +
  layout_sequences(align = "center", row_gap = 24) +
  coord_collinear(fit = "labels") +
  theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 17),
    plot.margin = margin(10, 14, 10, 14))

overview +
  geom_gene_label_repel(aes(label = label), data = labels, anchor_layer = "genes",
    label_layout = label_layout_repel(arrangement = "local", orientation = "horizontal",
      offset = 4, wrap = 22, max_lines = 3), size = 3.15) +
  geom_seq_label(aes(label = accver), anchor = "start",
    label_layout = label_layout_fixed(orientation = "horizontal", side = "positive", offset = 5),
    size = 3.7) +
  annotation_seq_scalebar(length = 10000) +
  scale_seq_position(breaks = scales::breaks_pretty(3), minor_breaks = NULL,
    labels = scales::label_comma()) +
  guides(fill = guide_coloursteps(position = "left", direction = "vertical",
    even.steps = FALSE, show.limits = TRUE,
    theme = theme(legend.key.height = grid::unit(68, "mm")))) +
  labs(title = "Sequence architecture, across four genomes",
    subtitle = "Selected labelled CDS · strand-separated arrows · routed cross-row comparisons")
