library(ggchord)
library(ggplot2)
library(dplyr)

comparison <- flagship_comparison_example$chord
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
  geom_link_ribbon(aes(fill = pident), data = comparison$links,
    alpha = 0.52, colour = "#53776A", linewidth = 0.16) +
  geom_seq(aes(colour = accver), linewidth = 1.1,
    arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed")) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = 1.25), thickness = 1.25, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_for("genes", scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"), guide = guide_legend(nrow = 1, order = 1))) +
  layout_sequences(
    orientation = -comparison$orientation, angular_gap = 21.6,
    curvature = setNames(c(0, 1, -0.65, 0.75), names(sequence_colours)),
    radius = setNames(c(1.08, 0.82, 1.18, 1.02), names(sequence_colours))
  ) +
  coord_chord(start_angle = 50, fit = "labels") +
  theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 17),
    plot.margin = margin(10, 14, 10, 14))

overview +
  geom_gene_label_repel(aes(label = label), data = labels, anchor_layer = "genes",
    label_layout = label_layout_repel(arrangement = "local", orientation = "horizontal",
      offset = 4, wrap = 22, max_lines = 3), size = 3.0) +
  geom_seq_label(aes(label = accver), anchor = "middle",
    label_layout = label_layout_fixed(orientation = "horizontal", side = "negative", offset = 5),
    size = 3.7) +
  scale_seq_position(breaks = scales::breaks_pretty(3), minor_breaks = NULL,
    labels = scales::label_comma()) +
  guides(fill = guide_coloursteps(position = "left", direction = "vertical",
    even.steps = FALSE, show.limits = TRUE,
    theme = theme(legend.key.height = grid::unit(68, "mm")))) +
  labs(title = "Four genomes, one connected view",
    subtitle = "Mixed backbone geometry · CDS sampled by source strand")
