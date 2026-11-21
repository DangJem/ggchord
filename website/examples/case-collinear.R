library(ggchord)
library(ggplot2)
library(dplyr)
library(scales)

comparison <- flagship_comparison_example$collinear
# Stable entity keys let the selected names follow their own CDS arrows.
identify_genes <- function(records) {
  records |> mutate(feature_id = paste(accver, start, end, strand, sep = "/"), label = anno)
}
genes <- identify_genes(comparison$features)
labels <- identify_genes(comparison$labels) |>
  arrange(accver, start) |>
  group_by(accver) |>
  mutate(label_side = if_else(row_number() %% 2L == 1L, "positive", "negative")) |>
  ungroup()
sequence_colours <- c(
  "MT108731.1" = "#E64B35", "MT118296.1" = "#233E4A",
  "OQ646790.1" = "#269B82", "OR222515.1" = "#D59B17"
)
strand_colours <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")

overview <- ggchord(comparison$seq_data) +
  geom_link_ribbon(aes(fill = pident), data = comparison$adjacent_links,
    alpha = 0.50, colour = "#53776A", linewidth = 0.18) +
  geom_link_ribbon(aes(fill = pident), data = comparison$cross_links, bulge = 0.1,
    alpha = 0.33, colour = "#547568", linewidth = 0.28, linetype = "dashed") +
  geom_seq(aes(colour = accver), linewidth = 1.1,
    arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = 2.2), thickness = 2.2, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_for("genes", scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"), guide = guide_legend(nrow = 1, order = 1))) +
  layout_sequences(align = "center", row_gap = 41) +
  coord_collinear(xlim = c(-2.05, 1.50), fit = "labels", clip = "off") +
  theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 17),
    plot.margin = margin(10, 14, 10, 14))

# Keep adjacent names and source coordinates readable on each of the four rows.
gene_names <- lapply(c("positive", "negative"), function(side) {
  geom_gene_label_repel(aes(label = label), data = filter(labels, label_side == side),
    anchor_layer = "genes", segment_colour = "#767A73", segment_linewidth = 0.18,
    label_layout = label_layout_repel(arrangement = "local", orientation = "horizontal",
      side = side, offset = 7, padding = 1.2, wrap = 22, max_lines = 3), size = 3.05)
})
overview + gene_names +
  geom_seq_label(aes(label = accver), anchor = "start", hjust = 1,
    label_layout = label_layout_fixed(orientation = "horizontal", along = -1.2, offset = 0),
    size = 3.7) +
  annotation_seq_scalebar(length = 10000) +
  scale_seq_position(name = NULL, breaks = breaks_pretty(2), minor_breaks = NULL,
    labels = label_comma()) +
  guides(fill = guide_coloursteps(position = "left", direction = "vertical",
    even.steps = FALSE, show.limits = TRUE,
    theme = theme(legend.key.height = unit(93, "mm")))) +
  labs(title = "Sequence architecture, across four genomes",
    subtitle = "Selected labelled CDS · strand-separated arrows · routed cross-row comparisons")
