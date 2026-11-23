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

cross_ribbons <- lapply(seq_len(nrow(comparison$cross_links)), function(i) {
  record <- comparison$cross_links[i, ]
  order <- match(c(record$qaccver, record$saccver), comparison$seq_data$accver)
  bulge <- 0.07 + 0.045 * abs(diff(order)) + 0.025 * (min(order) - 1)
  geom_link_ribbon(aes(fill = pident), data = record, bulge = bulge, avoid = "smooth",
    alpha = 0.33, colour = "#496B5D", linewidth = 0.33, linetype = "dashed")
})
overview <- ggchord(comparison$seq_data) +
  geom_link_ribbon(aes(fill = pident), data = comparison$adjacent_links, avoid = "smooth",
    alpha = 0.50) +
  cross_ribbons +
  geom_seq_arrow(aes(colour = accver), linewidth = 1.1) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = 2.2), thickness = 2.2, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_fill_pident() +
  scale_for(scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"),
    guide = guide_legend(nrow = 1, order = 1)), layers = "genes") +
  layout_sequences(align = "center", row_gap = 41) +
  coord_collinear(fit = "labels", clip = "off") +
  theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    legend.box.spacing = unit(2.5, "mm"),
    plot.title = element_text(face = "bold", size = 17),
    plot.margin = margin(10, 14, 10, 14))

# Keep adjacent names and source coordinates readable on each of the four rows.
gene_names <- lapply(c("positive", "negative"), function(side) {
  geom_gene_label_repel(aes(label = label), data = filter(labels, label_side == side),
    anchor_layer = "genes", segment_colour = "#767A73", segment_linewidth = 0.18,
    label_layout = label_layout_repel(arrangement = "local", orientation = "horizontal",
      side = side, offset = 7, padding = 1.2, fit = "wrap_shrink", min_size = 2.7,
      wrap = 23, max_lines = 3), size = 3.35)
})
overview + gene_names +
  geom_seq_label(aes(label = accver), anchor = "start", hjust = 1,
    label_layout = label_layout_fixed(orientation = "horizontal", along = -1.2, offset = 0),
    size = 3.7) +
  annotation_seq_scalebar(length = 10000) +
  scale_seq_continuous(name = NULL, minor_breaks = NULL,
    labels = label_comma(), guide = guide_seq_axis(side = "negative")) +
  guides(fill = guide_coloursteps(position = "left", direction = "vertical",
    even.steps = FALSE, show.limits = TRUE,
    theme = theme(legend.key.height = unit(93, "mm"), legend.key.width = unit(2.8, "mm")))) +
  labs(title = "Sequence architecture, across four genomes",
    subtitle = "Selected labelled CDS · strand-separated arrows · routed cross-row comparisons")
