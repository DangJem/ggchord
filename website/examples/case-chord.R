library(ggchord)
library(ggplot2)
library(dplyr)
library(scales)

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

# Follow the local tangent and reserve room inside each backbone.
sector_gap <- 21.6
start_angle <- 50
sequence_labels <- function(size, offset) {
  geom_seq_label(aes(label = accver), size = size,
    label_layout = layout_label_fixed(orientation = "tangent", side = "negative",
      offset = offset))
}
overview <- ggchord(comparison$seq_data) +
  geom_link_ribbon(aes(fill = pident), data = comparison$links, curvature = 1, avoid = "smooth",
    alpha = 0.46, linewidth = 0.16) +
  geom_seq_arrow(aes(colour = accver), linewidth = 1.1) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = c("+" = -1.65, "-" = 1.65)), thickness = 2.05, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_fill_pident() +
  scale_for(scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"),
    guide = guide_legend(nrow = 1, order = 1)), layers = "genes") +
  layout_sequences(
    orientation = comparison$orientation, angular_gap = sector_gap,
    curvature = setNames(c(0, 1, -0.65, 0.75), names(sequence_colours)),
    radius = setNames(c(1.08, 0.82, 1.18, 1.02), names(sequence_colours))
  ) +
  coord_chord(start_angle = start_angle, fit = "labels", clip = "off") +
  theme_ggchord_publication(base_size = 12) +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    legend.key.width = unit(11.5, "mm"), legend.key.height = unit(7.5, "mm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    ggchord.axis.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 21),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    plot.title.position = "plot", plot.margin = margin(12, 16, 10, 16))

# The first two labelled CDS on these backbones have adjacent anchors.
# Let the text extend away from the gap between their attachment points.
nearby_labels <- labels |>
  filter(accver %in% c("MT108731.1", "OQ646790.1")) |>
  group_by(accver) |>
  slice_min(start, n = 2, with_ties = FALSE) |>
  arrange(start, .by_group = TRUE) |>
  mutate(text_side = if_else(row_number() == 1L, "left", "right")) |>
  ungroup()
outside_labels <- function(records, hjust = 0.5) {
  geom_gene_label(aes(label = label), data = records, anchor_layer = "genes", hjust = hjust,
    segment_colour = "#666D64", segment_linewidth = 0.23,
    label_layout = layout_label_repel(arrangement = "local", orientation = "horizontal",
      side = "positive", offset = 10, padding = 1.1, fit = "wrap_shrink", min_size = 3.15,
      wrap = 18, max_lines = 3), size = 3.9)
}
overview +
  outside_labels(filter(nearby_labels, text_side == "left"), hjust = 1) +
  outside_labels(filter(nearby_labels, text_side == "right"), hjust = 0) +
  outside_labels(anti_join(labels, nearby_labels, by = "feature_id")) +
  sequence_labels(size = 4.35, offset = 17) +
  scale_seq_continuous(name = NULL, minor_breaks = NULL,
    labels = label_comma()) +
  guides(fill = guide_coloursteps(position = "left", direction = "vertical",
    even.steps = FALSE, show.limits = TRUE,
    theme = theme(legend.key.width = unit(3.6, "mm"),
      legend.key.height = unit(145, "mm"), legend.text = element_text(size = 11)))) +
  labs(title = "Four genomes, one connected view",
    subtitle = "Mixed backbone geometry · CDS sampled by source strand")
