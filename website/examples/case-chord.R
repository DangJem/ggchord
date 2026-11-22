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

# Rotate each accession along its sector midline. Place the text just inside the path.
sector_gap <- 21.6
start_angle <- 50
sectors <- comparison$seq_data |>
  mutate(span = (360 - n() * sector_gap) * length / sum(length),
    angle = start_angle + cumsum(span + sector_gap) - sector_gap - span / 2,
    angle = (angle + 90) %% 180 - 90)
sequence_labels <- function(size, offset) {
  lapply(seq_len(nrow(sectors)), function(i) {
    geom_seq_label(aes(label = accver), data = sectors[i, ], angle = sectors$angle[i],
      label_layout = label_layout_fixed(orientation = "horizontal", side = "negative",
        offset = offset), size = size)
  })
}
overview <- ggchord(comparison$seq_data) +
  geom_link_ribbon(aes(fill = pident), data = comparison$links, curvature = 1, avoid = "smooth",
    alpha = 0.52, colour = "#53776A", linewidth = 0.16) +
  geom_seq(aes(colour = accver), linewidth = 1.1,
    arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = c("+" = -1.65, "-" = 1.65)), thickness = 2.05, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_for("genes", scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"),
    guide = guide_legend(nrow = 1, order = 1, override.aes = list(strand = c("+", "-"))))) +
  layout_sequences(
    orientation = comparison$orientation, angular_gap = sector_gap,
    curvature = setNames(c(0, 1, -0.65, 0.75), names(sequence_colours)),
    radius = setNames(c(1.08, 0.82, 1.18, 1.02), names(sequence_colours))
  ) +
  coord_chord(start_angle = start_angle, fit = "labels", clip = "off") +
  theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 17),
    plot.margin = margin(10, 14, 10, 14))

overview +
  geom_gene_label_repel(aes(label = label), data = labels, anchor_layer = "genes",
    segment_colour = "#767A73", segment_linewidth = 0.18,
    label_layout = label_layout_repel(arrangement = "local", orientation = "horizontal",
      side = "both", offset = 12, padding = 1.2, fit = "wrap_shrink", min_size = 2.65,
      wrap = 24, max_lines = 3), size = 3.3) +
  sequence_labels(size = 4.0, offset = 11) +
  scale_seq_position(name = NULL, minor_breaks = NULL,
    labels = label_comma()) +
  guides(fill = guide_coloursteps(position = "left", direction = "vertical",
    even.steps = FALSE, show.limits = TRUE,
    theme = theme(legend.key.width = unit(2.8, "mm"),
      legend.key.height = unit(91, "mm"), legend.text = element_text(size = 10)))) +
  labs(title = "Four genomes, one connected view",
    subtitle = "Mixed backbone geometry · CDS sampled by source strand")
