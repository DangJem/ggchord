library(ggchord)
library(ggplot2)
library(dplyr)
library(scales)

comparison <- flagship_comparison_example$sankey
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
  geom_link_ribbon(aes(fill = pident), data = comparison$adjacent_links, avoid = "smooth",
    alpha = 0.50) +
  geom_link_ribbon(aes(fill = pident), data = comparison$cross_links, avoid = "smooth",
    alpha = 0.33, colour = "#547568", linewidth = 0.28) +
  geom_seq_arrow(aes(colour = accver), linewidth = 1.1) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = 2.1), thickness = 2.1, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_fill_pident() +
  scale_for(scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"),
    guide = guide_legend(nrow = 1, order = 1)), layers = "genes") +
  layout_sequences(column_gap = 145, row_gap = 20) +
  coord_sankey(fit = "geometry", clip = "off") +
  theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 17),
    plot.margin = margin(10, 14, 10, 14))

# Each display column points its CDS names away from the ribbons.
column_labels <- lapply(c("left", "right"), function(column) {
  rows <- labels |> mutate(instance_id = paste(column, accver, sep = ":"))
  geom_gene_label(aes(label = label), data = rows, anchor_layer = "genes", anchor = "tail",
    hjust = if (column == "left") 1 else 0, size = 3.25,
    segment_colour = "#767A73", segment_linewidth = 0.18,
    label_layout = layout_label_repel(arrangement = "local", orientation = "horizontal",
      side = if (column == "left") "negative" else "positive",
      offset = 7, padding = 1, wrap = 22, max_lines = 3))
})
sequence_labels <- function(size = 3.9) {
  lapply(c("left", "right"), function(column) {
    left <- column == "left"
    rows <- comparison$seq_data |> mutate(instance_id = paste(column, accver, sep = ":"))
    geom_seq_label(aes(label = accver), data = rows, anchor = "start",
      hjust = if (left) 1 else 0, vjust = 0,
      label_layout = layout_label_fixed(orientation = "horizontal", along = -2.5,
        side = if (left) "negative" else "positive", offset = 1.5), size = size)
  })
}
overview + column_labels + sequence_labels() +
  scale_seq_continuous(name = NULL, minor_breaks = NULL,
    labels = label_comma()) +
  guides(fill = guide_coloursteps(direction = "horizontal", even.steps = FALSE,
    show.limits = TRUE, theme = theme(legend.key.width = unit(115, "mm"), legend.key.height = unit(2.8, "mm")))) +
  theme(legend.box = "vertical") +
  labs(title = "Four genomes, mirrored relationships",
    subtitle = "Original Collinear sources; independent query and subject displays")
