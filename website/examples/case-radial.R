library(ggchord)
library(ggplot2)
library(dplyr)
library(scales)

comparison <- flagship_comparison_example$radial
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

spokes <- comparison$seq_data$accver
relationships <- comparison$links |>
  mutate(distance = abs(match(qaccver, spokes) - match(saccver, spokes)),
    route_class = if_else(distance %in% c(1, length(spokes) - 1L), "neighbour", "across"))
# Sequence names sit beyond the spoke tips, facing the nearest plot margin.
spoke_angles <- (45 - (seq_along(spokes) - 1) * 360 / length(spokes)) * pi / 180
spokes_on_left <- cos(spoke_angles) < 0
spoke_labels <- function(size, along = 0) {
  lapply(c("left", "right"), function(side) {
    indices <- spokes_on_left == (side == "left")
    geom_seq_label(aes(label = accver), data = comparison$seq_data[indices, ],
      anchor = "end", hjust = if (side == "left") 1.08 else -0.08, size = size,
      label_layout = label_layout_fixed(orientation = "horizontal", along = along, offset = 1))
  })
}
overview <- ggchord(comparison$seq_data) +
  geom_link_ribbon(aes(fill = pident, linetype = route_class), data = relationships, curvature = 1, bulge = 0.20,
    alpha = 0.48, colour = "#3F695D", linewidth = 0.27) +
  scale_linetype_manual(values = c(neighbour = "solid", across = "dashed"), guide = "none") +
  geom_seq(aes(colour = accver), linewidth = 1.1,
    arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  geom_gene(aes(fill = strand), data = genes, id = "genes",
    position = position_strand(offset = c("+" = -2.25, "-" = 2.25)), thickness = 2.25, colour = "#454842",
    linewidth = 0.15) +
  scale_colour_manual(name = "Sequence", values = sequence_colours,
    guide = guide_legend(nrow = 1, order = 2)) +
  scale_for("genes", scale_fill_manual(name = "Strand", values = strand_colours,
    breaks = c("+", "-"), guide = guide_legend(nrow = 1, order = 1))) +
  layout_sequences(orientation = comparison$orientation, inner_radius = 0.39) +
  coord_radial_synteny(start_angle = 45, fit = "labels", clip = "off") +
  theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = "#FFFCF5", colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 17),
    plot.margin = margin(10, 14, 10, 14))

# Alternate CDS names across a spoke, then align the nearer text edge with its leader.
gene_names <- labels |>
  arrange(accver, start) |>
  group_by(accver) |>
  mutate(side = if_else(row_number() %% 2L == 1L, "positive", "negative")) |>
  ungroup() |>
  mutate(normal_x = -sin(spoke_angles[match(accver, spokes)]),
    read_left = if_else(side == "positive", normal_x < 0, normal_x > 0)) |>
  group_by(side, read_left) |>
  group_map(function(records, placement) {
    geom_gene_label_repel(aes(label = label), data = records, anchor_layer = "genes",
      hjust = as.integer(placement$read_left), size = 3.0,
      segment_colour = "#767A73", segment_linewidth = 0.18,
      label_layout = label_layout_repel(arrangement = "local", orientation = "horizontal",
        side = placement$side, offset = 4, padding = 1, wrap = 20, max_lines = 3))
  })
overview + gene_names +
  spoke_labels(size = 3.7, along = 6.5) +
  scale_seq_position(name = NULL, minor_breaks = NULL,
    labels = label_comma()) +
  guides(fill = guide_coloursteps(direction = "horizontal", even.steps = FALSE,
    show.limits = TRUE, theme = theme(legend.key.width = unit(115, "mm")))) +
  theme(legend.box = "vertical") +
  labs(title = "Four genomes, a radial view",
    subtitle = "Original Chord sources; all 20 CDS annotations")
