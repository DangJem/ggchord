library(ggchord)
library(ggplot2)
library(dplyr)

reference <- flagship_circular_example$references[["pSpCas9(BB)-2A-GFP (PX458)"]]
sequence <- as_ggchord_seq(reference$seq_data, topology = "circular")
features <- reference$features
lanes <- position_feature_stack(spacing = 3.2, direction = "negative",
  base_position = position_offset(-3.5))
overview <- ggchord(sequence) +
  geom_seq(geometry = "double", thickness = 1.1, linewidth = 0.48, colour = "#292826") +
  geom_feature(aes(fill = feature_color, shape = feature_shape), data = features,
    id = "annotations", thickness = 2.8, position = lanes, colour = "#393635",
    linewidth = 0.21) +
  scale_fill_identity() +
  scale_shape_manual(values = setNames(unique(features$feature_shape), unique(features$feature_shape)),
    guide = "none") +
  geom_seq_center_label(aes(label = paste0(accver, "\n", scales::comma(length), " bp")),
    size = 4.3, fontface = "bold") +
  coord_circular(start_angle = 90, fit = "manual", xlim = c(-1.55, 1.55),
    ylim = c(-1.19, 1.19), clip = "off") +
  theme_ggchord_plasmid()

p <- overview +
  geom_feature_label_repel(aes(label = label), anchor_layer = "annotations", size = 2.7,
    label_layout = label_layout_repel(arrangement = "local", orientation = "tangent",
      side = "negative", offset = 4, fit = "shrink", min_size = 2.2)) +
  scale_seq_position(breaks = scales::breaks_pretty(5), minor_breaks = NULL,
    labels = scales::label_comma())
if (nrow(reference$sites)) {
  # Several enzyme records can share a cut. Keep every cut and group their labels.
  site_labels <- reference$sites |>
    group_by(accver, pos) |>
    summarise(label = paste(sort(unique(label)), collapse = ", "), .groups = "drop")
  p <- p +
    geom_restriction_site(data = reference$sites, size = 1.1, colour = "#514E4C") +
    geom_site_label_repel(aes(label = label), data = site_labels, size = 2.55, colour = "#373637",
      segment_colour = "#7D7973", segment_linewidth = 0.18,
      label_layout = label_layout_repel(arrangement = "columns", side = "positive",
        offset = 5, padding = 1.2, wrap = 26, max_lines = 2))
}
if (nrow(reference$primers)) {
  p <- p +
    geom_primer(data = reference$primers, id = "primers", colour = "#685876",
      position = position_offset(-6.8), thickness = 1) +
    geom_primer_label_repel(aes(label = label), anchor_layer = "primers", size = 2.3,
      label_layout = label_layout_repel(arrangement = "local", side = "negative", offset = 2))
}
p
