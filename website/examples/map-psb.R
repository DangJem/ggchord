library(ggchord)

reference <- flagship_circular_example$references[["pSB1C3"]]
lanes <- position_feature_stack(spacing = 0.085, base_position = position_offset(-0.1))
p <- ggchord(reference$seq_data) +
  geom_seq(seq_style = "double") +
  geom_feature(data = reference$features, preset = "plasmid", position = lanes) +
  geom_feature_label(data = reference$features, placement = "repel", position = lanes,
    overflow = "shrink", label_min_size = 2.5) +
  geom_seq_center_label() +
  coord_circular() +
  theme_ggchord_plasmid() +
  geom_restriction_site(data = reference$sites, label_size = 2.5)
if (nrow(reference$primers) > 0) {
  p <- p +
    geom_primer(data = reference$primers) +
    geom_primer_label_repel(data = reference$primers, size = 2.5)
}
p
