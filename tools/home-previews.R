# Homepage-only views of existing examples. Geometry inputs remain source-backed.
# Source after loading ggchord and the declared shared preparations; do not crop PNGs.
home_paper <- "#FFFCF5"
home_theme <- theme_ggchord_publication() + theme(
  legend.position = "none",
  plot.background = element_rect(fill = home_paper, colour = NA),
  panel.background = element_blank(),
  plot.margin = margin(7, 7, 7, 7),
  ggchord.axis.text = element_text(size = 8),
  ggchord.axis.minor.ticks = element_blank()
)
home_relationship_theme <- home_theme + theme(
  ggchord.axis.text = element_blank(), ggchord.axis.ticks = element_blank()
)
home_axis <- scale_seq_position_continuous(
  n.breaks = 2, minor_breaks = NULL,
  labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
)
home_colours <- list(
  scale_seq_colour_manual(values = seq_palette, guide = "none"),
  scale_gene_fill_manual(values = gene_palette, guide = "none")
)
# Same Chord parent: 12 observed HSPs, 20 CDS geometries, no CDS-name layer.
home_chord <- ggchord(show_sequences, show_links) +
  geom_link_ribbon(link_avoid = "smooth", alpha = 0.55, linewidth = 0.2) +
  geom_seq(
    seq_orientation = show_optimization$seq_orientation,
    seq_curvature = show_curvature, seq_radius = show_radius,
    seq_gap = 0.06, linewidth = 1.2
  ) +
  geom_gene(data = chord_genes, position = chord_position, width = 0.05) +
  geom_seq_label(size = 4, seq_label_radius = 0.75, seq_label_hjust = 0.5) +
  home_colours + home_axis + coord_chord(rotation = 50) + home_relationship_theme
# Same selected Collinear HSPs and 16 CDS; labels belong in the complete figure.
home_collinear <- ggchord(show_sequences) +
  geom_link_ribbon(data = show_adjacent, alpha = 0.55, linewidth = 0.2) +
  geom_link_ribbon(data = show_cross, adjacent_only = FALSE,
    link_ctrl_point = show_cross_controls, alpha = 0.35, linewidth = 0.25) +
  geom_seq(seq_gap = 0.35, linewidth = 1.2) +
  geom_gene(data = show_linear_features, position = position_strand(0.08), width = 0.065) +
  geom_seq_label(size = 4) +
  home_colours + home_axis + show_linear_coord + home_relationship_theme
home_radial <- ggchord(show_sequences, show_links) +
  geom_link_ribbon(adjacent_only = FALSE, alpha = 0.55, linewidth = 0.2) +
  geom_seq(seq_orientation = show_optimization$seq_orientation, linewidth = 1.2) +
  geom_gene(data = chord_genes, position = radial_position, width = 0.065) +
  geom_seq_label(size = 4, seq_label_orientation = "horizontal") +
  home_colours + home_axis +
  coord_radial_synteny(start_angle = 45, inner_radius = 0.7) + home_relationship_theme
home_sankey <- ggchord(show_sequences) +
  geom_link_ribbon(data = show_adjacent, alpha = 0.55, linewidth = 0.2) +
  geom_link_ribbon(data = show_cross, adjacent_only = FALSE, alpha = 0.35, linewidth = 0.25) +
  geom_seq(seq_gap = 0.22, linewidth = 1.2) +
  geom_gene(data = show_linear_features, position = position_strand(0.09), width = 0.075) +
  geom_seq_label(size = 3.8) +
  home_colours + home_axis + coord_sankey(column_gap = 5.2, ratio = 0.75) + home_relationship_theme
# Existing pUC19 source; all recognised features remain, only named CDS/origin text.
home_circular_features <- find_common_features(plasmid_example_pUC19)
home_circular_labels <- home_circular_features |>
  filter(type %in% c("CDS", "rep_origin"))
home_circular_position <- position_feature_stack(base_position = position_plasmid())
home_circular <- ggchord(plasmid_example_pUC19, feature_data = home_circular_features) +
  geom_seq(seq_style = "double", linewidth = 0.8) +
  geom_feature_plasmid(position = home_circular_position) +
  geom_feature_label_repel(data = home_circular_labels, position = home_circular_position,
    size = 4, label_fit = "shrink", label_min_size = 3) +
  geom_seq_center_label(size = 4.5) +
  scale_feature_fill_manual(
    values = setNames(home_circular_features$feature_color, home_circular_features$anno),
    guide = "none"
  ) +
  coord_circular() + home_axis + home_theme
# Same 32 observed HSPs; segment-centre sequence IDs replace dense coordinate text.
home_dot_axis <- export_ggchord_layout(flag_dotplot)$axis_layout
home_x <- home_dot_axis |> filter(axis == "x")
home_y <- home_dot_axis |> filter(axis == "y")
home_dotplot <- ggchord(flagship_sequences, flagship_hsp) +
  geom_alignment_segment(aes(link_colour = pident, link_linewidth = length), alpha = 0.9) +
  scale_link_linewidth_continuous(range = c(0.9, 1.8), guide = "none") +
  coord_dotplot(x_seqs = flagship_sequences$accver, y_seqs = flagship_sequences$accver,
    gap = 4000) +
  ggplot2::scale_x_continuous(breaks = (home_x$min + home_x$max) / 2,
    labels = home_x$accver, guide = ggplot2::guide_axis(angle = 25)) +
  ggplot2::scale_y_continuous(breaks = (home_y$min + home_y$max) / 2,
    labels = home_y$accver) +
  labs(x = "Query sequence", y = "Subject sequence") +
  home_theme + theme(
    axis.title = element_text(size = 11), axis.text = element_text(size = 10),
    panel.grid.major = element_blank()
  )
# Four-panel flagship source: all 32 HSPs and 266 CDS in each overview.
# Hide coordinate text for the preview, retaining names, geometry and real endpoints.
home_fusion_theme <- home_relationship_theme + theme(plot.margin = margin(2, 2, 2, 2))
home_fusion_base <- ggchord(flagship_sequences, flagship_hsp, feature_data = flagship_features) +
  geom_link_ribbon(adjacent_only = FALSE, alpha = 0.45, linewidth = 0.12, allow_repeated = TRUE) +
  geom_seq(aes(seq_colour = accver), linewidth = 1.1) +
  geom_gene(aes(gene_fill = strand), position = flagship_position, width = 0.025) +
  geom_seq_label(size = 2.6) +
  scale_seq_colour_manual(values = flagship_colours) +
  scale_gene_fill_manual(values = gene_palette) + home_fusion_theme
home_fusion <- ggchord_fusion(home_fusion_base + coord_chord(),
  center_box = c(0.01, 0.49, 0.52, 0.99)) +
  fusion_panel(home_fusion_base + coord_collinear(), id = "linear",
    box = c(0.54, 0.99, 0.52, 0.99)) +
  fusion_panel(home_fusion_base + coord_radial_synteny(start_angle = 45), id = "radial",
    box = c(0.01, 0.49, 0.01, 0.48)) +
  fusion_panel(flagship_detail + home_fusion_theme, id = "detail",
    box = c(0.54, 0.99, 0.01, 0.48)) +
  fusion_ribbon(flagship_hsp |> slice(2:5), qpanel = "center", spanel = "linear",
    fill = "#269B82", alpha = 0.22, allow_repeated = TRUE) +
  fusion_ribbon(flagship_hsp |> slice(6:9), qpanel = "center", spanel = "radial",
    fill = "#D59B17", alpha = 0.22, allow_repeated = TRUE) +
  fusion_ribbon(flagship_selected, qpanel = "center", spanel = "detail",
    fill = "#E64B35", alpha = 0.24, allow_repeated = TRUE) +
  theme(plot.background = element_rect(fill = home_paper, colour = NA))
