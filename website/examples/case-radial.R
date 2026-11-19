library(ggchord)

chord_genes <- flagship_comparison_example$chord$features
radial_cross_bulge <- 0.20
radial_labels <- flagship_comparison_example$radial$labels
radial_orientation <- flagship_comparison_example$radial$orientation
radial_position <- position_strand(c("+" = -0.055, "-" = 0.055))
gene_palette <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")
relationship_identity <- ggplot2::guide_coloursteps(
  available_aes = "ribbon_fill", even.steps = FALSE, show.limits = TRUE,
  position = "bottom", direction = "horizontal",
  theme = theme(
    legend.key.width = unit(130, "mm"), legend.key.height = unit(4, "mm"),
    legend.title = element_text(size = 10), legend.text = element_text(size = 9),
    legend.margin = margin(2, 2, 2, 2)
  )
)
relationship_guides <- guides(
  seq_colour = guide_ggchord_legend(nrow = 1, position = "bottom"),
  gene_fill = guide_ggchord_legend(
    nrow = 1, position = "bottom", override.aes = list(strand = c("+", "-"))
  ),
  ribbon_fill = relationship_identity
)
doc_background <- "#FFFCF5"
doc_theme <- theme_ggchord_publication(
  axis.text = element_text(size = 11),
  axis.text.check.overlap = TRUE
) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(9, "mm"), legend.key.height = unit(5, "mm"),
    legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
    legend.box.spacing = unit(6, "mm"),
    plot.background = element_rect(fill = doc_background, colour = NA),
    panel.background = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 10, colour = "#5C625A"),
    plot.margin = margin(14, 18, 14, 18)
  )
relationship_theme <- doc_theme + theme(
  legend.position = "bottom", legend.box = "vertical",
  legend.title = element_text(size = 10, face = "bold"),
  legend.text = element_text(size = 9),
  legend.box.spacing = unit(3, "mm"), legend.spacing.y = unit(2, "mm"),
  legend.margin = margin(2, 2, 2, 2),
  ggchord.axis.text = element_text(size = 9),
  plot.margin = margin(8, 10, 8, 10),
  plot.title = element_text(size = 17), plot.subtitle = element_text(size = 10)
)
seq_palette <- setNames(
  c("#E64B35", "#233E4A", "#269B82", "#D59B17"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
relationship_style <- list(
  scale_seq_colour_manual(values = seq_palette, name = "Sequence"),
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand"),
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)),
  relationship_theme, relationship_guides
)
show_links <- flagship_comparison_example$chord$links
show_sequences <- flagship_comparison_example$chord$seq_data
flag_radial <- ggchord(prepare_ggchord_data(show_sequences, show_links)) +
  geom_link_ribbon(adjacent_only = FALSE, link_avoid = "smooth", alpha = 0.46, linewidth = 0.20) +
  layout_sequences(orientation = radial_orientation, gap = 0.3) +
  geom_seq(linewidth = 1.1) +
  geom_feature(preset = 'gene', data = chord_genes, position = radial_position, width = 0.065) +
  geom_gene_label_repel(data = radial_labels, label_side = "auto", label_wrap = 20,
    label_max_lines = 3, position = radial_position, size = 3.1) +
  geom_seq_label(size = 4.4, seq_label_orientation = "horizontal",
    seq_label_vjust = c("MT108731.1" = 0.5, "MT118296.1" = 0.5,
      "OQ646790.1" = 1.5, "OR222515.1" = 0.5)) +
  relationship_style +
  coord_radial_synteny(start_angle = 45, inner_radius = 0.7,
    nonadjacent_bulge = radial_cross_bulge) +
  labs(title = "Four genomes, a radial view", subtitle = "Original Chord sources; all 20 CDS annotations")
flag_radial
