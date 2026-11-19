library(ggchord)

show_curvature <- c(
  "MT108731.1" = 0, "MT118296.1" = 1,
  "OQ646790.1" = -0.65, "OR222515.1" = 0.75
)
show_radius <- c(
  "MT108731.1" = 1.08, "MT118296.1" = 0.82,
  "OQ646790.1" = 1.18, "OR222515.1" = 1.02
)
chord_genes <- flagship_comparison_example$chord$features
chord_labels <- flagship_comparison_example$chord$labels
chord_position <- position_strand(c("+" = -0.030, "-" = 0.030))
show_links <- flagship_comparison_example$chord$links
show_optimization <- list(seq_orientation = flagship_comparison_example$chord$orientation)
show_sequences <- flagship_comparison_example$chord$seq_data
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
gene_palette <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")
identity_guide_theme <- theme(
  legend.key.height = unit(100, "mm"), legend.key.width = unit(5, "mm"),
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 11), legend.margin = margin(r = 7, unit = "mm")
)
identity_guide <- ggplot2::guide_coloursteps(
  available_aes = "ribbon_fill", even.steps = FALSE, show.limits = TRUE,
  position = "left", direction = "vertical", theme = identity_guide_theme
)
seq_palette <- setNames(
  c("#E64B35", "#233E4A", "#269B82", "#D59B17"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
doc_legend_theme <- theme(
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 11),
  legend.key.width = unit(9, "mm"), legend.key.height = unit(5, "mm"),
  legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
  legend.margin = margin(l = 4, r = 4, unit = "mm")
)
sequence_guide <- guide_ggchord_legend(nrow = 1, theme = doc_legend_theme)
strand_guide <- guide_ggchord_legend(
  nrow = 1, override.aes = list(strand = c("+", "-")), theme = doc_legend_theme
)
doc_axis <- function(coordinates = "relative") {
  scale_seq_position_continuous(
    coordinates = coordinates, n.breaks = 3, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
  )
}
showcase_style <- list(
  scale_seq_colour_manual(values = seq_palette, name = "Sequence", guide = sequence_guide),
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand"),
  guides(gene_fill = strand_guide, ribbon_fill = identity_guide),
  doc_axis(), doc_theme
)
p_show_chord <- ggchord(prepare_ggchord_data(seq_data = show_sequences, ribbon_data = show_links)) +
  geom_link_ribbon(link_avoid = "smooth", alpha = 0.46, linewidth = 0.14) +
  layout_sequences(orientation = show_optimization$seq_orientation, gap = 0.06) +
  geom_seq(seq_curvature = show_curvature, seq_radius = show_radius, linewidth = 1.1) +
  geom_feature(preset = 'gene', data = chord_genes, position = chord_position, width = 0.025) +
  geom_gene_label_repel(
    data = chord_labels, label_layout = "radial", label_wrap = 24,
    size = 2.8, label_max_lines = 3, position = chord_position
  ) +
  geom_seq_label(size = 4, seq_label_radius = 0.77, seq_label_hjust = 0.5) +
  showcase_style +
  coord_chord(rotation = 50) +
  labs(
    title = "Four genomes, one connected view",
    subtitle = "Mixed backbone geometry · CDS sampled by source strand"
  )
p_show_chord
