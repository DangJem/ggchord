library(ggchord)

show_linear_features <- flagship_comparison_example$collinear$features
show_adjacent <- flagship_comparison_example$collinear$adjacent_links
show_cross <- flagship_comparison_example$collinear$cross_links
show_linear_coord <- coord_collinear(align = "center", orientation = "auto", ratio = 0.7)
show_cross_controls <- flagship_comparison_example$collinear$route_controls
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
p_show_collinear <- ggchord(show_sequences) +
  geom_link_ribbon(
    aes(ribbon_fill = pident),
    data = show_adjacent,
    link_avoid = "smooth", alpha = 0.48, linewidth = 0.22
  ) +
  geom_link_ribbon(
    aes(ribbon_fill = pident),
    data = show_cross,
    adjacent_only = FALSE, link_ctrl_point = show_cross_controls,
    link_avoid = "smooth", alpha = 0.28, linewidth = 0.3
  ) +
  layout_sequences(gap = 0.42) +
  geom_seq(linewidth = 1) +
  geom_feature(preset = 'gene', data = show_linear_features, position = position_strand(0.09), width = 0.045) +
  geom_gene_label_repel(
    data = flagship_comparison_example$collinear$labels,
    label_layout = "auto", label_wrap = 24,
    label_side = "auto", position = position_strand(0.09), size = 3.3
  ) +
  geom_seq_label(size = 4) +
  geom_seq_scale(length = 10000, offset = unit(c(0, 10), "mm")) +
  showcase_style +
  show_linear_coord +
  labs(
    title = "Sequence architecture, across four genomes",
    subtitle = "Selected labelled CDS · strand-separated arrows · routed cross-row comparisons"
  )
p_show_collinear
