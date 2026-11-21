# Source tools/site-flagships.R first. The preview uses each case's original observation layers.
home_paper <- "#FFFCF5"
home_theme <- theme(
  legend.position = "none", plot.margin = margin(9, 10, 9, 10),
  plot.background = element_rect(fill = home_paper, colour = NA),
  panel.background = element_blank()
)
home_relationship <- function(view, names_layer) {
  flagship_cases[[view]]$overview + names_layer + scale_seq_position(guide = "none") + home_theme
}
home_chord <- home_relationship("chord",
  flagship_cases$chord$sequence_labels(size = 3.6, offset = 13)) +
  coord_chord(start_angle = 50, fit = "manual", xlim = c(-1.43, 1.20), ylim = c(-1.46, 1.14))
home_collinear <- home_relationship("collinear",
  geom_seq_label(aes(label = accver), anchor = "start", size = 3.45, hjust = 1,
    label_layout = label_layout_fixed(orientation = "horizontal", along = -1.2, offset = 0))) +
  layout_sequences(align = "center", row_gap = 16) +
  coord_collinear(xlim = c(-2.18, 1.50), fit = "labels", clip = "off")
home_radial <- home_relationship("radial", flagship_cases$radial$spoke_labels(size = 3.55)) +
  coord_radial_synteny(start_angle = 45, fit = "manual",
    xlim = c(-3.70, 3.35), ylim = c(-3.30, 3.20), clip = "off")
home_sankey <- home_relationship("sankey",
  geom_seq_label(aes(label = accver), anchor = "start", size = 3.1,
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 2))) +
  layout_sequences(column_gap = 80, row_gap = 5) +
  theme(plot.margin = margin(35, 9, 34, 9))
home_circular <- with(flagship_cases$circular,
  draw_reference_map(reference, overview = TRUE) + home_theme)
# Coordinates in the inter-sequence gaps have no source bases. Centre each accession on its own span.
home_dotplot <- with(flagship_cases$dotplot, {
  offsets <- c(0, head(cumsum(sequences$length + observed$axis_gap), -1))
  midpoints <- offsets + (1 + sequences$length) / 2
  overview +
    scale_x_continuous(breaks = midpoints, labels = sequences$accver, guide = guide_axis(angle = 25)) +
    scale_y_continuous(breaks = midpoints, labels = sequences$accver) +
    labs(x = "Query sequence", y = "Subject sequence") +
    home_theme + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10),
      axis.text.y.left = element_text(face = "plain"),
      panel.grid.major = element_blank())
})
home_fusion <- flagship_cases$fusion$draw_fusion(preview = TRUE)
