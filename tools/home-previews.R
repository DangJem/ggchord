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
  geom_seq_label(aes(label = accver), size = 3.2,
    label_layout = label_layout_fixed(orientation = "tangent", side = "negative", offset = 4)))
home_collinear <- home_relationship("collinear",
  geom_seq_label(aes(label = accver), anchor = "start", size = 3.1,
    label_layout = label_layout_fixed(orientation = "horizontal", along = -9, offset = 0))) +
  layout_sequences(align = "center", row_gap = 16)
home_radial <- home_relationship("radial",
  geom_seq_label(aes(label = accver), anchor = "end", size = 3.1,
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 2)))
home_sankey <- home_relationship("sankey",
  geom_seq_label(aes(label = accver), anchor = "start", size = 3.1,
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 2))) +
  layout_sequences(column_gap = 71, row_gap = 9)
home_circular <- with(flagship_cases$circular, overview +
  geom_feature_label_repel(aes(label = label), data = filter(features, type %in% c("CDS", "rep_origin")),
    anchor_layer = "annotations", size = 3.1, segment_colour = "#77736E", segment_linewidth = 0.22,
    label_layout = label_layout_repel(arrangement = "local", orientation = "horizontal", offset = 5)) +
  scale_seq_position(guide = "none") + home_theme)
# Coordinates in the inter-sequence gaps have no source bases. Centre each accession on its own span.
home_dotplot <- with(flagship_cases$dotplot, {
  offsets <- c(0, head(cumsum(sequences$length + observed$axis_gap), -1))
  midpoints <- offsets + (1 + sequences$length) / 2
  overview +
    scale_x_continuous(breaks = midpoints, labels = sequences$accver, guide = guide_axis(angle = 25)) +
    scale_y_continuous(breaks = midpoints, labels = sequences$accver) +
    labs(x = "Query sequence", y = "Subject sequence") +
    home_theme + theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10))
})
home_fusion <- flagship_cases$fusion$fusion
