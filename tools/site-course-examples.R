# Shared website course snippets; original figure recipes remain unchanged.
# ---- coursecomparison ----
comparison <- ggchord(radial_seqs, radial_hsp, feature_data = radial_features) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon(adjacent_only = FALSE) +
  geom_feature(shape = "arrow", width = 0.025)
p_coursechord <- comparison + coord_chord()
p_courselinear <- comparison + coord_collinear()
p_courseradial <- comparison + coord_radial_synteny()
p_coursesankey <- comparison + coord_sankey()
# ---- coursecomparisonreport ----
comparison_plots <- list(
  Chord = p_coursechord, Collinear = p_courselinear,
  Radial = p_courseradial, Sankey = p_coursesankey
)
comparison_reports <- lapply(comparison_plots, function(plot) {
  export_ggchord_layout(plot, original_data = TRUE)
})
names(comparison_reports)
# ---- coursesankeygap ----
p_coursesankeygap <- p_foundationsankey + coord_sankey(column_gap = 0.8)
p_coursesankeygap
# ---- courseradialangle ----
p_courseradialangle <- p_courseradial +
  coord_radial_synteny(start_angle = 45)
p_courseradialangle
# ---- coursefusion ----
course_center <- ggchord(fusion_seqs |> slice(1)) +
  geom_seq(show.legend = FALSE) +
  coord_circular()
course_detail <- ggchord(fusion_seqs) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  coord_collinear(seq_axis = "vertical")
p_coursefusion <- ggchord_fusion(course_center) +
  fusion_panel(course_detail, id = "detail", side = "right") +
  fusion_ribbon(
    fusion_links |> slice(1), qpanel = "center", spanel = "detail",
    fill = "#269B82", alpha = 0.35
  )
p_coursefusion
# ---- courseend ----
