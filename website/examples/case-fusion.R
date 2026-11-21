library(ggchord)
library(ggplot2)
library(dplyr)
library(scales)

study <- flagship_fusion_example
paper <- "#FFFCF5"
sequence_colours <- c(
  "MT108731.1" = "#E64B35", "MT118296.1" = "#233E4A",
  "OQ646790.1" = "#269B82", "OR222515.1" = "#D59B17"
)
strand_colours <- c("+" = "#CE1536", "-" = "#346B9D")
panel_theme <- theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = paper, colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    plot.title.position = "plot", plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9, colour = "#52646A"),
    plot.margin = margin(4, 5, 4, 5))

# Each panel uses the original source coordinates, including the windowed detail.
relationship_view <- function(sequences, links, genes) {
  ggchord(sequences) +
    geom_link_ribbon(aes(fill = pident), data = links, alpha = 0.35, linewidth = 0.1, curvature = 1) +
    geom_seq(aes(colour = accver), linewidth = 1, key_glyph = draw_key_path, show.legend = c(colour = TRUE),
      arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
    geom_gene(aes(fill = strand), data = genes, id = "genes",
      position = position_strand(1.3), thickness = 1.2, linewidth = 0.15, show.legend = c(fill = TRUE, colour = FALSE)) +
    scale_colour_manual(name = "Sequence", values = sequence_colours,
      limits = names(sequence_colours), drop = FALSE, guide = guide_legend(nrow = 2, order = 2)) +
    scale_for("genes", scale_fill_manual(name = "Strand", values = strand_colours,
      limits = names(strand_colours), drop = FALSE, guide = guide_legend(nrow = 1, order = 1))) +
    guides(fill = guide_coloursteps(position = "bottom", direction = "horizontal",
      even.steps = FALSE, show.limits = TRUE,
      theme = theme(legend.key.width = unit(59, "mm"), legend.key.height = unit(3, "mm")))) +
    scale_seq_position(guide = "none") +
    panel_theme
}
overview <- relationship_view(study$seq_data, study$links, study$features)
seq_names <- geom_seq_label(aes(label = accver), anchor = "start", hjust = 1, size = 2.7,
  label_layout = label_layout_fixed(orientation = "horizontal", along = -1.5, offset = 0))
chord <- overview +
  geom_seq_label(aes(label = accver), size = 2.6,
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 4)) +
  layout_sequences(angular_gap = 7, orientation = study$orientation) +
  coord_chord(start_angle = 45, fit = "labels") +
  labs(title = "A  |  Chord", subtitle = "6 HSPs; 3 CDS per genome")
linear <- overview + seq_names +
  layout_sequences(row_gap = 9, align = "center", orientation = study$orientation) +
  coord_collinear(xlim = c(-1.65, 1.12), fit = "labels", clip = "off") +
  scale_seq_position(name = NULL, breaks = breaks_pretty(2), labels = label_comma(), minor_breaks = NULL) +
  labs(title = "B  |  Collinear", subtitle = "The same intervals in linear coordinates")
radial <- overview +
  geom_seq_label(aes(label = accver), anchor = "end", size = 2.7,
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 3)) +
  layout_sequences(inner_radius = 0.26, orientation = study$orientation) +
  coord_radial_synteny(start_angle = 45, fit = "labels") +
  labs(title = "C  |  Radial", subtitle = "Shared intervals; A / D")

context <- relationship_view(study$context_sequences, study$selected, study$context_features) +
  seq_names + layout_sequences(row_gap = 12, align = "center") +
  scale_seq_position(name = NULL, breaks = breaks_pretty(2), labels = label_comma(), minor_breaks = NULL) +
  coord_collinear(xlim = c(-1.65, 1.12), fit = "labels", clip = "off") +
  labs(title = "D  |  Two-genome context", subtitle = "Full-length sequences; windows enlarged in E")
windows <- study$windows |> mutate(instance_id = accver)
detail_genes <- study$all_features |>
  semi_join(study$detail_features, by = c("accver" = "source_accver", "feature_id"))
# Place neighbouring CDS on alternate sides of each local backbone.
detail_genes <- detail_genes |>
  arrange(accver, start) |>
  group_by(accver) |>
  mutate(label_side = if_else(row_number() %% 2L == 1L, "negative", "positive")) |>
  ungroup()
detail_labels <- lapply(c("negative", "positive"), function(side) {
  geom_gene_label_repel(aes(label = anno), data = filter(detail_genes, label_side == side),
    anchor_layer = "genes", size = 2.45, hjust = 0, segment_colour = "#70766B",
    segment_linewidth = 0.18,
    label_layout = label_layout_repel(arrangement = "local", side = side,
      orientation = "horizontal", offset = 6, padding = 1.2, wrap = 19, max_lines = 3))
})
zoom <- relationship_view(study$context_sequences, study$selected, detail_genes) +
  seq_names +
  detail_labels +
  layout_sequences(windows = windows, row_gap = 19, align = "center") +
  coord_collinear(xlim = c(-1.65, 1.12), fit = "labels", clip = "off") +
  scale_seq_position(name = NULL, coordinates = "source", breaks = breaks_width(4000, offset = 2000),
    labels = label_comma(), minor_breaks = NULL) +
  theme(panel.border = element_rect(fill = NA, colour = "#A09AA3", linewidth = 0.3)) +
  labs(title = "E  |  Local zoom", subtitle = "Source windows + 800 bp flanks")

cells <- study$base_cells
bases <- ggplot(cells, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 5.1, fontface = "bold",
    overflow = "error", show.legend = FALSE) +
  scale_fill_residue("DNA") +
  scale_x_continuous(breaks = seq(1, 60, 10),
    labels = function(x) cells$source_position[match(x, cells$alignment_column)],
    expand = expansion(add = 0.6)) +
  scale_y_discrete(labels = "DNA", expand = expansion(add = 1.30)) +
  labs(title = "F  |  Nucleotide detail", x = "Source position (bp)", y = NULL,
    subtitle = paste0(study$base_cds$accver, " / ", study$base_cds$feature_id,
      "  |  first 60 coding bases (", study$base_cds$strand, ")  |  ",
      study$base_start, "-", study$base_end, " bp")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#A09AA3", linewidth = 0.3),
    plot.background = element_rect(fill = paper, colour = NA),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, colour = "#52646A", hjust = 0.5))
bases <- register_ggchord_alignment(bases, cells)

panels <- list(chord = chord, linear = linear, radial = radial,
  center = context, zoom = zoom, bases = bases)
areas <- tribble(
  ~panel,     ~xmin, ~xmax, ~ymin, ~ymax,
  "chord",    0.01,  0.46,  0.62,  0.99,
  "linear",   0.51,  0.99,  0.68,  0.99,
  "radial",   0.01,  0.46,  0.20,  0.60,
  "center",   0.48,  0.99,  0.50,  0.67,
  "zoom",     0.48,  0.99,  0.20,  0.49,
  "bases",    0.025, 0.98,  0.005, 0.17
)
connections <- bind_rows(
  mutate(study$relation_ab, qpanel_path = "chord", spanel_path = "linear"),
  mutate(study$relation_ac, qpanel_path = "chord", spanel_path = "radial"),
  mutate(study$relation_cd, qpanel_path = "radial", spanel_path = "center")
)
fusion <- ggchord_fusion(panels, layout = layout_fusion(areas), guides = "collect") +
  fusion_ribbon(aes(fill = relation), data = connections, curvature = 0.20,
    alpha = 0.50, colour = "#536571", linewidth = 0.22) +
  scale_fill_manual(name = "Correspondence",
    values = c("A - B" = "#527FB8", "A - C" = "#CE9741", "C - D" = "#9064AC")) +
  fusion_zoom("center", "zoom", study$windows, colour = "#8D8797",
    fill = NA, alpha = 0.45, linewidth = 0.22) +
  fusion_zoom("zoom", "bases", filter(study$focus_regions, detail_panel == "bases"),
    colour = "#8D8797", fill = NA, alpha = 0.45, linewidth = 0.22) +
  theme(legend.position = "bottom", legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.background = element_rect(fill = paper, colour = NA),
    plot.margin = margin(8, 12, 8, 12),
    plot.title = element_text(face = "bold", size = 20))
fusion +
  labs(title = "One dataset, six linked views",
    subtitle = "Six selected HSPs and twelve source CDS per overview; full inputs retained separately",
    caption = "Coloured ribbons: observed HSPs. Source windows connect D to E; the selected CDS bases connect E to F.")
