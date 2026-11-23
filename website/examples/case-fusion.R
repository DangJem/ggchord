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
detail_border <- element_rect(fill = NA, colour = "#927DA2", linewidth = 0.35)
panel_theme <- theme_ggchord_publication() +
  theme(plot.background = element_rect(fill = paper, colour = NA),
    panel.background = element_blank(), legend.position = "bottom",
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.5),
    legend.key.spacing.x = unit(2.0, "mm"),
    legend.key.spacing.y = unit(0.8, "mm"),
    legend.margin = margin(3, 4, 3, 4),
    plot.title.position = "plot", plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "#52646A"),
    plot.margin = margin(4, 5, 4, 5))

# Each panel uses the original source coordinates, including the windowed detail.
relationship_view <- function(sequences, links, genes) {
  ggchord(sequences) +
    geom_link_ribbon(aes(fill = pident), data = links, alpha = 0.42,
      linewidth = 0.22, colour = "#516960", curvature = 0.55, avoid = "smooth") +
    geom_seq(aes(colour = accver), linewidth = 1, key_glyph = draw_key_path, show.legend = c(colour = TRUE),
      arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
    geom_gene(aes(fill = strand), data = genes, id = "genes",
      position = position_strand(1.3), thickness = 1.55, linewidth = 0.17, show.legend = c(fill = TRUE, colour = FALSE)) +
    scale_colour_manual(name = "Sequence", values = sequence_colours,
      limits = names(sequence_colours), drop = FALSE, guide = guide_legend(nrow = 2, order = 2)) +
    scale_fill_pident() +
    scale_for(scale_fill_manual(name = "Strand", values = strand_colours,
      limits = names(strand_colours), drop = FALSE, guide = guide_legend(nrow = 1, order = 1,
        override.aes = list(strand = c("+", "-")))), layers = "genes") +
    scale_seq_continuous(guide = "none") +
    panel_theme +
    # The shared Identity ruler fits its labels without scaling to the full canvas.
    render_policy(mode = "fixed")
}
overview <- relationship_view(study$seq_data, study$links, study$features)
seq_names <- geom_seq_label(aes(label = accver), anchor = "start", hjust = 1, size = 3.05,
  label_layout = label_layout_fixed(orientation = "horizontal", along = -1.5, offset = 0))
chord <- overview +
  geom_seq_label(aes(label = accver), size = 3.05,
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 4)) +
  layout_sequences(angular_gap = 7, orientation = study$orientation) +
  coord_chord(start_angle = 45, fit = "labels") +
  labs(title = "A  |  Chord", subtitle = "6 HSPs; 3 CDS per genome")
include_endpoints <- function(breaks) {
  function(limits) {
    inner <- breaks(limits)
    sort(unique(c(limits, inner[inner > limits[1] & inner < limits[2]])))
  }
}
overview_axis <- scale_seq_continuous(name = NULL, breaks = include_endpoints(breaks_extended(n = 3)),
  labels = label_comma(), minor_breaks = NULL, guide = guide_seq_axis(side = "negative"))
overview_xlim <- c(-1.42, 1.12)
draw_linear <- function(row_gap) {
  overview + seq_names +
    layout_sequences(row_gap = row_gap, align = "center", orientation = study$orientation) +
    coord_collinear(xlim = overview_xlim, ylim = c(-0.62, 0.06), fit = "manual", clip = "off") +
    overview_axis +
    labs(title = "B  |  Collinear", subtitle = "The same intervals in linear coordinates")
}
spoke_names <- function(sequences, hjust) {
  geom_seq_label(aes(label = accver), data = sequences, anchor = "end", size = 3.05,
    hjust = hjust,
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 0, along = 4.5))
}
# With the 45-degree start, the first two spokes point right and the other two left.
radial <- overview +
  spoke_names(slice(study$seq_data, 1:2), hjust = 0) +
  spoke_names(slice(study$seq_data, 3:4), hjust = 1) +
  layout_sequences(inner_radius = 0.26, orientation = study$orientation) +
  coord_radial_synteny(start_angle = 45, fit = "labels") +
  labs(title = "C  |  Radial", subtitle = "Shared intervals; A / D")

draw_context <- function(row_gap) {
  relationship_view(study$context_sequences, study$selected, study$context_features) +
    seq_names + layout_sequences(row_gap = row_gap, align = "center") +
    overview_axis +
    coord_collinear(xlim = overview_xlim, ylim = c(-0.77, 0.16), fit = "manual", clip = "off") +
    labs(title = "D  |  Two-genome context", subtitle = "Full-length sequences; windows enlarged in E")
}
windows <- study$windows |> mutate(instance_id = accver)
detail_genes <- study$all_features |>
  semi_join(study$detail_features, by = c("accver" = "source_accver", "feature_id"))
# Place neighbouring CDS on alternate sides of each local backbone.
detail_genes <- detail_genes |>
  arrange(accver, start) |>
  group_by(accver) |>
  mutate(label_side = if_else(row_number() %% 2L == 1L, "negative", "positive")) |>
  ungroup()
detail_label <- function(features, side, hjust = 0) {
  geom_gene_label_repel(aes(label = anno), data = features,
    anchor_layer = "genes", size = 2.7, hjust = hjust, segment_colour = "#70766B",
    segment_linewidth = 0.18,
    label_layout = label_layout_repel(arrangement = "local", side = side,
      orientation = "horizontal", offset = 6, padding = 1.2, wrap = 19, max_lines = 3))
}
# Place the labels for the F-linked row above its backbone, clear of the zoom leaders.
linked_label <- detail_genes |>
  filter(accver == study$base_cds$accver, label_side == "negative")
detail_labels <- list(
  detail_label(filter(detail_genes, label_side == "negative", accver != study$base_cds$accver), "negative"),
  detail_label(filter(detail_genes, label_side == "positive"), "positive"),
  detail_label(linked_label, "positive", hjust = 1)
)
draw_zoom <- function(row_gap, top, axis_text_offset) {
  relationship_view(study$context_sequences, study$selected, detail_genes) +
    seq_names +
    detail_labels +
    layout_sequences(windows = windows, row_gap = row_gap, align = "center") +
    coord_collinear(ylim = c(-0.92, top), fit = "labels", clip = "off") +
    scale_seq_continuous(name = NULL, coordinates = "source",
      breaks = include_endpoints(breaks_width(4000, offset = 2000)), labels = label_comma(),
      minor_breaks = NULL, guide = guide_seq_axis(side = "negative")) +
    theme(panel.border = detail_border,
      ggchord.axis.text.offset = unit(axis_text_offset, "mm"),
      plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    labs(title = "E  |  Local zoom", subtitle = "Source windows + 800 bp flanks")
}

cells <- study$base_cells
base_range <- range(cells$alignment_column)
base_ticks <- cells |>
  filter(alignment_column %in% seq(base_range[1], base_range[2], by = 10))
draw_bases <- function(size) {
  # Fixed aspect keeps the bordered cells square. The source ruler shares their frame.
  result <- ggplot(cells, aes(alignment_column, accver, label = symbol)) +
    geom_residue(aes(fill = symbol), size = size, fontface = "bold",
      cell_colour = "#64566F", cell_linewidth = 0.24,
      overflow = "error", show.legend = FALSE) +
    scale_fill_residue("DNA", guide = "none") +
    annotate("text", x = base_range[1] - 2, y = 1, label = "DNA", size = 2.5, colour = "#25363B") +
    annotate("segment", x = base_ticks$alignment_column, xend = base_ticks$alignment_column,
      y = 0.5, yend = 0.41, colour = "#657177", linewidth = 0.22) +
    annotate("text", x = base_ticks$alignment_column, y = 0.33,
      label = base_ticks$source_position, vjust = 1, size = 2.35, colour = "#25363B") +
    annotate("text", x = mean(base_range), y = -0.67,
      label = "Source position (bp)", size = 2.5, colour = "#25363B") +
    scale_x_continuous(guide = "none") +
    scale_y_discrete(guide = "none") +
    coord_fixed(xlim = base_range + c(-3.75, 0.8), ylim = c(-1.10, 1.65), expand = FALSE) +
    labs(title = "F  |  Nucleotide detail", x = NULL, y = NULL,
      subtitle = paste0(study$base_cds$accver, " / ", study$base_cds$feature_id,
        "  |  first 60 coding bases (", study$base_cds$strand, ")  |  ",
        study$base_start, "-", study$base_end, " bp")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
      panel.border = detail_border,
      plot.background = element_rect(fill = paper, colour = NA),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, colour = "#52646A", hjust = 0.5),
      plot.margin = margin(3, 5, 3, 5))
  register_ggchord_alignment(result, cells)
}

areas <- tribble(
  ~panel,     ~xmin, ~xmax, ~ymin, ~ymax,
  "chord",    0.01,  0.46,  0.62,  0.99,
  "linear",   0.465, 0.99,  0.615, 0.99,
  "radial",   0.01,  0.46,  0.20,  0.60,
  "center",   0.465, 0.99,  0.405, 0.605,
  "zoom",     0.48,  0.99,  0.17,  0.395,
  "bases",    0.025, 0.98,  0.005, 0.126
)
connections <- bind_rows(
  mutate(study$relation_ab, qpanel_path = "chord", spanel_path = "linear"),
  mutate(study$relation_ac, qpanel_path = "chord", spanel_path = "radial"),
  mutate(study$relation_cd, qpanel_path = "radial", spanel_path = "center")
)
draw_fusion <- function(preview = FALSE) {
  layout <- areas
  views <- list(
    chord = chord,
    linear = draw_linear(row_gap = if (preview) 16 else 27),
    radial = radial,
    center = draw_context(row_gap = if (preview) 20 else 33),
    zoom = draw_zoom(row_gap = if (preview) 18 else 36,
      top = if (preview) 0.48 else 0.26, axis_text_offset = if (preview) 0.8 else 2.2),
    bases = draw_bases(size = if (preview) 4.15 else 5.4)
  )
  if (preview) {
    layout$ymin[layout$panel == "chord"] <- 0.60
    layout$ymax[layout$panel == "radial"] <- 0.58
    layout$ymin[layout$panel == "zoom"] <- 0.16
    layout$ymax[layout$panel == "bases"] <- 0.145
  }
  output <- ggchord_fusion(views, layout = layout_fusion(layout), guides = "collect") +
    fusion_ribbon(aes(fill = relation), data = filter(connections, relation == "A - B"),
      curvature = 0.09, alpha = 0.38, colour = "#536571", linewidth = 0.26) +
    fusion_ribbon(aes(fill = relation), data = filter(connections, relation != "A - B"),
      curvature = 0.20, alpha = 0.55, colour = "#536571", linewidth = 0.26) +
    scale_fill_manual(name = if (preview) "Links" else "Correspondence",
      values = c("A - B" = "#527FB8", "A - C" = "#CE9741", "C - D" = "#9064AC")) +
    fusion_zoom("center", "zoom", study$windows, colour = "#725094",
      fill = NA, alpha = 0.62, linewidth = 0.30, connector = "border") +
    fusion_zoom("zoom", "bases", filter(study$focus_regions, detail_panel == "bases"),
      colour = "#725094", fill = NA, alpha = 0.62, linewidth = 0.30, connector = "border") +
    theme(legend.position = "bottom", legend.title = element_text(size = 10.5, face = "bold"),
      legend.text = element_text(size = 9.5),
      legend.key.spacing.x = unit(2.0, "mm"), legend.margin = margin(3, 4, 3, 4),
      plot.background = element_rect(fill = paper, colour = NA),
      plot.margin = margin(8, 12, 8, 12),
      plot.title = element_text(face = "bold", size = 20))
  if (preview) return(output)
  output +
    labs(title = "One dataset, six linked views",
      subtitle = "Six selected HSPs and twelve source CDS per overview; full inputs retained separately",
      caption = "Coloured ribbons: observed HSPs. Source windows connect D to E; the selected CDS bases connect E to F.")
}

draw_fusion()
