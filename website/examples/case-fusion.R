library(ggchord)
library(ggplot2)

flagship_paper <- "#FFFCF5"
flagship_sequences <- flagship_dotplot_example$seq_data
flagship_hsp <- flagship_dotplot_example$links
flagship_features <- flagship_fusion_example$all_features
flagship_events <- flagship_fusion_example$events
flagship_colours <- c(
  "MT108731.1" = "#E64B35", "MT118296.1" = "#233E4A",
  "OQ646790.1" = "#269B82", "OR222515.1" = "#D59B17"
)
flagship_theme <- theme_ggchord_publication() +
  theme(
    plot.background = element_rect(fill = flagship_paper, colour = NA),
    panel.background = element_rect(fill = flagship_paper, colour = NA),
    legend.background = element_rect(fill = flagship_paper, colour = NA),
    legend.key = element_rect(fill = flagship_paper, colour = NA),
    legend.position = "right", legend.box = "vertical", legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.margin = margin(12, 16, 12, 16)
  )
flagship_position <- position_strand(0.045)
fusion_data <- flagship_fusion_example
fusion_orientation <- fusion_data$orientation
fusion_display_hsp <- fusion_data$links
fusion_display_features <- fusion_data$features
fusion_relation_ab <- fusion_data$relation_ab
fusion_relation_ac <- fusion_data$relation_ac
flagship_selected <- fusion_data$selected
flagship_windows <- fusion_data$windows
flagship_focused <- fusion_data$focused
fusion_detail_features <- fusion_data$detail_features
fusion_strand_colours <- c("+" = "#CE1536", "-" = "#346B9D")
fusion_shared_styles <- list(
  scale_seq_colour_manual(values = flagship_colours, limits = names(flagship_colours), drop = FALSE, name = "Sequence"),
  scale_gene_fill_manual(values = fusion_strand_colours, limits = c("+", "-"), drop = FALSE, name = "Strand"),
  guides(seq_colour = guide_ggchord_legend(nrow = 2, position = "bottom"),
    gene_fill = guide_ggchord_legend(nrow = 1, position = "bottom", override.aes = list(strand = c("+", "-"))),
    ribbon_fill = ggplot2::guide_coloursteps(available_aes = "ribbon_fill", even.steps = FALSE,
      show.limits = TRUE, position = "bottom", direction = "horizontal",
      theme = theme(legend.key.width = grid::unit(85, "mm"), legend.key.height = grid::unit(3, "mm")))),
  flagship_theme + theme(legend.position = "bottom", legend.box = "horizontal",
    legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
    legend.margin = margin(2, 2, 2, 2), legend.spacing.x = grid::unit(2, "mm"),
    plot.margin = margin(3, 5, 3, 5), plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "#52646A"),
    plot.title.position = "plot", ggchord.axis.text = element_text(size = 9))
)
make_fusion_view <- function(coord, title, subtitle, label_radius = NULL) {
  label_layer <- if (is.null(label_radius)) {
    geom_seq_label(size = 3.5, seq_label_orientation = "horizontal")
  } else {
    geom_seq_label(size = 3.5, seq_label_radius = label_radius,
      seq_label_orientation = "horizontal")
  }
  ggchord(prepare_ggchord_data(flagship_sequences, fusion_display_hsp, feature_data = fusion_display_features)) +
    geom_link_ribbon(adjacent_only = FALSE, alpha = 0.25, linewidth = 0.10,
      nonadjacent_linewidth = 0.22, allow_repeated = TRUE) +
    layout_sequences(orientation = fusion_orientation) +
    geom_seq(aes(seq_colour = accver), linewidth = 1.1) +
    geom_feature(preset = 'gene', aes(gene_fill = strand), position = flagship_position, width = 0.04) +
    label_layer +
    coord + fusion_shared_styles + labs(title = title, subtitle = subtitle)
}
fusion_chord_view <- make_fusion_view(coord_chord(rotation = 45,
    xlim = c(-1.22, 1.02), ylim = c(-1.02, 1.22), fit = "manual"),
  "A  |  Chord", "6 HSPs; 3 CDS per genome", label_radius = 1.12) +
  theme(ggchord.axis.text = element_blank(), ggchord.axis.ticks = element_blank())
fusion_linear_view <- make_fusion_view(coord_collinear(ratio = 0.78),
  "B  |  Collinear", "The same intervals in linear coordinates") +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL)
fusion_radial_view <- make_fusion_view(coord_radial_synteny(start_angle = 45, inner_radius = 0.7),
  "C  |  Radial", "Shared intervals; A / D") +
  theme(ggchord.axis.text = element_blank(), ggchord.axis.ticks = element_blank())
fusion_zoom_view <- ggchord(prepare_ggchord_data(flagship_focused$seq_data, flagship_focused$ribbon_data,
  feature_data = fusion_detail_features)) +
  geom_link_ribbon(allow_repeated = TRUE, alpha = 0.25, linewidth = 0.10) +
  layout_sequences(orientation = 1) +
  geom_seq(aes(seq_colour = accver), linewidth = 1.1) +
  geom_feature(preset = 'gene', aes(gene_fill = strand), position = flagship_position, width = 0.04) +
  geom_gene_label_repel(size = 3.05, label_wrap = 22, label_max_lines = 2, position = flagship_position) +
  geom_seq_label(size = 3.5, seq_label_orientation = "horizontal") +
  coord_collinear(ratio = 0.85, fit = "manual",
    xlim = c(-0.42, 2.08), ylim = c(-0.72, 0.22)) + fusion_shared_styles +
  scale_seq_position_continuous(n.breaks = 3, minor_breaks = NULL, coordinates = "source") +
  labs(title = "E  |  Local zoom", subtitle = "Source windows + 800 bp flanks") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
fusion_context_sequences <- fusion_data$context_sequences
fusion_base_cds <- fusion_data$base_cds
fusion_base_start <- fusion_data$base_start
fusion_base_end <- fusion_data$base_end
fusion_base_cells <- fusion_data$base_cells
fusion_base_view <- ggplot(fusion_base_cells, aes(alignment_column, accver, label = symbol)) +
  geom_residue(overflow = "error", show.legend = FALSE) +
  scale_x_continuous(breaks = seq(1, 60, 10),
    labels = function(x) fusion_base_cells$source_position[match(x, fusion_base_cells$alignment_column)],
    expand = expansion(add = 0.6)) +
  scale_y_discrete(labels = setNames("DNA", fusion_base_cds$feature_id), expand = expansion(add = 0.7)) +
  labs(title = "F  |  Nucleotide detail", x = "Source position (bp)", y = NULL,
    subtitle = paste0(fusion_base_cds$accver, " / ", fusion_base_cds$feature_id,
      "  |  first 60 coding bases (", fusion_base_cds$strand, ")  |  ",
      fusion_base_start, "-", fusion_base_end, " bp")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.background = element_rect(fill = flagship_paper, colour = NA),
    axis.text = element_text(colour = "#25363B"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, colour = "#52646A", hjust = 0.5), plot.margin = margin(4, 12, 4, 12))
fusion_base_view <- register_ggchord_alignment(fusion_base_view, fusion_base_cells)
fusion_focus_regions <- fusion_data$focus_regions
fusion_context_features <- fusion_data$context_features
fusion_detail_view <- ggchord(prepare_ggchord_data(fusion_context_sequences, flagship_selected, feature_data = fusion_context_features)) +
  geom_link_ribbon(allow_repeated = TRUE, alpha = 0.25, linewidth = 0.10) +
  layout_sequences(orientation = 1) +
  geom_seq(aes(seq_colour = accver), linewidth = 1.1) +
  geom_feature(preset = 'gene', aes(gene_fill = strand), position = flagship_position, width = 0.04) +
  geom_seq_label(size = 3.5, seq_label_orientation = "horizontal") +
  coord_collinear(ratio = 0.7) + fusion_shared_styles +
  scale_seq_position_continuous(n.breaks = 2, minor_breaks = NULL) +
  labs(title = "D  |  Two-genome context", subtitle = "Full-length sequences; outlined windows are enlarged in E")
flagship_detail <- fusion_zoom_view
fusion_connections <- list(
  fusion_ribbon(fusion_data$relation_ab, aes(fill = relation),
    qpanel = "chord", spanel = "linear", alpha = 0.65, colour = "#3E5763", linewidth = 0.35,
    curvature = 0.09, allow_repeated = TRUE),
  fusion_ribbon(fusion_data$relation_ac, aes(fill = relation),
    qpanel = "chord", spanel = "radial", alpha = 0.65, colour = "#3E5763", linewidth = 0.35,
    curvature = 0.2, allow_repeated = TRUE),
  fusion_ribbon(fusion_data$relation_cd, aes(fill = relation),
    qpanel = "radial", spanel = "center", alpha = 0.65, colour = "#3E5763", linewidth = 0.35,
    curvature = 0.2, route = "avoid", route_width = 1.5, allow_repeated = TRUE)
)
make_flagship_fusion <- function(preview = FALSE) {
  output_styles <- list(
    guides(ribbon_fill = ggplot2::guide_coloursteps(available_aes = "ribbon_fill",
      even.steps = FALSE, show.limits = TRUE, position = "bottom", direction = "horizontal",
      theme = theme(legend.key.width = grid::unit(if (preview) 90 else 85, "mm"),
        legend.key.height = grid::unit(3, "mm")))),
    theme(plot.title = element_text(size = if (preview) 15 else 16),
      plot.subtitle = element_text(size = if (preview) 9 else 10),
      legend.text = element_text(size = if (preview) 8 else 9))
  )
  lift <- function(box) c(box[1:2], 0.20 + 0.80 * box[3:4])
  out <- magnify_ggchord(fusion_detail_view + output_styles, region = fusion_focus_regions,
    detail_plot = list(zoom = fusion_zoom_view + output_styles, bases = fusion_base_view), guides = "collect",
    center_box = lift(c(0.465, 0.97, 0.385, 0.59)),
    detail_box = list(zoom = lift(c(0.49, 0.99, 0.01, 0.36)), bases = c(0.025, 0.985, 0, 0.17)),
    colour = "#725094", alpha = 0.20, highlight_width = 2.2,
    linewidth = 0.40, connector = "envelope", route = "avoid") +
    fusion_panel(fusion_chord_view + output_styles, id = "chord",
      box = lift(c(0.01, 0.455, if (preview) 0.52 else 0.54, 0.99))) +
    fusion_panel(fusion_linear_view + output_styles, id = "linear", box = lift(c(0.52, 0.995, 0.60, 0.99))) +
    fusion_panel(fusion_radial_view + output_styles, id = "radial",
      box = lift(c(0.01, 0.455, 0, if (preview) 0.50 else 0.52)))
  for (connection in fusion_connections) out <- out + connection
  out <- out + scale_fill_manual(name = if (preview) "Links" else "Correspondence", values = c("A - B" = "#527FB8",
    "A - C" = "#CE9741", "C - D" = "#9064AC")) +
    guides(fill = guide_legend(nrow = 1, position = "bottom",
      theme = if(preview) theme(legend.key.width = grid::unit(3, "mm"),
        legend.key.height = grid::unit(3, "mm")) else NULL)) +
    theme(plot.background = element_rect(fill = flagship_paper, colour = NA),
      legend.position = "bottom", legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8), plot.title = element_text(size = 22, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11, colour = "#52646A"),
      plot.caption = element_text(size = 10, colour = "#52646A", hjust = 0),
      plot.margin = margin(8, 10, 8, 10))
  if (!preview) out <- out + labs(title = "One dataset, six linked views",
    subtitle = "Six selected HSPs and twelve source CDS per overview; full inputs retained separately",
    caption = "Coloured ribbons: observed HSPs. Source-window leaders connect D to E; the selected CDS bases connect E to F.")
  out
}
flag_fusion <- make_flagship_fusion()
flag_fusion
