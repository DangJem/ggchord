# Run from the repository root using the current package source.
# See data-raw/ for the source records and selection rules of the bundled fixtures.
library(ggchord)
library(ggplot2)
library(dplyr)

comparison <- flagship_comparison_example$chord
references <- as_ggchord_seq(comparison$seq_data)
relationships <- comparison$links
genes <- comparison$features |>
  mutate(label = anno, feature_id = paste(accver, start, end, sep = "/"))
comparison_plot <- ggchord(references) +
  geom_link_ribbon(aes(fill = pident), data = relationships, alpha = 0.58) +
  geom_seq(colour = "#29444D", linewidth = 0.65) +
  geom_gene(data = genes, fill = "#266E69", colour = NA, thickness = 1.4, id = "genes") +
  geom_seq_label(size = 2.9, colour = "#29444D") +
  scale_seq_position(guide = "none") +
  theme_ggchord_publication()

chord <- comparison_plot +
  layout_sequences(angular_gap = 11) +
  coord_chord(start_angle = 90, fit = "labels") +
  theme(legend.position = "left") +
  labs(title = "Four related genomes")
linear <- comparison_plot +
  layout_sequences(orientation = comparison$orientation, row_gap = 10) +
  coord_collinear(fit = "labels") +
  theme(legend.position = "left") +
  labs(title = "Gene order and aligned intervals")
radial <- comparison_plot +
  layout_sequences(inner_radius = 0.26) +
  coord_radial_synteny(start_angle = 45, fit = "labels") +
  labs(title = "Relationships along radial reference paths")
sankey <- comparison_plot +
  layout_sequences(column_gap = 40, row_gap = 8) +
  coord_sankey(fit = "labels") +
  labs(title = "Query and subject intervals")

plasmid <- flagship_circular_example$dense
plasmid_ref <- as_ggchord_seq(plasmid$seq_data, topology = "circular")
# Represent the circular seam by two closed source segments.
annotation_segments <- plasmid$features$segments
wraps <- which(plasmid$features$start > plasmid$features$end)
for (i in wraps) {
  annotation_segments[[i]] <- data.frame(
    start = c(plasmid$features$start[i], 1),
    end = c(plasmid_ref$length, plasmid$features$end[i])
  )
}
annotations <- plasmid$features |>
  mutate(
    feature_id = as.character(match_id), label = anno,
    segments = annotation_segments,
    start = if_else(cross_origin, 1, start),
    end = if_else(cross_origin, plasmid_ref$length, end),
    class = if_else(type %in% c("CDS", "rep_origin"), type, "Other feature")
  )
labelled <- annotations |>
  filter(type %in% c("CDS", "rep_origin"))
cuts <- plasmid$sites |>
  filter(enzyme_site_count == 1) |>
  transmute(accver, pos = start) |>
  distinct()
circular <- ggchord(plasmid_ref) +
  geom_seq(colour = "#29444D", linewidth = 0.65) +
  geom_feature(
    aes(fill = class), data = annotations, thickness = 2.3,
    position = position_strand(offset = 2.2), id = "features"
  ) +
  geom_feature_label_repel(
    aes(label = label), data = labelled, anchor_layer = "features", size = 2.4
  ) +
  geom_site(data = cuts, size = 1, colour = "#C27644", position = position_offset(-4)) +
  geom_seq_center_label(aes(label = accver), size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c(CDS = "#347880", rep_origin = "#C27644",
    "Other feature" = "#A8C6BE"), name = "Annotation") +
  coord_circular(start_angle = 90, fit = "labels") +
  theme_ggchord_plasmid() +
  labs(title = "pETDuet-1 annotations and cut-site positions")

hsp <- flagship_dotplot_example
dotplot <- ggchord(hsp$seq_data) +
  geom_alignment_segment(
    aes(colour = pident, linewidth = length), data = hsp$links, alpha = 0.9
  ) +
  layout_dotplot(x_gap = hsp$axis_gap, y_gap = hsp$axis_gap) +
  coord_dotplot() +
  scale_linewidth_continuous(range = c(0.45, 1.15), name = "Alignment length") +
  theme(legend.position = "bottom") +
  labs(title = "High-scoring segments on query and subject axes")

focus <- flagship_fusion_example$focus_regions |>
  slice_head(n = 1) |>
  select(accver, start, end)
focus_genes <- flagship_fusion_example$all_features |>
  filter(accver == focus$accver, type == "CDS", end >= focus$start, start <= focus$end) |>
  mutate(label = anno)
focus_ref <- references |>
  filter(accver == focus$accver)
detail <- ggchord(focus_ref) +
  geom_seq(colour = "#29444D", linewidth = 0.7) +
  geom_gene(data = focus_genes, fill = "#266E69", colour = NA, id = "genes") +
  geom_gene_label_repel(aes(label = label), anchor_layer = "genes", size = 2.6) +
  layout_sequences(windows = mutate(focus, instance_id = "detail")) +
  coord_collinear(fit = "labels") +
  theme_ggchord_publication() +
  labs(title = "Source interval detail")
fusion <- ggchord_fusion(
  list(overview = chord, detail = detail),
  layout = layout_fusion_grid(matrix(c("overview", "detail"), nrow = 1),
    widths = c(1.15, 1))
) +
  fusion_zoom("overview", "detail", focus)

figures <- list(
  showcasechord = chord, showcasecollinear = linear, flagradial = radial,
  flagsankey = sankey, showcaseplasmid = circular, flagdotplot = dotplot,
  flagfusion = fusion
)
sizes <- list(
  showcasechord = c(7.2, 6.1), showcasecollinear = c(7.5, 4.8),
  flagradial = c(7.2, 6.1), flagsankey = c(7.5, 6),
  showcaseplasmid = c(7.2, 6.1), flagdotplot = c(7.2, 6.1),
  flagfusion = c(11, 5.4)
)
for (name in names(figures)) {
  size <- sizes[[name]]
  ggsave(file.path("man", "figures", paste0(name, ".png")), figures[[name]],
    width = size[1], height = size[2], dpi = 144, bg = "white")
}
