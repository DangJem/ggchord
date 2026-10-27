# Synthetic teaching data, not records of a named biological sequence.
# Run with the development version of ggchord loaded.
library(ggchord)

collinear_sequences <- data.frame(
  accver = c("cluster_A", "cluster_B", "cluster_C"),
  length = c(1000, 900, 1100)
)
collinear_genes <- data.frame(
  accver = rep(collinear_sequences$accver, each = 2),
  start = c(100, 350, 120, 310, 150, 510),
  end = c(280, 480, 260, 430, 320, 660),
  strand = rep(c("+", "-"), 3),
  anno = c("alpha", "beta", "alpha", "beta", "alpha", "gamma")
)
collinear_alignments <- data.frame(
  qaccver = c("cluster_A", "cluster_B", "cluster_A"),
  saccver = c("cluster_B", "cluster_C", "cluster_C"),
  qstart = c(100, 310, 500), qend = c(280, 430, 620),
  sstart = c(110, 600, 680), send = c(290, 480, 800),
  pident = c(94, 87, 81)
)
collinear_features <- data.frame(
  accver = "cluster_A", start = c(180, 220), end = c(330, 290),
  strand = "+", type = c("operator", "repeat"),
  label = c("operator", "repeat")
)
collinear_points <- data.frame(
  qaccver = "cluster_A", saccver = "cluster_B",
  qpos = c(200, 700), spos = c(200, 700)
)
collinear_lanes <- position_feature_stack(
  spacing = 0.09, base_position = position_strand(0.18)
)

collinear_example <- ggchord(
  collinear_sequences, collinear_alignments, collinear_genes,
  validate = "none"
) +
  geom_link_ribbon(alpha = 0.24) +
  geom_seq(
    seq_order = c("cluster_C", "cluster_A", "cluster_B"),
    seq_orientation = c(cluster_A = 1, cluster_B = -1, cluster_C = 1),
    seq_gap = 0.38
  ) +
  geom_seq_region(
    data = data.frame(accver = "cluster_C", start = 700, end = 900),
    region_offset = 0.10, region_side = "inside"
  ) +
  geom_gene(aes(gene_fill = anno), position = position_strand(0.08)) +
  geom_gene_label(
    position = position_strand(0.08),
    gene_label_orientation = "horizontal", gene_label_side = "auto",
    gene_label_overlap = "nudge"
  ) +
  geom_feature(
    aes(feature_fill = type, feature_shape = type),
    data = collinear_features, position = collinear_lanes,
    feature_width = 0.045
  ) +
  geom_feature_label_repel(
    data = collinear_features, position = collinear_lanes
  ) +
  geom_link_line(
    data = collinear_points, link_type = "straight", link_gap = 0.03,
    colour = "#626D78", alpha = 0.5,
    arrow = arrow(length = unit(1.5, "mm"))
  ) +
  geom_seq_label() +
  scale_gene_fill_manual(values = c(
    alpha = "#4477AA", beta = "#CC6677", gamma = "#228833"
  )) +
  scale_feature_fill_manual(values = c(
    operator = "#AA8833", `repeat` = "#66AABB"
  )) +
  scale_feature_shape_manual(values = c(
    operator = "block", `repeat` = "chevron"
  )) +
  scale_seq_position_continuous(
    breaks = c(0, 250, 500, 750, 1000),
    labels = function(x) paste0(x, " bp"), minor_breaks = NULL
  ) +
  coord_collinear() +
  labs(title = "Three synthetic gene clusters")

if (interactive()) print(collinear_example)

# Save using the standard ggplot2 interface:
# ggsave("collinear.png", collinear_example, width = 10, height = 7, dpi = 160)
# Inspect real linear positions and source rows:
# export_ggchord_layout(collinear_example, include = c("seq", "gene", "link"))
