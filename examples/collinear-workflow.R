# Synthetic intervals for a native ggplot2-style collinearity workflow.
library(ggchord)

workflow_sequences <- data.frame(
  accver = c("A", "B", "C"), length = c(1000, 800, 1200)
)
workflow_links <- data.frame(
  qaccver = c("A", "A", "B"), saccver = c("B", "C", "C"),
  qstart = c(100, 100, 200), qend = c(300, 300, 400),
  sstart = c(400, 100, 500), send = c(200, 300, 300)
)
workflow_features <- data.frame(
  feature_id = c("A_target", "B_target", "C_target"),
  accver = c("A", "B", "C"), start = c(100, 200, 300),
  end = c(300, 400, 500), strand = c("+", "-", "+"),
  type = "CDS", label = c("target A", "target B", "target C")
)

# Direction and anchor alignment are coordinate settings, composed with +.
# The non-adjacent A/C ribbon is omitted by default; input records stay intact.
workflow_aligned <- ggchord(workflow_sequences, workflow_links,
  validate = "none") +
  geom_link_ribbon(alpha = .25) + geom_seq() + geom_seq_label() +
  geom_feature(data = workflow_features,
    position = position_strand(.12)) +
  geom_feature_label_repel(data = workflow_features, label_side = "auto") +
  coord_collinear(anchor = workflow_features[c("accver", "start", "end")],
    align = "center", orientation = "auto")

# Explicitly include all comparisons, or select a named pair:
# geom_link_ribbon(adjacent_only = FALSE)
# geom_link_ribbon(pairs = data.frame(qaccver = "A", saccver = "C"))
# geom_link_ribbon(adjacent_only = FALSE, linetype = "solid")

# Focus retains source coordinates; end marks identify cropped sequence ends.
workflow_focused <- focus_ggchord_data(workflow_sequences,
  ribbon_data = workflow_links, gene_data = workflow_features,
  loci = data.frame(accver = c("A", "B", "C"),
    start = c(50, 150, 250), end = c(600, 650, 900)))
workflow_source_axis <- ggchord(workflow_focused$seq_data,
  workflow_focused$ribbon_data, workflow_focused$gene_data,
  validate = "none") +
  geom_link_ribbon(alpha = .25) + geom_seq() + geom_seq_break() +
  geom_seq_label() + geom_gene(position = position_strand(.12)) +
  scale_seq_position_continuous(coordinates = "source") +
  coord_collinear(orientation = "auto")

# Mapping helpers also support ordinary ggplot2 annotations.
workflow_mark <- map_ggchord_position(workflow_source_axis,
  "A", 400, source = TRUE)
workflow_source_axis <- workflow_source_axis +
  ggplot2::geom_point(data = workflow_mark, ggplot2::aes(x, y),
    inherit.aes = FALSE, colour = "#CC6677", size = 2)

if (interactive()) print(workflow_aligned)
# ggsave("focused-clusters.png", workflow_source_axis, width = 9, height = 6)
