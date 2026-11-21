# Synthetic intervals for a native ggplot2-style collinearity workflow.
library(ggchord)
library(ggplot2)
library(dplyr)

workflow_sequences <- tibble(
  accver = c("A", "B", "C"), length = c(1000, 800, 1200)
)
workflow_links <- tibble(
  qaccver = c("A", "A", "B"), saccver = c("B", "C", "C"),
  qstart = c(100, 100, 200), qend = c(300, 300, 400),
  sstart = c(400, 100, 500), send = c(200, 300, 300)
)
workflow_features <- tibble(
  feature_id = c("A_target", "B_target", "C_target"),
  accver = c("A", "B", "C"), start = c(100, 200, 300),
  end = c(300, 400, 500), strand = c("+", "-", "+"),
  type = "CDS", label = c("target A", "target B", "target C")
)

# Compare adjacent rows; all three source comparisons remain in workflow_links.
adjacent_links <- workflow_links |>
  filter(!(qaccver == "A" & saccver == "C"))
anchors <- workflow_features |>
  select(accver, start, end)

workflow_base <- ggchord(workflow_sequences) +
  geom_link_ribbon(data = adjacent_links, alpha = 0.25) +
  geom_seq() +
  geom_seq_label() +
  geom_gene(data = workflow_features, id = "genes") +
  geom_gene_label_repel(aes(label = label), anchor_layer = "genes") +
  coord_collinear()

workflow_aligned <- workflow_base +
  layout_sequences(anchor = anchors, align = "center")

# Draw windows against the original reference and label their source coordinates.
windows <- tibble(
  instance_id = c("A_window", "B_window", "C_window"),
  accver = c("A", "B", "C"),
  start = c(50, 150, 250), end = c(600, 650, 900)
)
workflow_source_axis <- workflow_base +
  layout_sequences(windows = windows) +
  geom_seq_break() +
  scale_seq_position(coordinates = "source")

# Place the mark at a source position; query its display location at export size.
mark <- tibble(accver = "A", pos = 400)
workflow_source_axis <- workflow_source_axis +
  geom_site(data = mark, shape = "point", colour = "#CC6677", size = 2)
workflow_mark <- map_ggchord_position(
  workflow_source_axis, transmute(mark, source_accver = accver, source_position = pos),
  width = 9, height = 6, device = "png"
)

if (interactive()) print(workflow_aligned)
# ggsave("focused-clusters.png", workflow_source_axis, width = 9, height = 6)
