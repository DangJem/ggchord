library(ggchord)
library(ggplot2)
library(dplyr)

record <- flagship_circular_example$references[["pUC19"]]
features <- record$features |>
  transmute(accver, start, end, strand, type, anno, feature_id = match_id)
captions <- features |>
  filter(type == "CDS")

ggchord(record$seq_data) +
  geom_seq(geometry = "double", thickness = 1.3, show.legend = FALSE) +
  geom_feature(
    aes(fill = type), data = features, id = "features",
    position = position_feature_stack(base_position = position_offset(-3), direction = "negative"), show.legend = FALSE
  ) +
  geom_feature_label_repel(
    aes(label = anno), data = captions, anchor_layer = "features",
    label_layout = label_layout_repel(arrangement = "columns"), size = 4
  ) +
  coord_circular(start_angle = 90, fit = "labels") +
  theme_ggchord_minimal()
