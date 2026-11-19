library(ggchord)
library(dplyr)

record <- flagship_circular_example$references[["pUC19"]]
features <- record$features |>
  transmute(accver, start, end, strand, type, anno, feature_id = match_id)
captions <- features |>
  filter(type == "CDS")

ggchord(record$seq_data) +
  geom_seq(show.legend = FALSE) +
  geom_feature(
    aes(fill = type), data = features, id = "features",
    position = position_feature_stack(), show.legend = FALSE
  ) +
  geom_feature_label(
    aes(label = anno), data = captions, anchor_layer = "features",
    placement = "repel", label_layout = "callout", size = 4
  ) +
  coord_circular()
