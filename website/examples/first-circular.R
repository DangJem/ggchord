library(ggchord)
library(ggplot2)

record <- flagship_circular_example$references[["pUC19"]]

ggchord(record$seq_data) +
  geom_seq(geometry = "double", thickness = 1.3, show.legend = FALSE) +
  geom_feature(
    aes(shape = feature_shape, fill = feature_color),
    data = record$features, thickness = 3,
    colour = "#647273", linewidth = 0.22,
    position = position_feature_stack(
      spacing = 2, base_position = position_offset(-5), direction = "negative"
    ),
    show.legend = FALSE
  ) +
  scale_shape_identity() +
  scale_fill_identity() +
  coord_circular() +
  theme_ggchord_minimal()
