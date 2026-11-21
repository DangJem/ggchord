library(ggchord)

record <- flagship_circular_example$references[["pUC19"]]

ggchord(record$seq_data) +
  geom_seq(geometry = "double", thickness = 1.3, show.legend = FALSE) +
  geom_feature(
    ggplot2::aes(shape = feature_shape, fill = feature_color),
    data = record$features, thickness = 3,
    colour = "#647273", linewidth = 0.22,
    position = position_feature_stack(
      spacing = 2, base_position = position_offset(-5), direction = "negative"
    ),
    show.legend = FALSE
  ) +
  ggplot2::scale_shape_identity() +
  ggplot2::scale_fill_identity() +
  coord_circular(start_angle = 90) +
  theme_ggchord_minimal()
