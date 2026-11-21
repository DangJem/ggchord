library(ggchord)
library(dplyr)

record <- flagship_circular_example$references[["pUC19"]]
features <- record$features |>
  filter(type %in% c("CDS", "primer_bind", "misc_feature")) |>
  select(accver, start, end, strand, type)

ggchord(record$seq_data) +
  geom_seq(geometry = "double", thickness = 1.3, show.legend = FALSE) +
  geom_feature(
    ggplot2::aes(fill = type), data = features,
    position = position_feature_stack(base_position = position_offset(-3), direction = "negative")
  ) +
  ggplot2::scale_fill_manual(
    name = "Feature", values = c(CDS = "#2B7580", primer_bind = "#C5753B",
                                 misc_feature = "#81729A")
  ) +
  coord_circular(start_angle = 90) +
  theme_ggchord_minimal()
