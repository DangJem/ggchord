library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)

ggchord(reference) +
  geom_seq(
    ggplot2::aes(colour = accver), show.legend = FALSE,
    arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm"))
  ) +
  ggplot2::scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    label_layout = label_layout_fixed(orientation = "feature", side = "negative", offset = 8)
  ) +
  layout_sequences(angular_gap = 9) +
  coord_chord() +
  theme_ggchord_minimal()
