library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])

ggchord(reference) +
  geom_seq(
    ggplot2::aes(colour = accver), show.legend = FALSE,
    arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm"))
  ) +
  ggplot2::scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 4)
  ) +
  geom_link_ribbon(ggplot2::aes(fill = pident), data = hits) +
  layout_sequences(row_gap = 35) +
  scale_seq_position(guide = guide_seq_axis(side = "negative")) +
  coord_collinear() +
  theme_ggchord_minimal()
