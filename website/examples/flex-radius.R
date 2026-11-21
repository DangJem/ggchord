library(ggchord)
library(ggplot2)

ggchord(seq_data_example) +
  layout_sequences(radius = c("MT118296.1" = 0.72), angular_gap = 9) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00",
    "OQ646790.1" = "#009E73", "OR222515.1" = "#CC79A7"
  )) +
  geom_seq_label(
    label_layout = label_layout_fixed(orientation = "feature", side = "negative", offset = 8)
  ) +
  coord_chord() +
  theme_ggchord_minimal()
