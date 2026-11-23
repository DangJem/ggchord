library(ggchord)
library(ggplot2)

display_order <- c("OQ646790.1", "MT108731.1", "MT118296.1", "OR222515.1")

backbones <- ggchord(seq_data_example) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00",
    "OQ646790.1" = "#009E73", "OR222515.1" = "#CC79A7"
  )) +
  geom_seq_label(
    label_layout = layout_label_fixed(orientation = "feature", side = "negative", offset = 8)
  ) +
  coord_chord() +
  theme_ggchord_minimal()

by_name <- backbones +
  layout_sequences(
    order = display_order, radius = c("MT118296.1" = 0.72), angular_gap = 9
  )
by_order <- backbones +
  layout_sequences(order = display_order, radius = c(1, 0.72, 1, 1), angular_gap = 9)

ggchord_fusion(list(by_name = by_name, by_order = by_order))
