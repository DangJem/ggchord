library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])

ggchord(reference) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 4)
  ) +
  geom_link_ribbon(data = hits, id = "observed", fill = "#3A7984", alpha = 0.65) +
  layout_sequences(row_gap = 35) +
  scale_seq_position(guide = guide_seq_axis(side = "negative")) +
  coord_collinear() +
  theme_ggchord_minimal()
