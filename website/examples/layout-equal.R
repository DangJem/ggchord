library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hit <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2]) |>
  slice_head(n = 1)

ggchord(reference) +
  layout_sequences(length_mode = "equal", row_gap = 35) +
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
  geom_link_ribbon(aes(fill = pident), data = hit) +
  scale_seq_position(
    breaks = list("MT108731.1" = c(1, 20000, 40000, 64323),
                  "MT118296.1" = c(1, 10000, 20000, 32090)),
    guide = guide_seq_axis(side = "negative")
  ) +
  coord_collinear() +
  theme_ggchord_minimal()
