library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)

ggchord(reference) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    label_layout = label_layout_fixed(
      orientation = "feature", side = "negative",
      offset = c("MT108731.1" = 18)
    )
  ) +
  layout_sequences(angular_gap = 9) +
  coord_chord() +
  theme_ggchord_minimal()
