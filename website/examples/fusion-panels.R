library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)

overview <- ggchord(filter(reference, accver == ids[1])) +
  geom_seq(geometry = "double", thickness = 1.3) +
  coord_circular(start_angle = 90) +
  theme_ggchord_minimal()

detail <- ggchord(reference) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    anchor = "end", label_layout = label_layout_fixed(orientation = "horizontal", offset = 3)
  ) +
  layout_sequences(row_gap = 24) +
  coord_collinear(axis = "y", fit = "labels") +
  theme_ggchord_minimal()

ggchord_fusion(list(overview = overview, detail = detail))
