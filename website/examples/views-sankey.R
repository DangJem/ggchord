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
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    anchor = "start",
    label_layout = layout_label_fixed(orientation = "vertical", along = -5)
  ) +
  geom_link_ribbon(aes(fill = pident), data = hit) +
  scale_fill_pident() +
  layout_sequences(column_gap = 75, row_gap = 11) +
  coord_sankey(fit = "labels") +
  theme_ggchord_minimal()
