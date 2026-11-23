library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1", "OQ646790.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver %in% ids[2:3]) |>
  group_by(saccver) |>
  slice_head(n = 1) |>
  ungroup()

ggchord(reference) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00",
    "OQ646790.1" = "#009E73"
  )) +
  geom_seq_label(
    anchor = "end", label_layout = label_layout_fixed(orientation = "feature", along = 9)
  ) +
  geom_link_ribbon(aes(fill = pident), data = hits) +
  scale_fill_pident() +
  coord_radial_synteny(start_angle = 90, fit = "labels") +
  theme_ggchord_minimal()
