library(ggchord)
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
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  geom_link_ribbon(data = hits, adjacent_only = FALSE) +
  coord_radial_synteny(start_angle = 90)
