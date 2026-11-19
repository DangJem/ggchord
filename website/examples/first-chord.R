library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hit <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2]) |>
  slice_head(n = 1)

ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  geom_link_ribbon(data = hit) +
  coord_chord()
