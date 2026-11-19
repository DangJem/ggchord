library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)

ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label(seq_label_radius = c("MT108731.1" = 0.72, "MT118296.1" = 1)) +
  coord_chord()
