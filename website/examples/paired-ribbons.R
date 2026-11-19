library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])

ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  geom_link_ribbon(data = hits, id = "observed", fill = "#3A7984", alpha = 0.65) +
  coord_collinear()
