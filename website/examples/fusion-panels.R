library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)

overview <- ggchord(filter(reference, accver == ids[1])) +
  geom_seq(show.legend = FALSE) +
  coord_circular()

detail <- ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  coord_collinear(seq_axis = "vertical")

ggchord_fusion(overview) +
  fusion_panel(detail, id = "detail", side = "right")
