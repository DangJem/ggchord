library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])

plot <- ggchord(reference) +
  geom_seq() +
  geom_link_ribbon(data = hits, id = "observed") +
  coord_collinear()
layout <- export_ggchord_layout(plot, data_stage = "source")
layout$layer_data$observed$source |>
  select(qaccver, saccver, qstart, qend, sstart, send)
