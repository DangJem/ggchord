library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hit <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2]) |>
  slice_head(n = 1)

overview <- ggchord(filter(reference, accver == ids[1])) +
  geom_seq(show.legend = FALSE) +
  coord_circular()

detail <- ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  coord_collinear(seq_axis = "vertical")

ggchord_fusion(overview) +
  fusion_panel(detail, id = "detail", side = "right") +
  fusion_ribbon(
    hit, qpanel = "center", spanel = "detail",
    fill = "#269B82", alpha = 0.55
  )
