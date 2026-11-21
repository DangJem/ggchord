library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hit <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2]) |>
  slice_head(n = 1)

ggchord(reference) +
  geom_alignment_segment(data = hit, colour = "#315F73", linewidth = 1.2, alpha = 1) +
  layout_dotplot(x_order = ids[1], y_order = ids[2]) +
  coord_dotplot() +
  ggplot2::theme(plot.margin = ggplot2::margin(8, 22, 8, 8))
