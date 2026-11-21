library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])

ggchord(reference) +
  geom_alignment_segment(data = hits, colour = "#315F73", linewidth = 1.2, alpha = 1) +
  layout_dotplot(x_order = ids[1], y_order = ids[2]) +
  coord_dotplot() +
  theme(plot.margin = margin(8, 22, 8, 8))
