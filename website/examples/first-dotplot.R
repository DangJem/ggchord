library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hit <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2]) |>
  slice_head(n = 1)

ggchord(reference) +
  geom_alignment_segment(data = hit, linewidth = 1) +
  coord_dotplot(x_seqs = ids[1], y_seqs = ids[2]) +
  ggplot2::theme(plot.margin = ggplot2::margin(8, 22, 8, 8))
