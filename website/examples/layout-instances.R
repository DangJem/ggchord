library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])

ggchord(reference) +
  layout_sequences(
    order = c(top = "MT108731.1", middle = "MT118296.1", bottom = "MT108731.1")
  ) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  geom_link_ribbon(data = slice_head(hits, n = 1), between = c("top", "middle")) +
  geom_link_ribbon(data = slice_tail(hits, n = 2), between = c("bottom", "middle")) +
  coord_collinear()
