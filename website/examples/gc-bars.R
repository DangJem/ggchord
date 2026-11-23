library(ggchord)
library(ggplot2)
library(dplyr)

path <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
windows <- read_bedgraph(path) |>
  mutate(midpoint = (start + end) / 2, window_length = end - start + 1)
reference <- plasmid_example_pUC19 |>
  select(accver, length)

ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  geom_track_col(
    aes(x = midpoint, y = value, width = window_length),
    data = windows, bounds = "clip", fill = "#99C9C2"
  ) +
  layout_track(thickness = 15) +
  scale_y_continuous(name = "GC (%)", limits = c(0, 100)) +
  coord_collinear()
