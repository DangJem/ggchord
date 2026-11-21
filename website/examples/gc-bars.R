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
  track_band("gc", thickness = 15) +
  geom_track_col(
    aes(x = midpoint, y = value, width = window_length),
    data = windows, track = "gc", bounds = "clip", fill = "#99C9C2"
  ) +
  scale_track_value("gc", name = "GC (%)", limits = c(0, 100)) +
  coord_collinear()
