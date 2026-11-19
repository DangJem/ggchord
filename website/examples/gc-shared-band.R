library(ggchord)
library(dplyr)

path <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
windows <- read_bedgraph(path) |>
  mutate(midpoint = (start + end) / 2)
reference <- plasmid_example_pUC19 |>
  select(accver, length)

ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  track_band("gc", width = grid::unit(15, "mm")) +
  geom_track_col(aes(y = value), data = windows, track = "gc", fill = "#99C9C2") +
  geom_track_line(
    aes(x = midpoint, y = value), data = windows, track = "gc", colour = "#1F5961"
  ) +
  scale_track_value_continuous("gc", name = "GC (%)", limits = c(0, 100)) +
  coord_collinear()
