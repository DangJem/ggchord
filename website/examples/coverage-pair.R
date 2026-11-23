library(ggchord)
library(dplyr)

reference <- seq_data_example |>
  filter(accver == "MT108731.1")
intervals <- ribbon_data_example |>
  filter(qaccver == reference$accver, saccver == "MT118296.1") |>
  transmute(accver = qaccver,
            start = pmin(qstart, qend), end = pmax(qstart, qend))

ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  layout_track("coverage", thickness = 12) +
  geom_track_col(
    data = intervals, stat = "coverage", binwidth = 10000,
    track = "coverage", fill = "#327B83"
  ) +
  scale_y_continuous(limits = c(0, 1), name = "Coverage") +
  coord_collinear()
