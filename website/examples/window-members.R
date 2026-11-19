library(ggchord)
library(dplyr)

reference <- seq_data_example |>
  filter(accver == "MT108731.1")
intervals <- ribbon_data_example |>
  filter(qaccver == reference$accver) |>
  transmute(accver = qaccver, subject = saccver,
            start = pmin(qstart, qend), end = pmax(qstart, qend))

plot <- ggchord(reference) +
  geom_seq() +
  stat_coverage(data = intervals, binwidth = 10000, id = "windows") +
  ggplot2::facet_wrap(ggplot2::vars(subject)) +
  coord_collinear()

export_ggchord_layout(plot, include = "track")$track_statistics$windows |>
  select(subject, window_start, window_end, covered, denominator, coverage, source_rows) |>
  filter(subject == "MT118296.1")
