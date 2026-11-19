library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])
midpoints <- hits |>
  transmute(qaccver, saccver,
            qpos = floor((qstart + qend) / 2),
            spos = floor((sstart + send) / 2))

ggchord(reference) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  geom_link_line(data = midpoints, colour = "#316774", linewidth = .85, alpha = 1) +
  coord_collinear()
