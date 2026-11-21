library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hits <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2])

windows <- tibble(
  instance_id = c("top", "middle", "bottom"),
  accver = c(ids[1], ids[2], ids[1])
) |>
  left_join(select(reference, accver, length), by = "accver") |>
  transmute(instance_id, accver, start = 1L, end = length)

upper_hit <- hits |>
  slice_head(n = 1) |>
  mutate(qinstance = "top", sinstance = "middle")
lower_hits <- hits |>
  slice_tail(n = 2) |>
  mutate(qinstance = "bottom", sinstance = "middle")

ggchord(reference) +
  layout_sequences(order = windows$instance_id, windows = windows, row_gap = 29) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    label_layout = label_layout_fixed(orientation = "horizontal", offset = 4)
  ) +
  geom_link_ribbon(aes(fill = pident), data = upper_hit) +
  geom_link_ribbon(aes(fill = pident), data = lower_hits) +
  scale_seq_position(guide = guide_seq_axis(side = "negative")) +
  coord_collinear() +
  theme_ggchord_minimal()
