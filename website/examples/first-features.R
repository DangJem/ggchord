library(ggchord)
library(ggplot2)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
hit <- ribbon_data_example |>
  filter(qaccver == ids[1], saccver == ids[2]) |>
  slice_head(n = 1)
features <- gene_data_example |>
  filter(accver %in% ids) |>
  transmute(accver, start, end, strand = source_strand)

ggchord(reference) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00"
  )) +
  geom_seq_label(
    label_layout = layout_label_fixed(orientation = "feature", side = "negative", offset = 8)
  ) +
  geom_link_ribbon(aes(fill = pident), data = hit) +
  scale_fill_pident() +
  geom_gene(data = features, fill = "#29766E", show.legend = FALSE) +
  layout_sequences(angular_gap = 9) +
  coord_chord() +
  theme_ggchord_minimal()
