library(ggchord)
library(ggplot2)
library(dplyr)

cds <- gene_data_example |>
  filter(
    (accver == "OR222515.1" & start %in% c(2624, 19246)) |
      (accver == "OQ646790.1" & start == 43765)
  ) |>
  mutate(strand = source_strand)

offsets <- list(
  "OR222515.1" = c("+" = 5, "-" = -9),
  "OQ646790.1" = 4
)
strand_position <- position_strand(offsets)

ggchord(seq_data_example) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00",
    "OQ646790.1" = "#009E73", "OR222515.1" = "#CC79A7"
  )) +
  geom_gene(
    data = filter(cds, strand == "+"), position = strand_position,
    thickness = 2, fill = "#CE1536"
  ) +
  geom_gene(
    data = filter(cds, strand == "-"), position = strand_position,
    thickness = 5, fill = "#346B9D"
  ) +
  geom_seq_label(
    label_layout = layout_label_fixed(orientation = "feature", side = "negative", offset = 8)
  ) +
  layout_sequences(angular_gap = 9) +
  coord_chord() +
  theme_ggchord_minimal()
