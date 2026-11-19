library(ggchord)
library(dplyr)

cds <- gene_data_example |>
  filter(
    (accver == "OR222515.1" & start %in% c(2624, 19246)) |
      (accver == "OQ646790.1" & start == 43765)
  ) |>
  mutate(strand = source_strand)

offsets <- list(
  "OR222515.1" = c("+" = 0.08, "-" = -0.14),
  "OQ646790.1" = 0.06
)

ggchord(seq_data_example) +
  geom_seq() +
  geom_feature(
    data = cds, preset = "gene",
    position = position_strand(offsets),
    width = c("+" = 0.03, "-" = 0.08)
  ) +
  geom_seq_label() +
  coord_chord()
