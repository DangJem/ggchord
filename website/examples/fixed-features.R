library(ggchord)
library(dplyr)

record <- flagship_circular_example$references[["pUC19"]]
features <- record$features |>
  filter(type %in% c("CDS", "primer_bind", "misc_feature")) |>
  select(accver, start, end, strand, type)

ggchord(record$seq_data) +
  geom_seq(show.legend = FALSE) +
  geom_feature(data = features, fill = "#2B7580", position = position_feature_stack()) +
  coord_circular()
