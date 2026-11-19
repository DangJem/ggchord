library(ggchord)

record <- flagship_circular_example$references[["pUC19"]]

ggchord(record$seq_data) +
  geom_seq(show.legend = FALSE) +
  geom_feature(data = record$features, preset = "plasmid") +
  coord_circular()
