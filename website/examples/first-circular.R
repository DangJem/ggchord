library(ggchord)
puc <- flagship_circular_example$references[["pUC19"]]
ggchord(puc$seq_data) +
  geom_seq() +
  geom_feature(data = puc$features, preset = "plasmid") +
  coord_circular()
