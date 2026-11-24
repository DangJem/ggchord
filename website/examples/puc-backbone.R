library(ggchord)

record <- flagship_circular_example$references[["pUC19"]]

ggchord(record$seq_data) +
  geom_seq(geometry = "double", thickness = 1.3, show.legend = FALSE) +
  coord_circular() +
  theme_ggchord_minimal()
