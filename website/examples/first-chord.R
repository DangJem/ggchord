library(ggchord)
hit <- ribbon_data_example[1, ]
ggchord(seq_data_example) +
  geom_seq() +
  geom_seq_label() +
  geom_link_ribbon(data = hit) +
  coord_chord()
