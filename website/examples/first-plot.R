library(ggchord)
ggchord(seq_data_example) +
  geom_seq() +
  geom_link_ribbon(data = head(ribbon_data_example, 4)) +
  coord_chord()
