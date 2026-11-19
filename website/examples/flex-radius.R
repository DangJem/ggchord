library(ggchord)

ggchord(seq_data_example) +
  geom_seq(seq_radius = c("MT118296.1" = 0.72)) +
  geom_seq_label() +
  coord_chord()
