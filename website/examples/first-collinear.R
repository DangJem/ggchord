library(ggchord)
refs <- seq_data_example[1:2, ]
hits <- subset(ribbon_data_example, qaccver == refs$accver[1] & saccver == refs$accver[2])
ggchord(refs) +
  geom_link_ribbon(data = hits) +
  geom_seq() +
  geom_seq_label() +
  coord_collinear()
