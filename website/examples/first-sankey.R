library(ggchord)
hits <- ribbon_data_example[1:3, ]
ggchord(seq_data_example[1:2, ]) +
  geom_link_ribbon(data = hits) +
  geom_seq() +
  geom_seq_label() +
  coord_sankey()
