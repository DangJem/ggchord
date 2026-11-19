library(ggchord)
observed <- flagship_comparison_example$chord
ggchord(observed$seq_data) +
  geom_link_ribbon(data = observed$links, adjacent_only = FALSE, alpha = 0.45) +
  geom_seq() +
  geom_seq_label() +
  coord_radial_synteny()
