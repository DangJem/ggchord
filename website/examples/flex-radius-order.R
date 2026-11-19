library(ggchord)

display_order <- c("OQ646790.1", "MT108731.1", "MT118296.1", "OR222515.1")

by_name <- ggchord(seq_data_example) +
  layout_sequences(order = display_order) +
  geom_seq(seq_radius = c("MT118296.1" = 0.72), show.legend = FALSE) +
  geom_seq_label() +
  coord_chord()

by_order <- ggchord(seq_data_example) +
  layout_sequences(order = display_order) +
  geom_seq(seq_radius = c(1, 0.72, 1, 1), show.legend = FALSE) +
  geom_seq_label() +
  coord_chord()

ggchord_fusion(by_name, center_box = c(0, 0.5, 0, 1)) +
  fusion_panel(by_order, id = "ordered", box = c(0.5, 1, 0, 1))
