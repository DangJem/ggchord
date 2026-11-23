library(ggchord)
library(ggplot2)

plot <- ggchord(seq_data_example) +
  geom_seq(
    aes(colour = accver), show.legend = FALSE,
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_colour_manual(values = c(
    "MT108731.1" = "#0072B2", "MT118296.1" = "#D55E00",
    "OQ646790.1" = "#009E73", "OR222515.1" = "#CC79A7"
  )) +
  geom_seq_label(
    label_layout = layout_label_fixed(orientation = "feature", side = "negative", offset = 8)
  ) +
  layout_sequences(angular_gap = 9) +
  coord_chord() +
  theme_ggchord_minimal()

output_file <- file.path(tempdir(), "sequence-backbones.pdf")
ggsave(output_file, plot = plot, width = 160, height = 110, units = "mm")
plot
