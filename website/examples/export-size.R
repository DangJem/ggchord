library(ggchord)

plot <- ggchord(seq_data_example) +
  geom_seq(show.legend = FALSE) +
  geom_seq_label() +
  coord_chord()

output_file <- file.path(tempdir(), "sequence-backbones.pdf")
ggsave(output_file, plot = plot, width = 160, height = 110, units = "mm")
plot
