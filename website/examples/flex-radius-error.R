library(ggchord)

invalid_radius <- ggchord(seq_data_example) +
  layout_sequences(radius = c("+" = 0.72, "-" = 1)) +
  geom_seq() +
  coord_chord()

tryCatch(
  ggplot2::ggplot_build(invalid_radius),
  error = function(error) conditionMessage(error)
)
