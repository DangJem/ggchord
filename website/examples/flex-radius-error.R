library(ggchord)

invalid_radius <- ggchord(seq_data_example) +
  geom_seq(seq_radius = c("+" = 0.72, "-" = 1)) +
  coord_chord()

tryCatch(
  get_chord_layout(invalid_radius),
  error = function(error) conditionMessage(error)
)
