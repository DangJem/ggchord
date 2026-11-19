library(ggchord)
library(dplyr)

path <- system.file("extdata", "tutorial", "MT108731.1.gff3", package = "ggchord")
features <- read_gff3(path, unstranded = "keep")
features |>
  select(accver, start, end, strand, type, anno, .source_line) |>
  head(5)
