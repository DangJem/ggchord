library(ggchord)
library(dplyr)

path <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
windows <- read_bedgraph(path)
windows |>
  select(accver, chromStart, chromEnd, start, end, value) |>
  head(5)
