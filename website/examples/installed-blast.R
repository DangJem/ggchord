library(ggchord)
library(dplyr)

path <- system.file("extdata", "tutorial", "MT108731.1__MT118296.1.tsv",
                    package = "ggchord")
hits <- read_blast(path)
hits |>
  select(qaccver, saccver, qstart, qend, sstart, send, pident) |>
  head(4)
