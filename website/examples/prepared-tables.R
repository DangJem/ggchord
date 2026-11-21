library(ggchord)
library(dplyr)

ids <- c("MT108731.1", "MT118296.1")
reference <- seq_data_example |>
  filter(accver %in% ids)
features <- gene_data_example |>
  filter(accver %in% ids) |>
  transmute(accver, start, end, strand = source_strand, type = "CDS", anno)

prepared <- ggchord_data(reference, features = list(cds = features))
prepared$features$cds |>
  select(accver, start, end, strand, type, anno) |>
  slice_head(n = 4)
