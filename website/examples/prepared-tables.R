library(ggchord)
library(dplyr)

reference <- seq_data_example |>
  select(accver, length)
features <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand, type = "CDS", anno)
prepared <- prepare_ggchord_data(reference, feature_data = features)
prepared$feature_data |>
  select(accver, start, end, strand, type, anno) |>
  head(4)
