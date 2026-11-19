library(ggchord)
library(dplyr)

features <- flagship_circular_example$references[["pUC19"]]$features
coding <- features |>
  filter(type == "CDS", anno == "lacZα") |>
  transmute(accver, start, end, strand, feature_id = "lacZ_alpha", phase = 0L)
model <- build_ggchord_cds(coding, genetic_code = 11)
residues <- tibble(feature_id = "display_span", parent_id = "lacZ_alpha",
                   start = 2L, end = 5L)

project_ggchord_features(residues, model) |>
  select(accver, start, end, strand)
