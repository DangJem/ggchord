library(ggchord)
library(dplyr)

path <- system.file("extdata", "plant-feature-groups.txt", package = "ggchord")
groups <- read_feature_relations(path, format = "named_groups")
groups$membership |>
  filter(group_id == first(group_id)) |>
  select(group_id, feature_id, .source_line)
