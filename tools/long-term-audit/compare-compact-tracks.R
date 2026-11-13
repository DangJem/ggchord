# Inspect the measured saved layouts, without rebuilding or changing inputs.
ids <- names(jsonlite::read_json('tools/long-term-audit/v0.27.0-compact-final-results.json')$cases)
stopifnot(length(ids)==13L)
observe <- function(prefix, id) {
 l <- readRDS(paste0('/tmp/',prefix,'-',id,'.rds'))
 labels <- l$gene_labels
 polygon <- l$gene_polys
 data.frame(case=id,stage=prefix,
  feature_label_count=nrow(labels),
  inner_source_count=sum(labels$lane>0L),
  inner_source_exterior=sum(labels$lane>0L & labels$feature_label_mode=='external'),
  total_feature_exterior=sum(labels$feature_label_mode=='external'),
  feature_min_radius=min(sqrt(polygon$x^2+polygon$y^2)),
  label_min_radius=min(sqrt(labels$text_x^2+labels$text_y^2)))
}
rows<-do.call(rbind,lapply(ids,function(id) do.call(rbind,lapply(
 c('v0.27.0-cluster-final','v0.27.0-compact-final'),observe,id=id))))
write.csv(rows,'tools/long-term-audit/v0.27.0-compact-tracks.csv',row.names=FALSE)
print(rows)
