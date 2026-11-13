# Reuse one measured layout: change drawing order only, without rebuilding it.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional='tools/long-term-audit/check-leader-draw-order.R')
pkgload::load_all(quiet=TRUE)
x <- readRDS('/tmp/ggchord-association-after-pCAMBIA1300.rds')
built <- x$built
stopifnot(!is.null(built))
data_before <- built$data
owners_before <- vapply(built$plot$layers,function(l)l$ggchord_layer_id %||% '',character(1))
clone_layer <- function(parent) {
 force(parent)
 ggplot2::ggproto(NULL,parent)
}
files <- list()
for(order in c('below','above')) {
 trial <- built
 for(i in seq_along(trial$plot$layers)) {
  layer <- trial$plot$layers[[i]]
  if (!isTRUE(layer$ggchord_params$is_feature_label)) next
  copy <- clone_layer(layer)
  copy$ggchord_params <- layer$ggchord_params
  copy$ggchord_params$label_segment_layer <- order
  trial$plot$layers[[i]] <- copy
 }
 path <- paste0('tools/long-term-audit/v0.27.0-leaders-',order,'.png')
 ggchord_diagnostic_device(path,'png',10,8,110)
 table <- ggchord_base_gtable(trial)
 grid::grid.newpage();grid::grid.draw(table); grDevices::dev.off()
 stopifnot(identical(trial$data,data_before),
  identical(vapply(trial$plot$layers,function(l)l$ggchord_layer_id %||% '',character(1)),owners_before))
 files[[order]] <- path
}
stopifnot(identical(built$data,data_before))
a <- png::readPNG(files$below)
b <- png::readPNG(files$above)
stopifnot(identical(dim(a),dim(b)))
different_pixels <- sum(apply(abs(a-b),c(1,2),max)>0)
stopifnot(different_pixels>0L)
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(files=files,trained_data_and_owner_ids_preserved=TRUE,
 different_pixels=different_pixels,same_dimensions=TRUE,
 layout='Reused final measured pCAMBIA1300, same 10x8in/110dpi device',
 runtime_sha256=as.list(audit$engine),input_sha256=as.list(audit$inputs)),
 'tools/long-term-audit/v0.27.0-leader-draw-order-results.json',auto_unbox=TRUE,pretty=TRUE)
cat('Below/above draw-order controls completed\n')
