# Inspect completed measured layouts; no plot rebuild or release checks.
source('tools/long-term-audit/feature-association-metrics.R')
pkgload::load_all(quiet=TRUE)
observations <- jsonlite::read_json('tools/long-term-audit/v0.27.0-compact-final-results.json')
ids <- names(observations$cases)
stopifnot(length(ids)==13L)
rows <- list()
for (id in ids) {
 current <- readRDS(paste0('/tmp/v0.27.0-compact-final-',id,'.rds'))
 previous <- readRDS(paste0('/tmp/v0.27.0-cluster-final-',id,'.rds'))
 # Compare the actual input rows behind every feature, primer and site owner.
 for(owner in names(previous$layer_inputs)) for(component in names(previous$layer_inputs[[owner]])) {
   a<-previous$layer_inputs[[owner]][[component]]
   b<-current$layer_inputs[[owner]][[component]]
   if(!is.data.frame(a))next
   stopifnot(is.data.frame(b))
   stopifnot(identical(a,b))
 }
 before <- association_metrics(previous)$summary
 after <- association_metrics(current)$summary
 fragments <- 0L; escaped <- 0L; degenerate <- 0L; contour_contacts <- 0L
 for (owner in current$layer_geometry) {
  d <- owner$gene_poly
  if (!is.data.frame(d) || !'.component' %in% names(d)) next
  marks <- split(d[d$.component %in% 'boundary',],d$group[d$.component %in% 'boundary'])
  for(mark in marks) {
   if(nrow(mark)<2L)next
   fragments <- fragments+1L
   fragment_length <- sqrt(diff(mark$x)^2+diff(mark$y)^2)
   if(any(fragment_length<=64*.Machine$double.eps*max(1,abs(c(mark$x,mark$y))))) {
     degenerate<-degenerate+1L
     next
   }
   own <- d[d$.component %in% 'polygon' & d$source_row==mark$source_row[1L],]
   polygons <- split(own,own$group)
   on_contour <- FALSE
   for(t in c(.1,.5,.9)) {
    point<-unlist(mark[1L,c('x','y')]+t*(mark[2L,c('x','y')]-mark[1L,c('x','y')]))
    on_or_inside <- function(p) {
      if(ggchord_point_in_polygon(point[1L],point[2L],p))return(TRUE)
      # Source segment boundaries can coincide with a closed contour edge.
      a<-cbind(p$x,p$y);b<-a[c(2:nrow(a),1L),,drop=FALSE];edge<-b-a
      delta<-matrix(point,nrow(a),2,byrow=TRUE)-a
      t<-pmin(1,pmax(0,rowSums(delta*edge)/pmax(rowSums(edge^2),1e-20)))
      nearest<-a+edge*t
      touch <- min(sqrt(rowSums((nearest-matrix(point,nrow(a),2,byrow=TRUE))^2))) <=
        64*.Machine$double.eps*max(1,abs(c(point,a)))
      if(touch) on_contour <<- TRUE
      touch
    }
    if(!any(vapply(polygons,on_or_inside,logical(1)))) {
     escaped<-escaped+1L;break
    }
   }
   if(on_contour)contour_contacts<-contour_contacts+1L
  }
 }
 labels <- current$gene_labels
 old_labels <- previous$gene_labels
 endpoint_rows <- Filter(function(row)isTRUE(row$endpoint),observations$cases[[id]]$axis_labels)
 rows[[id]] <- data.frame(case=id,matching_source_input_rows=TRUE,
  endpoint_strings=length(endpoint_rows),
  hidden_endpoint_strings=sum(vapply(endpoint_rows,function(row)!isTRUE(row$visible),logical(1))),
  anchors_off_source_before=before$anchors_off_own_feature,
  anchors_off_source_after=after$anchors_off_own_feature,
  internal_glyphs_off_tangent=after$internal_glyphs_off_tangent,
  boundary_fragments=fragments,degenerate_boundary_fragments=degenerate,
  boundary_contour_contacts=contour_contacts,
  boundary_fragments_outside_source=escaped,
  inner_source_exterior_before=sum(old_labels$lane>0L & old_labels$feature_label_mode=='external'),
  inner_source_exterior_after=sum(labels$lane>0L & labels$feature_label_mode=='external'))
}
result <- do.call(rbind,rows)
if(any(result$anchors_off_source_after!=0L | result$internal_glyphs_off_tangent!=0L |
 result$boundary_fragments_outside_source!=0L | result$degenerate_boundary_fragments!=0L))
 write.csv(result,'tools/long-term-audit/v0.27.0-compact-contours-candidate.csv',row.names=FALSE)
stopifnot(all(result$anchors_off_source_after==0L),
 all(result$internal_glyphs_off_tangent==0L),all(result$boundary_fragments_outside_source==0L),
 all(result$degenerate_boundary_fragments==0L),all(result$hidden_endpoint_strings==0L))
write.csv(result,'tools/long-term-audit/v0.27.0-compact-contours.csv',row.names=FALSE)
print(result)
