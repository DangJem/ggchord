# M3 explicit multi-detail development checks; no package/release test suite.
pkgload::load_all(quiet=TRUE)
library(dplyr)
seqs <- seq_data_example |> select(accver,length)
hit <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hit,accver=qaccver,start=pmin(qstart,qend),end=pmax(qstart,qend),detail_panel='query'),
  transmute(hit,accver=saccver,start=pmin(sstart,send),end=pmax(sstart,send),detail_panel='subject'))
regions$region_row <- c('query metadata','subject metadata')
regions$detail_panel_note <- 'preserve metadata'
details <- lapply(1:2,function(i) {
  local <- focus_ggchord_data(seqs,loci=regions[i,])
  ggchord(local$seq_data)+geom_seq(show.legend=FALSE)+coord_collinear(seq_axis='vertical')
}) |> setNames(c('query','subject'))
details$subject <- details$subject+ggplot2::scale_y_reverse()
whole <- ggchord(seqs,hit)+geom_seq(show.legend=FALSE)+geom_link_ribbon(show.legend=FALSE)
x <- magnify_ggchord(whole,regions,details,side=c(query='left',subject='right'))
checks <- lapply(list(c(4,4),c(8,8),c(12,8)),function(size) {
  e <- export_ggchord_fusion(x,width=size[1],height=size[2])
  stopifnot(identical(e$focus_region_data,regions),nrow(e$leaders)==4L,
    nrow(e$connections)==0L,length(e$relationship_data)==0L,
    identical(e$focus_regions$detail_panel,regions$detail_panel))
  for(i in 1:2) {
    fronts <- filter(e$focus_fronts,region_row==i)
    stopifnot(setequal(unique(fronts$panel_id),c('center',regions$detail_panel[i])))
    for(id in unique(fronts$panel_id)) stopifnot(isTRUE(all.equal(
      range(fronts$position_boundary[fronts$panel_id==id]),c(regions$start[i]-1,regions$end[i]))))
    id <- regions$detail_panel[i]
    backbone <- e$source_layouts[[id]]$seq_arcs[[regions$accver[i]]]
    frame <- e$diagnostics[[id]]$render_context
    expected <- sapply(c('x','y'),function(axis) {
      values <- backbone[[axis]][c(1,nrow(backbone))]
      scale <- x$panels[[id]]$plot$scales$get_scales(axis)
      if(!is.null(scale)) values <- scale$transform(values)
      bounds <- frame$panel_bounds_mm[if(axis=='x') 1:2 else 3:4]
      bounds[1]+(values-frame[[paste0(axis,'_range')]][1])/
        diff(frame[[paste0(axis,'_range')]])*diff(bounds)
    })
    actual <- fronts |> filter(panel_id==id) |> slice(c(1,n()))
    stopifnot(max(sqrt(rowSums((as.matrix(actual[c('x','y')])-expected)^2))) < .15)
    ends <- fronts |> filter(panel_id==regions$detail_panel[i]) |> slice(c(1,n()))
    leader <- e$leaders |> filter(region_row==i)
    stopifnot(isTRUE(all.equal(leader$xend,ends$x)),isTRUE(all.equal(leader$yend,ends$y)))
  }
  list(size_in=size,panels=e$panels,regions=e$focus_regions)
})
reject <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
reject(magnify_ggchord(whole,select(regions,-detail_panel),details,side=c(query='left',subject='right')))
bad <- regions;bad$detail_panel[2] <- 'unknown'
reject(magnify_ggchord(whole,bad,details,side=c(query='left',subject='right')))
reject(magnify_ggchord(whole,regions,details,side=c(query='left',subject='left')))
reject(export_ggchord_fusion(magnify_ggchord(whole,regions,details,
  side=c(query='left',subject='right'),center_size=c(.8,.8))))
reject(export_ggchord_fusion(magnify_ggchord(whole,regions,details,
  side=c(query='left',subject='right'),detail_size=list(query=c(.9,.9),subject=c(.9,.9)))))
outside <- tryCatch(export_ggchord_fusion(magnify_ggchord(whole,regions[1,],
  list(query=details$query),side=c(query='right'),center_size=c(.99,.99))),error=identity)
stopifnot(inherits(outside,'error'),grepl('outside canvas',conditionMessage(outside)))
old <- select(regions[1,],-detail_panel)
one <- export_ggchord_fusion(magnify_ggchord(whole,old,details$query))
stopifnot(nrow(one$leaders)==2L,identical(one$focus_region_data,old))
# Same source in two independent details requires panel identity, not guessed routing.
repeat_regions <- regions[c(1,1),];repeat_regions$detail_panel <- c('first','second')
rep <- export_ggchord_fusion(magnify_ggchord(whole,repeat_regions,
  list(first=details$query,second=details$query),side=c(first='top',second='bottom')))
stopifnot(identical(rep$focus_regions$detail_panel,c('first','second')))
ggplot2::ggsave('tools/long-term-audit/v0.25.0-multiple-details.png',x,width=12,height=8,dpi=120)
jsonlite::write_json(list(version='0.25.0.9000',checks=checks,
  single_compatible=TRUE,raw_metadata_preserved=TRUE,repeated_source_explicit=TRUE,
  errors_checked=TRUE,no_release_matrix=TRUE),'tools/long-term-audit/v0.25.0-multiple-details.json',
  pretty=TRUE,auto_unbox=TRUE)
cat('PASS: independent details, source boundaries, three canvas sizes and constraints\n')
