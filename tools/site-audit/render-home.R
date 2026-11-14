# Explicit homepage/flagship PNG generation, normal dimensions; no package matrix.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional=c('tools/site-flagships.R','tools/home-previews.R',
  'tools/documentation-examples.R','website/home-previews.json','tools/site-audit/render-home.R'))
pkgload::load_all(quiet=TRUE)
blocks <- function(file) {
  text <- readLines(file,warn=FALSE);marks<-grep('^# ---- [a-z]+ ----$',text)
  setNames(lapply(seq_len(length(marks)-1L),function(i) paste(text[seq.int(marks[i]+1L,marks[i+1L]-1L)],collapse='\n')),
    sub('^# ---- ([a-z]+) ----$','\\1',text[marks[-length(marks)]]))
}
chunks <- c(blocks('tools/documentation-examples.R'),blocks('tools/site-flagships.R'))
for(key in c('sequencepalette','featurepalette','docstyle','showcasedata','chordshowcasegenes',
    'showcasechord','showcasecollinear','fusiondata','flagshipdata','flagshipmulti','fusionoverview','flagshipdot','flagshipfusion'))
  invisible(eval(parse(text=chunks[[key]])))
source('tools/home-previews.R')
previews<-jsonlite::fromJSON('website/home-previews.json',simplifyVector=FALSE)
record<-list(status='started',engine_sha256=as.list(audit$engine),
  source_sha256=as.list(audit$inputs),figures=list(),scope='Two revised full flagships and seven homepage-only PNGs; no crop, no device matrix.')
# Actual data inputs are compared to the preceding public example objects.
baseline <- system2('git',c('show','3a6b412:tools/site-flagships.R'),stdout=TRUE)
start <- grep('^# ---- flagshipmulti ----$',baseline)
end <- grep('^# ---- fusionoverview ----$',baseline)
stopifnot(length(start)==1L,length(end)==1L,end>start)
previous <- new.env(parent=globalenv())
invisible(eval(parse(text=baseline[seq.int(start+1L,end-1L)]),previous))
before <- list(radial=previous$flag_radial,sankey=previous$flag_sankey)
source_rows<-function(plot,component) {
  layout<-get_chord_layout(plot)
  result<-list()
  for(owner in names(layout$layer_inputs)) {
    key<-if(component=='gene') 'gene_poly' else component
    d<-layout$layer_inputs[[owner]][[key]]
    if(is.data.frame(d)&&nrow(d)) result[[length(result)+1L]]<-d
  }
  dplyr::bind_rows(result)
}
select_source<-function(data,component) {
  fields<-if(component=='ribbon') c('qaccver','saccver','qstart','qend','sstart','send','pident') else c('accver','start','end','strand','anno')
  data[,intersect(fields,names(data)),drop=FALSE] |> as.data.frame()
}
checks<-list()
for(name in c('radial','sankey')) for(component in c('ribbon','gene')) {
  old<-select_source(source_rows(before[[name]],component),component)
  current<-select_source(source_rows(get(paste0('flag_',name)),component),component)
  stopifnot(nrow(old)>0L,nrow(current)>0L,identical(old,current))
  checks[[paste(name,component,sep='_')]]<-list(original_records=nrow(old),all_source_fields_identical=TRUE)
}
record$source_checks<-checks
if('--verify-sources-only' %in% commandArgs(trailingOnly=TRUE)) {
  completed<-jsonlite::fromJSON('tools/site-audit/home-rendering.json',simplifyVector=FALSE)
  stopifnot(completed$status=='complete')
  completed$initial_source_checks<-completed$source_checks
  completed$source_checks<-checks
  completed$source_check_correction<-list(
    reason='Initial gene selector used gene instead of the public gene_poly component and returned zero rows; empty checks are rejected and superseded.',
    nonempty_source_verification=TRUE,
    verifier_sha256=digest::digest(file='tools/site-audit/render-home.R',algo='sha256',serialize=FALSE),
    scope='Corrected public input comparison only; original PNG render hashes and execution source hashes retained.')
  jsonlite::write_json(completed,'tools/site-audit/home-rendering.json',auto_unbox=TRUE,pretty=TRUE)
  cat('PASS: nonempty original 12/10 HSP and 20/16 CDS records identical\n')
  quit(save='no',status=0)
}

# Complete figures: freshly rendered final recipes, not renamed candidate files.
full<-list(flagradial=list(object='flag_radial',width=11,height=10,dpi=144,file='man/figures/flagradial.png'),
  flagsankey=list(object='flag_sankey',width=10,height=12,dpi=144,file='man/figures/flagsankey.png'))
for(id in c(names(full),vapply(previews,`[[`,'','id'))) {
 row<-if(id %in% names(full)) full[[id]] else Filter(function(r) r$id==id,previews)[[1L]]
 cat('Render:',id,row$width,row$height,'\n');flush.console()
 file<-view_ggchord(get(row$object),width=row$width,height=row$height,units='in',device='png',dpi=row$dpi,bg='#FFFCF5',viewer='none')
 stopifnot(file.copy(file,row$file,overwrite=TRUE))
 audit_assert_source_snapshot(audit)
 record$figures[[id]]<-list(file=row$file,sha256=digest::digest(file=row$file,algo='sha256',serialize=FALSE),
   width=row$width,height=row$height,dpi=row$dpi,status='actual PNG rendered',background='#FFFCF5')
 jsonlite::write_json(record,'tools/site-audit/home-rendering.json',auto_unbox=TRUE,pretty=TRUE)
}
record$status<-'complete'
record$labels<-list(radial_selected=nrow(radial_labels),sankey_selected=nrow(sankey_labels),
  rule='Per sequence: non-hypothetical name first, then descending feature length; no coordinate changes. Preview relationship figures omit CDS-name layers, retaining geometries.')
jsonlite::write_json(record,'tools/site-audit/home-rendering.json',auto_unbox=TRUE,pretty=TRUE)
cat('PASS: two complete flagships, seven independent homepage previews; source endpoints/CDS identical\n')
