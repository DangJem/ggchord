# Directed source matrix: one PNG recipe, no package tests or device matrix.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional=c('tools/site-flagships.R','tools/home-previews.R',
  'tools/documentation-examples.R','website/manual-figures.json','website/home-previews.json',
  'tools/site-audit/render-dotplot.R'))
pkgload::load_all(quiet=TRUE)
blocks <- function(file) {
  lines<-readLines(file,warn=FALSE);at<-grep('^# ---- [a-z]+ ----$',lines)
  setNames(lapply(seq_len(length(at)-1L),function(i)
    paste(lines[seq.int(at[i]+1L,at[i+1L]-1L)],collapse='\n')),
    sub('^# ---- ([a-z]+) ----$','\\1',lines[at[-length(at)]]))
}
chunks<-c(blocks('tools/documentation-examples.R'),blocks('tools/site-flagships.R'))
for(key in c('sequencepalette','featurepalette','docstyle','showcasedata','chordshowcasegenes',
  'showcasechord','showcasecollinear','flagshipdata','flagshipmulti','fusionoverview','flagshipdot','flagshipfusion'))
  invisible(eval(parse(text=chunks[[key]])))
source('tools/home-previews.R')
# The executable API/signatures/projection remain unchanged: only preparation changes.
old<-parse(text=system2('git',c('show','3a56988:R/coord-dotplot.R'),stdout=TRUE),keep.source=FALSE)
current<-parse('R/coord-dotplot.R',keep.source=FALSE)
stopifnot(length(old)==length(current))
changed<-which(!vapply(seq_along(old),function(i) identical(old[[i]],current[[i]]),logical(1)))
stopifnot(length(changed)==1L,identical(old[[changed]][[2]],as.name('ggchord_prepare_dotplot')),
  identical(old[[changed]][[3]][[2]],current[[changed]][[3]][[2]]))
original<-new.env(parent=asNamespace('ggchord'));invisible(eval(old,original))
source_layout<-original$ggchord_dotplot_geometry(flag_dotplot)
now<-export_ggchord_layout(flag_dotplot,original_data=TRUE)
old_alignment<-dplyr::bind_rows(source_layout$alignment)
stopifnot(nrow(now$alignment)==32L,identical(old_alignment,now$alignment),
  identical(source_layout$axis_layout,now$axis_layout),
  sum(dotplot_cells$HSP)==32L,nrow(dotplot_pair_counts)==6L,
  !any(flagship_hsp$qaccver==flagship_hsp$saccver))
# Reversed axis and inverse mapping keep real native positions; explicit scales win.
reversed<-ggchord(flagship_sequences,flagship_hsp)+geom_alignment_segment()+
  coord_dotplot(x_seqs=flagship_sequences$accver,y_seqs=flagship_sequences$accver,x_orientation=-1)
d<-export_ggchord_layout(reversed)
stopifnot(nrow(d$alignment)==32L,all(d$axis_layout$orientation[d$axis_layout$axis=='x']==-1))
id<-flagship_hsp$qaccver[1];position<-flagship_hsp$qstart[1]
p<-map_ggchord_position(reversed,id,position,axis='x')
back<-map_ggchord_position(reversed,id,p$x,inverse=TRUE,axis='x')
stopifnot(back$position==position)
custom<-flag_dotplot+ggplot2::scale_x_continuous(breaks=100,labels='custom')
built<-ggplot2::ggplot_build(custom)
stopifnot(identical(built$plot$scales$get_scales('x')$breaks,100),
  identical(built$plot$scales$get_scales('x')$labels,'custom'))
# Existing registered single-pair Dotplot remains a real fresh single-PNG output.
for(key in c('dotplotdata','dotplotview')) invisible(eval(parse(text=chunks[[key]])))
plots<-read.csv('website/assets/plots.csv',stringsAsFactors=FALSE)
fingerprints<-jsonlite::fromJSON(system2('python3','tools/documentation-fingerprints.py',stdout=TRUE))
record<-list(status='started',baseline='3a56988',engine_sha256=as.list(audit$engine),
  source_sha256=as.list(audit$inputs),changed_function='ggchord_prepare_dotplot',
  other_expressions_and_signatures_identical=TRUE,
  checks=list(original_HSPs=32L,alignment_and_axis_fields_identical=TRUE,
    directed_pairs=6L,self_HSPs=0L,reciprocal_pairs=0L,
    reversed_axis_retains_32=TRUE,inverse_mapping=TRUE,explicit_position_scale_wins=TRUE),
  pair_counts=dotplot_pair_counts,figures=list(),
  scope='Multi-sequence cell/axis decoration only; real directed records retained; no browser, package or release checks.')
full<-jsonlite::fromJSON('website/manual-figures.json',simplifyVector=FALSE)
home<-jsonlite::fromJSON('website/home-previews.json',simplifyVector=FALSE)
rows<-c(Filter(function(r) r$id=='flagdotplot',full),Filter(function(r) r$id=='homedotplot',home))
i<-match('viewdotplot',plots$id)
rows[[3]]<-as.list(plots[i,]);rows[[3]]$object<-'p_viewdotplot'
for(row in rows) {
  cat('Render:',row$id,'\n');flush.console()
  file<-view_ggchord(get(row$object),width=row$width,height=row$height,
    units='in',device='png',dpi=row$dpi,bg=row$background %||% '#FFFCF5',viewer='none')
  stopifnot(file.copy(file,row$file,overwrite=TRUE))
  audit_assert_source_snapshot(audit)
  record$figures[[row$id]]<-list(file=row$file,sha256=digest::digest(file=row$file,algo='sha256',serialize=FALSE),
    md5=unname(tools::md5sum(row$file)),width=row$width,height=row$height,dpi=row$dpi,status='actual PNG rendered')
}
plots$md5[i]<-record$figures$viewdotplot$md5
plots$code_md5[i]<-fingerprints$code
plots$dependency_md5[i]<-fingerprints$assets$viewdotplot
plots$generator[i]<-'tools/site-audit/render-dotplot.R'
record$figures$viewdotplot$dependency_md5<-plots$dependency_md5[i]
write.csv(plots,'website/assets/plots.csv',row.names=FALSE)
record$status<-'complete'
jsonlite::write_json(record,'tools/site-audit/dotplot-rendering.json',auto_unbox=TRUE,pretty=TRUE)
cat('PASS: 32 original HSPs, six directed pairs, mapping/scale precedence, three PNGs\n')
