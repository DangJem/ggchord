# Explicit synthetic scalability profile, not biological data or release checks.
pkgload::load_all('.',quiet=TRUE,export_all=TRUE)
arg<-grep('^--counts=',commandArgs(trailingOnly=TRUE),value=TRUE)
counts<-if(length(arg))as.integer(strsplit(sub('^--counts=','',arg),',',fixed=TRUE)[[1]]) else c(10L,100L)
stopifnot(all(is.finite(counts)),all(counts>0))
result<-list(scope='Synthetic geometry; one PNG device, no fixed-time acceptance threshold; R heap high-water mark is not whole-process RSS',
 source_parent_sha=system2('git',c('rev-parse','HEAD'),stdout=TRUE),
 source_sha256=as.list(stats::setNames(vapply(list.files('R','[.]R$',full.names=TRUE),function(path)digest::digest(file=path,algo='sha256'),character(1)),list.files('R','[.]R$',full.names=TRUE))),
 device=list(name='png',width_inches=6,height_inches=4,dpi=100),R=R.version.string,results=list())
write<-function()jsonlite::write_json(result,'tools/long-term-audit/ribbon-profile-v0.17.0.json',pretty=TRUE,auto_unbox=TRUE)
for(n in counts) {
 cat('Profiling synthetic ribbons:',n,'\n');flush.console()
 result$results[[as.character(n)]]<-list(status='running');write()
 s<-data.frame(accver=c('A','B'),length=1000)
 i<-seq_len(n)
 h<-data.frame(qaccver='A',saccver='B',qstart=(i*17)%%500+1,
   qend=(i*17)%%500+101,sstart=(i*31)%%500+201,send=(i*31)%%500+351)
 reverse<-i%%2==0;temp<-h$sstart[reverse];h$sstart[reverse]<-h$send[reverse];h$send[reverse]<-temp
 p<-ggchord(s,h,validate='none')+geom_seq()+geom_link_ribbon(gradient=TRUE)+coord_collinear()
 warnings<-character();path<-tempfile(fileext='.png');gc(reset=TRUE)
 value<-tryCatch({
   elapsed<-system.time({
     grDevices::png(path,width=6,height=4,units='in',res=100)
     tryCatch(withCallingHandlers({
       table<-ggplot2::ggplotGrob(p);grid::grid.newpage();grid::grid.draw(table)
       layout<-attr(table,'ggchord_render_layout')
     },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
       finally=grDevices::dev.off())
   })[['elapsed']]
   retained<-length(unique(layout$ribbon_polys$source_row))
   stopifnot(retained==n)
   memory<-gc()
   list(status='rendered',input_rows=n,retained_rows=retained,elapsed_seconds=elapsed,
     PNG_bytes=unname(file.info(path)$size),peak_R_heap_mb=sum(memory[,match('max used',colnames(memory))+1L]),warnings=unique(warnings))
 },error=function(e)list(status='failed',error=conditionMessage(e),warnings=unique(warnings)))
 result$results[[as.character(n)]]<-value;write()
 unlink(path);rm(p,h);if(exists('table',inherits=FALSE))rm(table);if(exists('layout',inherits=FALSE))rm(layout);gc()
 cat(n,value$status,'\n');flush.console()
}
