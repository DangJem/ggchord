# Explicit normal-PNG regression for the reported native reference replay defect.
# No package tests or release device matrix are run.
pkgload::load_all(quiet=TRUE)
folder <- tempfile('ggchord-replay-');dir.create(folder)
results <- list()
for (topic in c('ggchord_fusion','fusion_panel','fusion_ribbon','magnify_ggchord')) {
  code_file <- tempfile(fileext='.R')
  tools::Rd2ex(file.path('man',paste0(topic,'.Rd')),out=code_file)
  plots <- list()
  handler <- evaluate::new_output_handler(value=function(x,visible) {
    if(visible && inherits(x,'ggchord_fusion')) plots[[length(plots)+1L]] <<- x
    invisible(NULL)
  })
  evaluate::evaluate(paste(readLines(code_file,warn=FALSE),collapse='\n'),
    envir=new.env(parent=globalenv()),new_device=FALSE,output_handler=handler,stop_on_error=2)
  stopifnot(length(plots)>0)
  for(i in seq_along(plots)) {
    stem <- paste0(topic,'-',i)
    ragg::agg_png(file.path(folder,paste0(stem,'-capture.png')),width=480,height=480,res=72,bg='white')
    grDevices::dev.control(displaylist='enable')
    print(plots[[i]]);captured <- grDevices::recordPlot();grDevices::dev.off()
    replay <- file.path(folder,paste0(stem,'-replay.png'))
    ragg::agg_png(replay,width=1400,height=866,res=192,bg='white')
    grDevices::replayPlot(captured);grDevices::dev.off()
    direct <- file.path(folder,paste0(stem,'-direct.png'))
    ggplot2::ggsave(direct,plots[[i]],width=1400/192,height=866/192,dpi=192,bg='white')
    stopifnot(unname(tools::md5sum(replay))==unname(tools::md5sum(direct)))
    results[[stem]] <- list(recorded_replay_matches_direct_export=TRUE,
      replay_sha256=digest::digest(file=replay,algo='sha256'))
    cat('PASS:',stem,'replay matches direct PNG\n');flush.console()
  }
}
record <- list(status='passed',scope='Actual Rd examples on one capture/export device pair; no release matrix',
  source_sha256=digest::digest(file='R/fusion.R',algo='sha256'),cases=results)
jsonlite::write_json(record,'tools/site-audit/fusion-replay-v0.28.0.json',auto_unbox=TRUE,pretty=TRUE)
