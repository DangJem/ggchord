# Small synthetic geometry controls, not biological observations or release checks.
source('tools/audit-source-snapshot.R')
audit<-audit_source_snapshot(additional='tools/long-term-audit/check-feature-leader-context.R')
pkgload::load_all(quiet=TRUE)
seqs<-data.frame(accver=c('one','two'),length=1000)
features<-data.frame(accver=rep(seqs$accver,each=2),start=rep(c(100,500),2),
 end=rep(c(350,800),2),strand=rep(c('+','-'),2),type='CDS',
 anno=c('first','second','third','fourth'))
coordinates<-list(collinear=coord_collinear(),chord=coord_chord(),circular=coord_circular())
cases<-list()
for(name in names(coordinates)) {
 s<-if(name=='circular')seqs[1,,drop=FALSE] else seqs
 f<-features[features$accver %in% s$accver,,drop=FALSE]
 plot<-ggchord(s,feature_data=f)+geom_seq()+geom_feature()+
   geom_feature_label_repel()+coordinates[[name]]
 built<-ggplot2::ggplot_build(plot)
 owner<-which(vapply(built$plot$layers,function(layer)
   isTRUE(layer$ggchord_params$is_feature_label),logical(1)))
 d<-built$data[[owner]]
 stopifnot(nrow(d)>0L,'.feature_internal_exit' %in% names(d))
 if(name!='circular')stopifnot(!any(d$.feature_internal_exit %in% TRUE))
 cases[[name]]<-list(build_passed=TRUE,inside_to_exterior_flag_present=TRUE)
}
empty<-ggchord(seqs,feature_data=features)+geom_seq()+geom_feature()+
 geom_feature_label_repel(data=features[FALSE,,drop=FALSE])+coord_collinear()
built<-ggplot2::ggplot_build(empty)
owner<-which(vapply(built$plot$layers,function(layer)
 isTRUE(layer$ggchord_params$is_feature_label),logical(1)))
stopifnot(nrow(built$data[[owner]])==0L)
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(cases=cases,empty_label_layer_passed=TRUE,
 runtime_sha256=as.list(audit$engine)),
 'tools/long-term-audit/v0.27.0-feature-leader-context-results.json',pretty=TRUE,auto_unbox=TRUE)
cat('Feature label context and empty-data builds passed\n')
