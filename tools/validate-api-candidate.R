# Selected v0.16.4 migration/error/combination gaps. This is not a release run.
pkgload::load_all(quiet=TRUE)
results<-list()
record<-function(case,detail) results[[length(results)+1L]]<<-
 data.frame(case=case,result='PASS',detail=detail)
fails<-function(expr,pattern,absent=NULL) {
 e<-tryCatch(force(expr),error=identity)
 stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e),fixed=TRUE))
 if(!is.null(absent)) stopifnot(!grepl(absent,conditionMessage(e),fixed=TRUE))
 invisible(e)
}
fails(geom_gene_label_repel(label_rotation=90),'geom_gene_label(label_rotation')
fails(geom_feature_label_repel(label_orientation='horizontal'),'geom_feature_label(label_orientation')
fails(geom_feature_label_repel(gene_label_rotation=90),'geom_feature_label(label_rotation')
fails(geom_gene_label_repel(gene_label_mystery=1),'unknown label argument','-> label_mystery')
fails(geom_gene_label_repel(gene_label_layout='aligned'),"label_layout = 'auto'")
fails(geom_gene_label(gene_label_layout='aligned'),'unknown label argument',"-> label_layout")
record('F51_manual_and_unknown_label_controls','Ignored manual geometry now gives an executable fixed-label route; unknown legacy prefixes do not invent parameters')
# Reuse one small build for current migration targets and normal text styles.
s<-data.frame(accver='A',length=1000)
g<-data.frame(accver='A',start=100,end=300,strand='+',anno='gene')
f<-data.frame(accver='A',start=600,end=800,strand='.',type='misc',label='feature')
p<-ggchord(s)+geom_seq()+geom_gene(data=g,width=.05)+
 geom_feature(data=f,shape='block',width=.06,position=position_strand(.1))+
 geom_gene_label(data=g,label_rotation=90)+
 geom_feature_label(data=f,label_orientation='horizontal')+coord_collinear()
b<-ggplot2::ggplot_build(p)
stopifnot(all(vapply(b$data[2:5],nrow,integer(1))>0))
record('current_migration_targets','Unified width/shape/Position and fixed gene/unknown-strand feature labels build together')
# Compare default automatic behavior with the actual previous public wrapper.
ns<-asNamespace('ggchord');prior<-new.env(parent=ns)
source<-system2('git',c('show','6d8e94a:R/geom-gene.R'),stdout=TRUE)
eval(parse(text=source),prior)
automatic<-function(layer) export_ggchord_layout(ggchord(s)+geom_seq()+layer+coord_collinear(),include='labels')
old<-automatic(prior$geom_gene_label_repel(data=g,colour='#123456',size=3))
current<-automatic(geom_gene_label_repel(data=g,colour='#123456',size=3))
stopifnot(identical(old,current))
record('automatic_default_preserved','Normal automatic labels and explicit supported text styles match the previous implementation')
fails(geom_gene(gene_offset=.1),'position_strand')
fails(geom_feature(feature_shape='block'),'shape')
fails(geom_link_ribbon(ribbon_gap=.1),'link_gap')
fails(geom_restriction_site(data=data.frame(accver='A',position=100,enzyme='e'),leader='trunk'),'radial')
fails(ggchord(s,rotation=0),'coord_chord')
fails(theme_ggchord(base_size=10),'element_text')
record('retired_guards','Representative constructor/geom/Position/leader/theme migrations retain clear replacement errors')
# One error model per responsibility; unknown geoms remain delegated upstream.
fails(ggchord(s,candidate_unknown=1),'unused argument')
fails(geom_seq(aes(seq_id=accver)),'accver')
ordinary<-withCallingHandlers(geom_seq(candidate_unknown=1),warning=function(w) invokeRestart('muffleWarning'))
stopifnot(inherits(ordinary,'LayerInstance'))
record('error_responsibilities','Constructor and retired identifiers stop; ordinary geom extras retain ggplot2 delegation')
# Alias and reexport checks are identity checks, not plotting/serialization matrices.
for(pair in list(c('scale_seq_color_manual','scale_seq_colour_manual'),
 c('scale_link_color_manual','scale_link_colour_manual'),
 c('guide_ggchord_colorbar','guide_ggchord_colourbar'))) {
 stopifnot(identical(getExportedValue('ggchord',pair[1]),getExportedValue('ggchord',pair[2])))
}
re<-read.csv('design-roadmap/audit-generated/api_reexports.csv')
stopifnot(all(vapply(seq_len(nrow(re)),function(i) identical(
 getExportedValue('ggchord',re$symbol[i]),getExportedValue(re$package[i],re$symbol[i])),logical(1))))
record('aliases_and_reexports','Representative spelling aliases and all 22 upstream reexports retain identity')
out<-do.call(rbind,results)
write.csv(out,'design-roadmap/audit-generated/candidate-validation-v0.16.4.csv',row.names=FALSE)
print(out,row.names=FALSE)
