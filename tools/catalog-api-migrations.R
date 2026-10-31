# Current migration catalogue; preserves the v0.16.0 historical review.
# Literal retired-argument maps are read from R AST, never executed as user code.
pkgload::load_all(quiet=TRUE)
out<-'design-roadmap/audit-generated'
prior<-read.csv('design-roadmap/audit-reviewed-details/retired-entries.csv',
 na.strings=character(),check.names=FALSE,stringsAsFactors=FALSE)
revision<-system2('git',c('rev-parse','HEAD'),stdout=TRUE)
rows<-data.frame(function_name=prior$function_name,old_entry=prior$old_entry,
 entry_kind=ifelse(grepl('^aes\\(|boundaries',prior$old_entry),'mapping_or_field',
  ifelse(grepl('=',prior$old_entry,fixed=TRUE),'value','argument')),
 current_status='REMOVED',replacement=prior$replacement,
 removed_by_version=ifelse(grepl('v0.10.0',prior$removed_or_planned_version),'0.10.0',
  ifelse(grepl('v0.13.1',prior$removed_or_planned_version),'0.13.1','0.11.0')),
 version_evidence=prior$removed_or_planned_version,
 guard=prior$evidence,forwarded_by='',decision='KEEP_MIGRATION_ERROR; no compatibility execution',
 prior_replacement=prior$replacement,prior_status=prior$current_status,
 prior_revision=prior$source_revision,source_revision=revision,
 runtime_evidence='SOURCE_REVIEW; representative finite triggers only',stringsAsFactors=FALSE)
# A current name in prose can be replaced without rewriting geom_* function names.
rows$replacement<-gsub('(?<!geom_)gene_label_(?!repel)','label_',rows$replacement,perl=TRUE)
add<-function(fn,old,replacement,version,guard,kind='argument',status='REMOVED',forwarded='') {
 key<-paste(rows$function_name,rows$old_entry,sep='\t');at<-match(paste(fn,old,sep='\t'),key)
 if(is.na(at)) {
  at<-nrow(rows)+1L
  rows[at,]<<-list(fn,old,kind,status,replacement,version,
   if(version=='0.16.1') 'NEWS v0.16.1 and implemented named migration guard' else
    if(version=='0.16.4') 'F51: previously ignored unsupported manual control; not a removed supported formal' else
     'Current guard/source and existing NEWS category; first removal tag not re-established',
   guard,forwarded,'KEEP_MIGRATION_ERROR; no compatibility execution','','',
   'new current-source record',revision,'SOURCE_REVIEW; representative finite triggers only')
 } else {
  rows$replacement[at]<<-replacement;rows$guard[at]<<-guard
  rows$current_status[at]<<-status;rows$forwarded_by[at]<<-forwarded
  if(prior$current_status[match(paste(fn,old,sep='\t'),paste(prior$function_name,prior$old_entry,sep='\t'))] %in% 'DEPRECATED_STILL_ACCEPTED') {
   rows$removed_by_version[at]<<-'0.16.1'
   rows$version_evidence[at]<<-'NEWS v0.16.1 supersedes former v1.0.0 removal date'
  }
 }
}
literal<-function(x,constants) {
 if(is.character(x)) return(x)
 if(is.symbol(x)) return(constants[[as.character(x)]])
 if(is.call(x) && identical(x[[1]],as.name('c')) &&
    all(vapply(as.list(x)[-1],is.character,logical(1)))) return(eval(x,baseenv()))
 NULL
}
for(file in list.files('R',pattern='[.]R$',full.names=TRUE)) for(expr in parse(file)) {
 if(!is.call(expr)||!identical(expr[[1]],as.name('<-'))||!is.symbol(expr[[2]])||
    !is.call(expr[[3]])||!identical(expr[[3]][[1]],as.name('function'))) next
 owner<-as.character(expr[[2]]);constants<-list()
 walk<-function(x) {
  if(!is.call(x)) return(invisible(NULL))
  if(identical(x[[1]],as.name('<-')) && is.symbol(x[[2]])) {
   val<-literal(x[[3]],constants)
   if(!is.null(val)) constants[[as.character(x[[2]])]]<<-val
  }
  if(identical(x[[1]],as.name('ggchord_reject_retired'))) {
   map<-literal(x[[4]],constants)
   caller<-literal(x[[3]],constants)
   fns<-if(owner=='ggchord_theme_from_environment') grep('^theme_ggchord',getNamespaceExports('ggchord'),value=TRUE) else
    if(!is.null(caller)) sub('\\(\\)$','',caller) else owner
   if(!is.null(map)) for(fn in fns) for(old in names(map)) {
    version<-if(old %in% c('gene_width','feature_width','feature_shape','gene_offset','feature_offset','ribbon_gap','ribbon_ctrl_point')) '0.16.1' else '0.11.0'
    forwarded<-switch(fn,geom_feature='geom_feature_plasmid',
     geom_link_ribbon='stat_ribbon_bundle;stat_ribbon_density',
     geom_gene_label_repel='geom_feature_label_repel;geom_primer_label_repel','')
    add(fn,old,unname(map[[old]]),version,paste0(file,'::',owner,'/ggchord_reject_retired'),forwarded=forwarded)
   }
  }
  for(value in as.list(x)[-1]) if(!rlang::is_missing(value)) walk(value)
 }
 walk(expr[[3]][[3]])
}
# Prefix migrations are generated from actual public formals and the explicit exceptions.
for(fn in c('geom_gene_label','geom_gene_label_repel','geom_feature_label','geom_feature_label_repel')) {
 names<-names(formals(getExportedValue('ggchord',fn)))
 current<-grep('^label_',names,value=TRUE)
 for(name in current) add(fn,paste0('gene_',name),name,'0.16.1',
  'R/geom-gene.R::ggchord_reject_label_prefix; current public formals')
 add(fn,'gene_label_size','size','0.11.0','R/geom-gene.R::ggchord_reject_label_prefix')
 if('max_overlaps' %in% names) add(fn,'gene_label_max_overlaps','max_overlaps','0.16.1',
  'R/geom-gene.R::ggchord_reject_label_prefix')
 if(grepl('repel',fn)) {
  fixed<-if(grepl('feature',fn)) 'geom_feature_label' else 'geom_gene_label'
  for(name in c('label_rotation','label_radial_offset','label_circum_offset','label_circum_limit','label_orientation'))
   add(fn,paste0('gene_',name),paste0(fixed,'(',name,' = ...) for manual placement'),
    '0.16.1','R/geom-gene.R::ggchord_reject_label_prefix')
  for(name in c('label_rotation','label_radial_offset','label_circum_offset','label_circum_limit','label_orientation','label_overlap'))
   add(fn,name,paste0(fixed,'(',name,' = ...), or label_layout for automatic placement'),
    '0.16.4','R/geom-gene.R::ggchord_reject_manual_label_args',status='UNSUPPORTED_MANUAL_CONTROL',
    forwarded=if(grepl('feature',fn)) 'geom_primer_label_repel' else '')
  add(fn,'gene_label_segment','select label_layout; use label_segment_linetype only for supported leader style',
   '0.11.0','R/geom-gene.R::ggchord_reject_label_prefix / automatic layout guard')
 }
}
add('geom_gene_label_repel','gene_label_layout="aligned"','label_layout = "auto"','0.11.0',
 'R/geom-gene.R::ggchord_reject_label_prefix / aligned value migration',kind='value')
add('geom_gene_label_repel','label_layout="aligned"','label_layout = "auto"','0.11.0',
 'R/geom-gene.R::ggchord_label_repel_layer aligned value guard',kind='value')
add('geom_restriction_site','leader="trunk"','leader = "radial"','0.16.1',
 'R/restriction-sites.R::geom_restriction_site explicit value guard',kind='value')
# Earlier removed public functions remain migration facts, not restored aliases.
for(item in list(
 c('ggchord','geom_axis','automatic axes + theme_ggchord(axis.*); scale_seq_position_continuous() for breaks','0.11.0'),
 c('theme_ggchord','theme_ggchord_elements','theme_ggchord() / theme() / element_*()','0.11.0'),
 c('geom_seq_label','geom_seq_group_label','remove sequence grouping; compose explicit sequence/label layers','0.10.0'),
 c('scale_seq_colour_manual','scale_group_colour_manual','remove grouping; explicitly map sequence categories with aes(seq_colour = ...) and scale_seq_colour_manual()','0.10.0'),
 c('scale_seq_color_manual','scale_group_color_manual','remove grouping; explicitly map sequence categories with aes(seq_colour = ...) and scale_seq_colour_manual()','0.10.0'),
 c('geom_link_ribbon','geom_ribbon_highlight','add another geom_link_ribbon(data = subset, fill = ...) layer','0.12.0'),
 c('geom_link_ribbon','geom_ribbon','geom_link_ribbon(); ribbon aesthetics and scales keep their meanings','0.13.1'))) {
 stopifnot(!item[2] %in% getNamespaceExports('ggchord'))
 add(item[1],paste0('ggchord::',item[2],'()'),item[3],item[4],
  'NEWS historical removal entry and current NAMESPACE absence; no compatibility implementation',kind='function')
 at<-which(rows$function_name==item[1] & rows$old_entry==paste0('ggchord::',item[2],'()'))
 rows$version_evidence[at]<-if(item[2]=='geom_ribbon')
  'NEWS v0.12 deprecation and v0.13 development removal; v0.13.1 is the published successor, v0.13.0 was never released' else
  paste0('NEWS ',item[4],' explicit removed-function entry; current namespace confirms absence')
}
# No v1.0.0 wait remains: the two offsets were actually removed in v0.16.1.
for(fn in c('geom_gene','geom_feature')) {
 old<-if(fn=='geom_gene') 'gene_offset' else 'feature_offset'
 at<-which(rows$function_name==fn & rows$old_entry==old)
 rows$replacement[at]<-'position = position_strand(offset = ...); explicit minus-strand values are signed'
}
# Guards remain even though the implementation is gone, to prevent partial matching/ignored dots.
rows$first_removed_version<-ifelse(grepl('未逐tag|未逐tag|first removal',rows$version_evidence),
 'NOT_REESTABLISHED; removed_by_version is an upper bound',rows$removed_by_version)
rows$first_removed_version[rows$current_status=='UNSUPPORTED_MANUAL_CONTROL']<-
 'NOT_APPLICABLE; unsupported manual controls were previously ignored'
rows$decision[rows$current_status=='UNSUPPORTED_MANUAL_CONTROL']<-
 'REJECT_UNSUPPORTED_CONTROL; choose a fixed label or automatic layout'
rows$rejection_change_version<-ifelse(rows$current_status=='UNSUPPORTED_MANUAL_CONTROL','0.16.4','retained prior migration error')
rows$removed_by_version[rows$current_status=='UNSUPPORTED_MANUAL_CONTROL']<-'NOT_APPLICABLE'
rows$first_removed_version[rows$old_entry=='ggchord::geom_ribbon()']<-
 'v0.13 development removal (unreleased); published successor v0.13.1'
# A parent's v0.11 cleanup does not establish the first tag of a later wrapper.
wrapped<-grepl('^theme_ggchord_',rows$function_name) & rows$removed_by_version=='0.11.0'
rows$removed_by_version[wrapped]<-'0.16.0'
rows$first_removed_version[wrapped]<-'NOT_REESTABLISHED; parent syntax cleanup is not a wrapper first-removal tag'
rows$version_evidence[wrapped]<-'v0.16.0 review confirms rejection; v0.11.0 describes the parent syntax family, not the later wrapper creation/removal tag'
feature_prefix<-rows$function_name %in% c('geom_feature_label','geom_feature_label_repel') &
 grepl('^gene_label_',rows$old_entry)
rows$removed_by_version[feature_prefix]<-'0.16.1'
rows$first_removed_version[feature_prefix]<-'NOT_REESTABLISHED; shared legacy name, per-wrapper first supported usage/removal not asserted'
rows$version_evidence[feature_prefix]<-'Current prefix guard and v0.16.1 canonical naming; inherited old gene names are not asserted to have been supported feature formals'

rows<-rows[order(rows$function_name,rows$old_entry),]
stopifnot(!anyDuplicated(paste(rows$function_name,rows$old_entry)),
 all(rows$function_name %in% c(getNamespaceExports('ggchord'),'all public data/mapping boundaries')),
 !any(rows$current_status=='DEPRECATED_STILL_ACCEPTED'))
# A removed named argument cannot accidentally be a current formal again.
for(i in which(rows$entry_kind=='argument' & rows$current_status=='REMOVED')) {
 fn<-rows$function_name[i]
 if(fn=='all public data/mapping boundaries') next
 stopifnot(!rows$old_entry[i] %in% names(formals(getExportedValue('ggchord',fn))))
}
write.csv(rows,file.path(out,'candidate-migrations-v0.16.4.csv'),row.names=FALSE)
cat(nrow(rows),'current migration/error records; historical review retained separately.\n')
