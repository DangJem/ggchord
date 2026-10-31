# Consolidate reviewed v0.16.0-3 contracts into the v0.16.4 candidate baseline.
# snapshot writes the reviewed candidate; check detects public signature drift.
# No plotting, package tests, builds/checks, vignettes or device matrix.
pkgload::load_all(quiet=TRUE)
mode <- commandArgs(trailingOnly=TRUE)
if(!length(mode)) mode <- 'check'
stopifnot(length(mode)==1L,mode %in% c('snapshot','check'))
out <- 'design-roadmap/audit-generated'
ns <- asNamespace('ggchord')
revision <- system2('git',c('rev-parse','HEAD'),stdout=TRUE)
exports <- sort(getNamespaceExports('ggchord'))
inputs<-sort(unique(c(list.files('R',pattern='[.]R$',full.names=TRUE),
 list.files('man',pattern='[.]Rd$',full.names=TRUE),'NAMESPACE','DESCRIPTION','NEWS.md',
 'tools/freeze-api.R','tools/catalog-api-migrations.R','tools/validate-api-candidate.R',
 'design-roadmap/audit-generated/api_arguments.csv','design-roadmap/audit-generated/api_exports.csv',
 'design-roadmap/audit-generated/api_reexports.csv',
 'design-roadmap/audit-generated/data-parameters-v0.16.2.csv',
 'design-roadmap/audit-generated/engine-parameters-v0.16.3.csv',
 'design-roadmap/audit-generated/grammar-parameters-v0.16.1.csv',
 'design-roadmap/audit-generated/grammar-compatibility-v0.16.1.csv',
 'design-roadmap/audit-generated/grammar-aesthetics-v0.16.1.csv',
 'design-roadmap/audit-reviewed-details/retired-entries.csv',
 'design-roadmap/API-CANDIDATE-v0.16.4.md')))
before<-tools::md5sum(inputs)

reexports <- read.csv(file.path(out,'api_reexports.csv'),na.strings=character())
old_exports <- read.csv(file.path(out,'api_exports.csv'),na.strings=character())
base <- read.csv(file.path(out,'api_arguments.csv'),na.strings=character(),check.names=FALSE)
base <- base[base$is_exported,]
data <- read.csv(file.path(out,'data-parameters-v0.16.2.csv'),na.strings=character(),check.names=FALSE)
engine <- read.csv(file.path(out,'engine-parameters-v0.16.3.csv'),na.strings=character(),check.names=FALSE)
grammar <- read.csv(file.path(out,'grammar-parameters-v0.16.1.csv'),na.strings=character())
rd_text <- function(x) gsub('[[:space:]]+',' ',trimws(paste(unlist(x,use.names=FALSE),collapse='')))
docs <- list()
for(file in list.files('man',pattern='[.]Rd$',full.names=TRUE)) {
 rd <- tools::parse_Rd(file)
 tags<-vapply(rd,function(x) attr(x,'Rd_tag') %||% '',character(1))
 args<-list()
 for(block in rd[tags=='\\arguments']) for(item in block)
  if(identical(attr(item,'Rd_tag'),'\\item')) for(name in trimws(strsplit(rd_text(item[[1]]),',',fixed=TRUE)[[1]])) args[[name]]<-rd_text(item[[2]])
 value<-if(any(tags=='\\value')) rd_text(rd[tags=='\\value'][[1]]) else 'Upstream object; see its original package documentation'
 for(name in vapply(rd[tags=='\\alias'],rd_text,character(1))) docs[[name]]<-list(file=file,args=args,value=value)
}
# Direct binding aliases, not wrapper names that merely have similar behavior.
aliases<-character()
for(file in list.files('R',pattern='[.]R$',full.names=TRUE)) for(x in parse(file))
 if(is.call(x) && identical(x[[1]],as.name('<-')) && is.symbol(x[[2]]) && is.symbol(x[[3]]) &&
    as.character(x[[2]]) %in% exports) aliases[as.character(x[[2]])]<-as.character(x[[3]])
for(name in names(aliases)) stopifnot(identical(getExportedValue('ggchord',name),
 getExportedValue('ggchord',aliases[[name]])))
default <- function(args,i) {
 if(names(args)[i]=='...') return('<ellipsis>')
 text<-paste(deparse(args[[i]],width.cutoff=500),collapse=' ')
 if(nzchar(text)) text else '<required>'
}
family <- function(name) {
 if(name %in% reexports$symbol) return('upstream reexport')
 if(name %in% c('get_chord_layout','export_ggchord_layout','view_ggchord','diagnose_ggchord')) return('output')
 if(name %in% unique(data[['function']]) && !startsWith(name,'stat_')) return('data workflow')
 if(name=='ggchord') return('constructor')
 sub('_.*$','',name)
}
manifest <- do.call(rbind,lapply(exports,function(name) {
 value<-getExportedValue('ggchord',name)
 upstream<-match(name,reexports$symbol)
 if(!is.na(upstream)) stopifnot(identical(value,getExportedValue(reexports$package[upstream],name)))
 args<-if(is.function(value)) formals(value) else list()
 signature<-if(is.function(value)) paste(paste0(names(args),'=',vapply(seq_along(args),function(i) default(args,i),character(1))),collapse='; ') else '<upstream ggproto object>'
 row<-old_exports[old_exports$symbol==name,,drop=FALSE]
 doc<-docs[[name]]
 data.frame(name=name,family=family(name),kind=if(is.function(value)) 'function' else 'upstream object',
  ownership=if(is.na(upstream)) 'ggchord' else reexports$package[upstream],
  alias_of=if(name %in% names(aliases)) unname(aliases[name]) else '',
  signature=signature,documentation=doc$file %||% 'upstream namespace documentation',
  current_decision='KEEP_CANDIDATE',candidate_version='0.16.4',source_revision=revision,
  responsibility=if(is.na(upstream)) 'Existing public task retained; current Rd and family contract are authoritative' else
   'Exact original upstream object; upstream signatures/behavior follow the installed supported dependency',
  runtime_evidence='Inherited finite evidence from v0.16.1-3; not every entry executed in v0.16.4')
}))
if(mode=='check') {
 path<-file.path(out,'candidate-exports-v0.16.4.csv')
 if(!file.exists(path)) stop('No candidate snapshot; review contracts before explicitly taking a snapshot')
 prior<-read.csv(path,stringsAsFactors=FALSE,na.strings=character())
 added<-setdiff(manifest$name,prior$name);removed<-setdiff(prior$name,manifest$name)
 shared<-intersect(prior$name,manifest$name)
 fields<-c('signature','ownership','alias_of','kind')
 changed<-shared[vapply(shared,function(name) {
  selected<-if(manifest$ownership[manifest$name==name]=='ggchord') fields else c('ownership','alias_of')
  a<-prior[prior$name==name,selected,drop=FALSE];b<-manifest[manifest$name==name,selected,drop=FALSE]
  !identical(unname(as.list(a)),unname(as.list(b)))
 },logical(1))]
 if(length(c(added,removed,changed))) stop('Candidate API drift: added=',paste(added,collapse=','),
  '; removed=',paste(removed,collapse=','),'; changed=',paste(changed,collapse=','),
  '. Update the decision/migration record before explicitly refreshing the candidate.')
 prior_s3<-read.csv(file.path(out,'candidate-s3-v0.16.4.csv'),stringsAsFactors=FALSE)
 stopifnot(setequal(prior_s3$registration,grep('^S3method',readLines('NAMESPACE'),value=TRUE)))
 cat(length(exports),'candidate exports retain their signatures, aliases and upstream identity.\n')
 quit(save='no',status=0L)
}
stopifnot(setequal(exports,old_exports$symbol))
# Latest reviewed semantics take precedence. Historical evidence is kept separately.
select_semantics <- function(name,parameter) {
 for(item in list(list(table=engine,path='engine-parameters-v0.16.3.csv',version='0.16.3'),
                  list(table=data,path='data-parameters-v0.16.2.csv',version='0.16.2'))) {
  at<-which(item$table[['function']]==name & item$table$parameter==parameter)
  if(length(at)) return(list(row=item$table[at,,drop=FALSE],path=item$path,version=item$version,
   prior_name=name,prior_parameter=parameter))
 }
 legacy<-parameter;legacy_name<-name
 if(parameter=='width' && name %in% c('geom_gene','geom_feature','geom_feature_plasmid')) legacy<-if(name=='geom_gene') 'gene_width' else 'feature_width'
 if(startsWith(parameter,'label_') && name %in% c('geom_gene_label','geom_gene_label_repel')) legacy<-paste0('gene_',parameter)
 if(name=='geom_link_ribbon' && parameter %in% c('link_gap','link_ctrl_point')) legacy<-sub('link_','ribbon_',parameter)
 at<-which(base[['function']]==legacy_name & base$parameter==legacy)
 if(!length(at) && name %in% c('geom_feature_label','geom_feature_label_repel') && startsWith(parameter,'label_')) {
  legacy_name<-if(name=='geom_feature_label') 'geom_gene_label' else 'geom_gene_label_repel'
  legacy<-paste0('gene_',parameter)
  at<-which(base[['function']]==legacy_name & base$parameter==legacy)
 }
 if(length(at)!=1L) stop('Missing or ambiguous semantic review for ',name,'/',parameter)
 list(row=base[at,,drop=FALSE],path='api_arguments.csv',version='0.16.0',
  prior_name=legacy_name,prior_parameter=legacy)
}
parameters<-do.call(rbind,lapply(setdiff(exports,reexports$symbol),function(name) {
 args<-formals(getExportedValue('ggchord',name))
 do.call(rbind,lapply(seq_along(args),function(i) {
  parameter<-names(args)[i];review<-select_semantics(name,parameter);row<-review$row
  doc<-docs[[name]];text<-doc$args[[parameter]]
  if(is.null(text) || !nzchar(text)) stop('Current Rd argument missing: ',name,'/',parameter)
  responsibility<-grammar$responsibility[grammar$function_name==name & grammar$parameter==parameter]
  if(!length(responsibility)) responsibility<-if(family(name)=='output') 'result selection / target render / cache ownership' else 'explicit data workflow parameter'
  unit<-row$unit;auto_none<-row$auto_none_semantics
  null<-row$null_semantics;na<-row$na_semantics
  type<-row$type_contract;allowed<-row$allowed_values
  missing<-if(parameter=='...') 'May be empty; current forwarding/retired/error policy is entry-specific' else
   paste0('Omitted syntax default: ',default(args,i),'; contextual resolution follows current Rd/null_semantics')
  override<-''
  if(parameter=='shape' && name %in% c('geom_feature','geom_feature_plasmid')) {
   type<-'One character fixed glyph name; mapped feature_shape is separate'
   allowed<-'arrow/compact_arrow/promoter_arrow/primer_arrow/marker/block/capped_line/primer_arc/chevron/lollipop; explicit fixed shape overrides its mapping'
   null<-'Explicit NULL is not a fixed shape and is rejected; it is not the retired optional alias'
   na<-'Explicit NA is rejected as an invalid fixed shape'
   missing<-'Default arrow only when no mapped feature_shape exists; omitted fixed shape preserves mapping'
   override<-'v0.16.1 removed fixed feature_shape and made shape canonical; old alias semantics do not apply'
   auto_none<-'No auto/none fixed glyph; omitted shape retains a mapping; explicit NULL/NA are invalid fixed shapes'
  }
  if(name=='geom_link_ribbon' && parameter=='alpha') {
   allowed<-'Finite base alpha in [0,1]; a mapping does not bypass structural validation; explicit fixed alpha wins'
   override<-'v0.16.1 explicit fixed style priority supersedes old mapping-first wording'
  }
  note<-paste('Current Rd overrides renamed historical wording; family contract:',
    if(name %in% data[['function']]) 'API-DATA-CONTRACTS-v0.16.2.md' else
     if(family(name)=='output') 'API-ENGINE-CONTRACTS-v0.16.3.md' else 'API-GRAMMAR-v0.16.1.md')
  if(parameter=='shape' && name %in% c('geom_feature','geom_feature_plasmid'))
   note<-paste(note,'Fixed shape only; mapped feature_shape remains independent; explicit fixed style wins.')
  if(name %in% c('geom_gene_label_repel','geom_feature_label_repel') && parameter=='...')
   note<-paste(note,'v0.16.4 rejects known manual label controls; use the corresponding fixed-label layer.')
  data.frame(function_name=name,parameter=parameter,index=i,default_expression=default(args,i),
   responsibility=responsibility[1],type_contract=type,allowed_values=allowed,
   unit=unit,null_semantics=null,na_semantics=na,missing_semantics=missing,prior_missing_semantics=row$missing_semantics,
   auto_none_semantics=auto_none,prior_auto_none_semantics=row$auto_none_semantics,current_documentation=text,
   applies_to_coords=row$applies_to_coords,current_contract_note=note,current_semantic_override=override,
   current_decision='KEEP_CANDIDATE; no public rename/removal/default change in v0.16.4',
   semantic_source=review$path,semantic_version=review$version,
   prior_function=review$prior_name,prior_parameter=review$prior_parameter,
   prior_review_revision=row$review_source_revision %||% row$prior_review_revision %||% row$source_revision,
   current_rd=doc$file,source_revision=revision,
   runtime_evidence='Source consolidation; inherited limited evidence, not per-parameter execution')
 }))
}))
returns<-do.call(rbind,lapply(exports,function(name) {
 doc<-docs[[name]]
 data.frame(function_name=name,documented_return=doc$value %||% 'Upstream original object',
  documentation=doc$file %||% 'original dependency namespace documentation',
  stability='Retain documented members and meanings; additional/internal solver fields are not frozen',
  detailed_contract=if(family(name)=='output') 'engine-return-contracts-v0.16.3.csv; ggchord_layout_contract' else
    if(name %in% data[['function']]) 'data-return-contracts-v0.16.2.csv; ggchord_data_contract' else
      if(name %in% reexports$symbol) paste0(reexports$package[match(name,reexports$symbol)],' original object') else
       'Current entry Rd; API-GRAMMAR-v0.16.1.md',source_revision=revision)
}))
stopifnot(!anyDuplicated(paste(parameters$function_name,parameters$parameter)))
write.csv(manifest,file.path(out,'candidate-exports-v0.16.4.csv'),row.names=FALSE)
write.csv(parameters,file.path(out,'candidate-parameters-v0.16.4.csv'),row.names=FALSE)
write.csv(returns,file.path(out,'candidate-returns-v0.16.4.csv'),row.names=FALSE)
write.csv(data.frame(registration=grep('^S3method',readLines('NAMESPACE'),value=TRUE),
 policy='KEEP_DISPATCH; original generic protocol, not private implementation fields',source_revision=revision),
 file.path(out,'candidate-s3-v0.16.4.csv'),row.names=FALSE)
for(pair in list(c('grammar-compatibility-v0.16.1.csv','candidate-compatibility-v0.16.4.csv'),
 c('grammar-aesthetics-v0.16.1.csv','candidate-aesthetics-v0.16.4.csv'))) {
 value<-read.csv(file.path(out,pair[1]),na.strings=character(),check.names=FALSE)
 value$inherited_contract_source<-pair[1]
 if(pair[2]=='candidate-aesthetics-v0.16.4.csv') {
  extra<-data.frame(aesthetic_or_role=c('feature_label','track_value','track_fill','track_colour','group',
    'colour','fill','alpha','linewidth','linetype','size','family','fontface','lineheight',
    'seq_color','ribbon_color','link_color'),
   responsibility=c('input_biology','input_numeric_track_value',rep('scoped_visual_style',2),
    'track_series',rep('entry_specific_standard_aesthetic',9),rep('spelling_alias',3)),
   evaluation=c('Input-stage feature text; late stages forbidden',
    'Input numeric height; shared per-track domain, not colour',
    rep('Visual mapping on a scoped track; generated internal keys are opaque',2),
    'Track input series within each accession/window; never bridges NA or windows',
    rep('Only where the entry/geom supports it; standard ggplot2 semantics and fixed-value priority',9),
    'Alias of seq_colour; reject both spellings together',
    'Alias of ribbon_colour; reject both spellings together',
    'Alias of link_colour; reject both spellings together'),
   inherited_contract_source='Current R/geom-track.R, geom-feature.R and zzz.R; no new role added',
   stringsAsFactors=FALSE)
  value<-rbind(value,extra[!extra$aesthetic_or_role %in% value$aesthetic_or_role,])
 }
 value$candidate_status<-'KEEP; evidence level retained, not upgraded to runtime validation'
 write.csv(value,file.path(out,pair[2]),row.names=FALSE)
}
stopifnot(identical(before,tools::md5sum(inputs)))
write.csv(data.frame(path=inputs,md5=unname(before),source_revision=revision,
 source_state='Snapshot input bytes; revision is HEAD, hashes identify any working changes'),
 file.path(out,'candidate-source-hashes-v0.16.4.csv'),row.names=FALSE)
cat(length(exports),'exports;',nrow(parameters),'self-owned parameters;',nrow(returns),'returns consolidated.\n')
