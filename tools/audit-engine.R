# Rebuild v0.16.3 structural evidence, not runtime/release validation.
pkgload::load_all(quiet=TRUE)
ns <- asNamespace('ggchord')
outdir <- 'design-roadmap/audit-generated'
revision <- system2('git',c('rev-parse','HEAD'),stdout=TRUE)
files <- sort(list.files('R',pattern='\\.R$',full.names=TRUE))
parsed <- lapply(files,parse,keep.source=TRUE)
names(parsed) <- files
fingerprint <- tools::md5sum(files)
exported <- getNamespaceExports('ggchord')
head_name <- function(x) {
  if (is.symbol(x)) return(as.character(x))
  if (is.call(x) && as.character(x[[1]]) %in% c('::',':::') && length(x)==3L)
    return(as.character(x[[3]]))
  ''
}
definitions <- list();calls <- list();symbols <- character();strings <- character()
walk <- function(x,caller,file,quoted=FALSE) {
  if(is.symbol(x)) {symbols <<- c(symbols,as.character(x));return(invisible(NULL))}
  if(is.character(x)) {strings <<- c(strings,x);return(invisible(NULL))}
  if(!is.call(x) && !is.expression(x) && !is.pairlist(x)) return(invisible(NULL))
  if(is.call(x)) {
    callee <- head_name(x[[1]])
    if(nzchar(callee)) calls[[length(calls)+1L]] <<- data.frame(
      caller=caller,callee=callee,source_file=file,quoted=quoted,
      evidence='SYNTAX_ONLY; calls/callbacks are not a complete runtime graph')
    quoted <- quoted || callee %in% c('quote','substitute','expression','bquote')
  }
  for (value in as.list(x)) if(!rlang::is_missing(value)) walk(value,caller,file,quoted)
}
for(file in files) for(index in seq_along(parsed[[file]])) {
  expr <- parsed[[file]][[index]]
  caller <- '<top-level>'
  if(is.call(expr) && head_name(expr[[1]]) %in% c('<-','=') && is.symbol(expr[[2]])) {
    caller <- as.character(expr[[2]])
    if(is.call(expr[[3]]) && identical(head_name(expr[[3]][[1]]),'function')) {
      sr <- attr(parsed[[file]],'srcref')[[index]]
      definitions[[length(definitions)+1L]] <- data.frame(name=caller,source_file=file,
        line=if(length(sr)) as.integer(sr[1]) else NA_integer_)
    }
  }
  walk(expr,caller,file)
}
definitions <- do.call(rbind,definitions)
calls <- do.call(rbind,calls)
dependencies <- unique(calls[calls$callee %in% definitions$name,])
dependencies$source_revision <- revision
write.csv(dependencies,file.path(outdir,'engine-internal-dependencies-v0.16.3.csv'),row.names=FALSE)

# Stage signatures and literal result names are executable contracts, not text guesses.
short <- c('sequence','axis','ribbon','annotation','transform','label','finalize')
functions <- paste0('ggchord_layout_',short,'_step')
producers <- character();stage_rows<-list();field_rows<-list()
for(i in seq_along(functions)) {
  f <- get(functions[i],ns);input<-names(formals(f));tail<-utils::tail(as.list(body(f)),1L)[[1L]]
  stopifnot(is.call(tail),identical(head_name(tail[[1]]),'list'))
  output <- names(as.list(tail))[-1L]
  upstream <- unique(unname(producers[intersect(input,names(producers))]))
  own <- unique(definitions$source_file[definitions$name==functions[i]])
  stage_rows[[i]]<-data.frame(stage=short[i],source_file=own,inputs=paste(input,collapse=';'),
    outputs=paste(output,collapse=';'),depends_on=paste(upstream,collapse=';'),
    owner='one compute_chord_layout call; reference cache belongs to outer geometry call',
    baseline='723575a',source_revision=revision,
    evidence='function formals/results plus manual read-before-write/closure/data-mask review')
  field_rows[[i]]<-data.frame(stage=short[i],field=input,
    producer=ifelse(input %in% names(producers),unname(producers[input]),'normalized_input'),
    invalidation='new layout computation; device/context changes invalidate estimates')
  producers[output] <- short[i]
}
write.csv(do.call(rbind,stage_rows),file.path(outdir,'engine-stages-v0.16.3.csv'),row.names=FALSE)
write.csv(do.call(rbind,field_rows),file.path(outdir,'engine-stage-fields-v0.16.3.csv'),row.names=FALSE)

# Low-use candidates require positive protocol/string/callback evidence or manual review.
namespace <- readLines('NAMESPACE')
s3 <- sub('^S3method\\(([^,]+),([^)]*)\\).*','\\1.\\2',namespace[grepl('^S3method\\(',namespace)])
s3 <- gsub('"', '', s3, fixed=TRUE)
counts <- table(symbols);direct <- table(calls$callee)
candidates <- definitions[!definitions$name %in% exported & vapply(definitions$name,
  function(x) sum(calls$callee==x)==0L || x %in% c('if_null_else','breakPointsFunc'),logical(1)),]
candidates$direct_calls <- as.integer(direct[candidates$name]);candidates$direct_calls[is.na(candidates$direct_calls)]<-0L
candidates$symbol_uses <- as.integer(counts[candidates$name])-1L
candidates$string_mentions <- vapply(candidates$name,function(x) sum(grepl(x,strings,fixed=TRUE)),integer(1))
candidates$s3_registration <- candidates$name %in% s3
candidates$quoted_mentions <- vapply(candidates$name,function(x) any(dependencies$callee==x & dependencies$quoted),logical(1))
candidates$decision <- ifelse(candidates$s3_registration,'KEEP_S3',
  ifelse(candidates$symbol_uses>0L | candidates$string_mentions>0L,'KEEP_REFERENCED','REVIEW_REQUIRED'))
candidates$rationale <- ifelse(candidates$s3_registration,'NAMESPACE dispatch protocol',
  'Positive symbol-value references reviewed with ggproto/callback/string/quoted syntax; no deletion inferred')
# Reviewed exceptions from source paths. No zero-count auto deletion.
for(i in seq_len(nrow(candidates))) {
 name<-candidates$name[i]
 if(name=='.onLoad') {candidates$decision[i]<-'KEEP_LIFECYCLE';candidates$rationale[i]<-'R namespace lifecycle hook registers theme elements'}
 if(name=='ggchord_parse_rebase') {candidates$decision[i]<-'KEEP_BUILD_TIME';candidates$rationale[i]<-'data-raw/generate_rebase_database.R calls this parser to rebuild bundled source data'}
 if(name=='ggchord_label_box_conflicts') {candidates$decision[i]<-'KEEP_VALIDATION';candidates$rationale[i]<-'Used by tools/validate-release.R and the existing radial-geometry regression; internal detector is not a public field contract'}
 if(name=='if_null_else') candidates$rationale[i]<-'Used by draw-key callback in helpers_elements.R; repeated NULL semantics are not dead code'
 if(name=='breakPointsFunc') candidates$rationale[i]<-'Axis calls consume it; unused n/tol retain the documented internal signature; no numerical rewrite'
 if(name %in% functions) candidates$rationale[i]<-'Passed as function values in the ordered stage list; do.call dispatches declared arguments'
}
review_file<-'design-roadmap/audit-reviewed-details/engine-unused-review-v0.16.3.csv'
if(file.exists(review_file)) {
 review<-read.csv(review_file,stringsAsFactors=FALSE,check.names=FALSE,na.strings=character())
 for(i in seq_len(nrow(candidates))) {
  at<-match(candidates$name[i],review$name)
  if(!is.na(at)) {candidates$decision[i]<-review$decision[at];candidates$rationale[i]<-review$rationale[at]}
 }
 removed<-review[review$decision=='REMOVE',]
 for(i in seq_len(nrow(removed))) {
  stopifnot(!removed$name[i] %in% definitions$name)
  candidates<-rbind(candidates,data.frame(name=removed$name[i],source_file=removed$source_file[i],line=NA_integer_,
   direct_calls=0L,symbol_uses=0L,string_mentions=0L,s3_registration=FALSE,quoted_mentions=FALSE,
   decision='REMOVE',rationale=removed$rationale[i]))
 }
}
candidates$source_revision<-revision
write.csv(candidates,file.path(outdir,'engine-unused-candidates-v0.16.3.csv'),row.names=FALSE)

contract_file <- 'design-roadmap/audit-reviewed-details/engine-output-contracts-v0.16.3.csv'
if(file.exists(contract_file)) {
  returns<-read.csv(contract_file,stringsAsFactors=FALSE,na.strings=character(),check.names=FALSE)
  stopifnot(all(returns$function_name %in% exported))
  returns$source_revision<-revision
  write.csv(returns,file.path(outdir,'engine-return-contracts-v0.16.3.csv'),row.names=FALSE)
}
# Preserve historical semantic review separately from current formals/ownership.
public <- c('get_chord_layout','export_ggchord_layout','view_ggchord','diagnose_ggchord')
prior <- read.csv('design-roadmap/audit-reviewed/output-diagnostics.csv',
  check.names=FALSE,na.strings=character(),stringsAsFactors=FALSE)
parameters<-do.call(rbind,lapply(public,function(name) {
 args<-formals(get(name,ns))
 do.call(rbind,lapply(seq_along(args),function(i) {
  at<-which(prior[['function']]==name & prior$parameter==names(args)[i])
  stopifnot(length(at)==1L)
  row<-prior[at,,drop=FALSE]
  row$prior_review_revision<-row$source_revision
  row$current_default<-paste(deparse(args[[i]],width.cutoff=400),collapse=' ')
  if(!nzchar(row$current_default)) row$current_default<-'<required>'
  row$current_ownership<-switch(name,
   get_chord_layout='Plot-owned cache; all layouts invalidate on device/context change',
   export_ggchord_layout='Consistent local layer registration; getter may fill plot cache',
   view_ggchord='Private render copy preserves source layer parameters and cache',
   diagnose_ggchord='Private render copy and restored temporary context/device')
  row$current_review<-'SOURCE_REVIEWED_v0.16.3; documented fields retained; see ggchord_layout_contract'
  row$runtime_evidence<-'LIMITED representative cases in engine state validation; not every parameter executed'
  row$source_revision<-revision
  row
 }))
}))
write.csv(parameters,file.path(outdir,'engine-parameters-v0.16.3.csv'),row.names=FALSE)
stopifnot(identical(fingerprint,tools::md5sum(files)))
tracked_inputs<-sort(unique(c(files,'NAMESPACE','DESCRIPTION','NEWS.md','tools/audit-engine.R',
 'tools/validate-engine.R','tools/validate-engine-state.R','tools/validate-engine-registry.R',
 'tools/validate-tracks.R',
 'design-roadmap/audit-reviewed/output-diagnostics.csv',
 'design-roadmap/audit-reviewed-details/engine-output-contracts-v0.16.3.csv',
 'design-roadmap/audit-reviewed-details/engine-unused-review-v0.16.3.csv',
 'man/ggchord_layout_contract.Rd','man/get_chord_layout.Rd',
 'man/export_ggchord_layout.Rd','man/view_ggchord.Rd','man/diagnose_ggchord.Rd')))
hashes <- data.frame(path=tracked_inputs,md5=unname(tools::md5sum(tracked_inputs)),source_revision=revision)
write.csv(hashes,file.path(outdir,'engine-source-hashes-v0.16.3.csv'),row.names=FALSE)
cat(nrow(dependencies),'syntax dependencies;',length(functions),'explicit stages;',
  nrow(candidates),'reviewed low-use candidates;',nrow(parameters),'public output parameters.\n')
print(candidates[c('name','direct_calls','symbol_uses','string_mentions','s3_registration','decision')],row.names=FALSE)
