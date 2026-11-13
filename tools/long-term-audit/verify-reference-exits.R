# Verify completed direct renders and exact source inputs; no plot rebuild or release checks.
source('tools/audit-source-snapshot.R');audit<-audit_source_snapshot()
pkgload::load_all(quiet=TRUE)
image_path <- function(path) {
  if (file.exists(path)) return(path)
  retention <- jsonlite::read_json('tools/long-term-audit/v0.27.0-exit-image-retention.json')
  reuse <- retention$reused_images[[basename(path)]]
  stopifnot(!is.null(reuse))
  file.path('tools/long-term-audit',reuse$file)
}
before<-jsonlite::read_json('tools/long-term-audit/v0.27.0-compact-final-results.json')
report<-jsonlite::read_json('tools/long-term-audit/v0.27.0-exit-final-results.json')
stopifnot(identical(report$runtime_sha256,as.list(audit$engine)))
cases<-report$cases
stopifnot(!anyDuplicated(names(cases)),setequal(names(cases),names(before$cases)))
cases<-cases[names(before$cases)];summary<-list()
for(id in names(cases)) {
 suffix<-if(grepl('exit-final-a-',cases[[id]]$file,fixed=TRUE))'a' else 'b'
 a<-readRDS(paste0('/tmp/v0.27.0-compact-final-',id,'.rds'))
 b<-readRDS(paste0('/tmp/v0.27.0-exit-final-',suffix,'-',id,'.rds'))
 stopifnot(identical(a$layer_inputs,b$layer_inputs))
 v<-cases[[id]];q<-v$quality
 stopifnot(isTRUE(v$converged),!length(v$warnings),all(vapply(v$visibility,function(z)isTRUE(z$visible),logical(1))),
  q$n_text_overlaps==0L,q$n_leader_crossings==0L,q$n_leaders_through_text==0L,q$n_label_outside_canvas==0L,
  v$association$anchors_off_own_feature==0L,v$association$leaders_over_one_bend==0L,
  v$association$internal_glyphs_off_tangent==0L,v$association$inside_labels_outside_own_feature==0L,
  all(vapply(Filter(function(z)isTRUE(z$endpoint),v$axis_labels),function(z)isTRUE(z$visible),logical(1))))
 summary[[id]]<-list(identical_complete_layer_inputs=TRUE,
  source_labels=length(v$visibility),through_text_before=before$cases[[id]]$quality$n_leaders_through_text,
  through_text_after=q$n_leaders_through_text,converged=TRUE,file=v$file,
  png_sha256=digest::digest(file=image_path(v$file),algo='sha256'))
}
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(status='passed measured controls for all 13 real cases; not release or arbitrary-input certification',
 cases=summary,source_labels=sum(vapply(summary,`[[`,integer(1),'source_labels')),
 runtime_sha256=as.list(audit$engine),remaining='Separate matched/compact flagship contacts and small-canvas failures; stroke/ink/portability and final documentation still need their own evidence'),
 'tools/long-term-audit/v0.27.0-exit-final-verification.json',pretty=TRUE,auto_unbox=TRUE)
cat('All 13 unchanged complete inputs, source labels and measured simple exits: passed\n')
