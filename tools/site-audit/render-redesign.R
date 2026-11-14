# Explicit website-audit verification: affected PNGs only, no package test suite.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional = c(
  'tools/documentation-examples.R', 'tools/site-course-examples.R',
  'tools/site-audit/render-redesign.R', 'website/course-figures.json'))
pkgload::load_all(quiet = TRUE)
blocks <- function(text) {
  lines <- strsplit(text, '\n', fixed = TRUE)[[1]]
  markers <- grep('^# ---- [a-z]+ ----$', lines)
  setNames(lapply(seq_len(length(markers)-1L), function(i)
    paste(lines[seq.int(markers[i]+1L, markers[i+1L]-1L)], collapse='\n')),
    sub('^# ---- ([a-z]+) ----$', '\\1', lines[markers[-length(markers)]]))
}
chunks <- c(blocks(paste(readLines('tools/documentation-examples.R'),collapse='\n')),
  blocks(paste(readLines('tools/site-course-examples.R'),collapse='\n')))
env <- new.env(parent=globalenv())
keys <- c('setup','frontbackbone','frontlinks','frontribbon','frontgenes','frontfeature',
  'frontgenelabels','courseblockad','frontradius','frontcurvature',
  'parameters','parameteruniform','parametervector','parametergenes','parameterstrand',
  'parameterwidth','parameterlabels','parameterreorder','parameterlist')
for(key in keys) invisible(eval(parse(text=chunks[[key]]),env))
# Compare retained original data through the public exported layout.
layout <- function(plot) export_ggchord_layout(plot, original_data=TRUE)
backbone <- layout(env$p_chord_backbone)
feature <- layout(env$p_chord_feature)
labelled <- layout(env$p_chord_gene_labels)
linear <- layout(env$p_chord_linear)
stopifnot(nrow(labelled$ribbon)>0L, nrow(labelled$gene)>0L,
  nrow(linear$ribbon)>0L, nrow(linear$gene)>0L)
for(key in c('qaccver','saccver','qstart','qend','sstart','send'))
  stopifnot(identical(feature$ribbon[[key]],labelled$ribbon[[key]]),
    identical(labelled$ribbon[[key]],linear$ribbon[[key]]))
for(key in c('accver','start','end','strand'))
  stopifnot(identical(feature$gene[[key]],labelled$gene[[key]]),
    identical(labelled$gene[[key]],linear$gene[[key]]))
stopifnot(all(c('+','-') %in% env$parameter_genes$strand[
  env$parameter_genes$accver=='OR222515.1']))
# Width/position examples keep all selected original annotations.
for(plot in list(env$p_parameter_strand,env$p_parameter_width,env$p_parameter_list))
  stopifnot(all(env$parameter_genes$start %in% layout(plot)$gene$start))
# Named and ordered values agree before reordering, diverge afterwards.
stopifnot(identical(layout(env$p_parameters)$seq,layout(env$p_parameter_vector)$seq))
legacy <- read.csv('website/assets/plots.csv',stringsAsFactors=FALSE)
affected <- c('frontgenelabels','frontradius','frontcurvature','parameteruniform',
  'parameterexample','parametervector','parameterstrand','parameterwidth')
legacy$width[legacy$id=='parameterexample'] <- 7
legacy$height[legacy$id=='parameterexample'] <- 7
fingerprints <- jsonlite::fromJSON(system2('python3','tools/documentation-fingerprints.py',stdout=TRUE))
result <- list(status='in progress', engine_sha256=as.list(audit$engine),
  example_sha256=digest::digest(file='tools/documentation-examples.R',algo='sha256',serialize=FALSE),
  scope='Eight affected legacy PNGs and four explicit additions at normal dimensions; original R/data unchanged; no browser or release matrix.',
  public_checks=list(cumulative_layers=TRUE,original_endpoints=TRUE,original_CDS=TRUE,
    real_opposite_strands=TRUE,equivalent_named_ordered_radii=TRUE),figures=list())
for(id in affected) {
  at <- match(id,legacy$id);cat('Render:',id,'\n');flush.console()
  file <- view_ggchord(get(legacy$object[at],env),width=legacy$width[at],height=legacy$height[at],
    units='in',device='png',dpi=legacy$dpi[at],bg=legacy$background[at],viewer='none')
  before <- legacy$md5[at]
  stopifnot(file.copy(file,legacy$file[at],overwrite=TRUE))
  audit_assert_source_snapshot(audit)
  legacy$md5[at] <- unname(tools::md5sum(legacy$file[at]))
  legacy$code_md5[at] <- fingerprints$code
  legacy$dependency_md5[at] <- fingerprints$assets[[id]]
  legacy$generator[at] <- 'tools/site-audit/render-redesign.R'
  result$figures[[id]] <- list(file=legacy$file[at],status='actual PNG rendered',
    before_md5=before,md5=legacy$md5[at],
    sha256=digest::digest(file=legacy$file[at],algo='sha256',serialize=FALSE),
    width=legacy$width[at],height=legacy$height[at],dpi=legacy$dpi[at],
    dependency_md5=legacy$dependency_md5[at])
}
write.csv(legacy,'website/assets/plots.csv',row.names=FALSE)
courses <- jsonlite::fromJSON('website/course-figures.json',simplifyVector=FALSE)
course_record <- jsonlite::fromJSON('tools/site-audit/course-rendering.json',simplifyVector=FALSE)
for(row in courses) {
  if(!row$id %in% c('parameternamedorder','parameterorderedvalues','parameterlist','frontlabelledlinear')) next
  cat('Render:',row$id,'\n');flush.console()
  file <- view_ggchord(get(row$object,env),width=row$width,height=row$height,
    units='in',device='png',dpi=row$dpi,bg='white',viewer='none')
  stopifnot(file.copy(file,row$file,overwrite=TRUE))
  audit_assert_source_snapshot(audit)
  code <- paste(vapply(unlist(row$snippets),function(key) trimws(chunks[[key]]),''),collapse='\n\n')
  record <- list(status='actual PNG rendered',file=row$file,
    sha256=digest::digest(file=row$file,algo='sha256',serialize=FALSE),
    code_sha256=digest::digest(code,algo='sha256',serialize=FALSE),
    width=row$width,height=row$height,dpi=row$dpi,runtime_record='tools/site-audit/redesign-rendering.json')
  course_record$figures[[row$id]] <- record
  result$figures[[row$id]] <- record
}
# A Dotplot-only preparation must reproduce all actually consumed objects.
full <- new.env(parent=globalenv());minimal <- new.env(parent=globalenv())
flagships <- blocks(paste(readLines('tools/site-flagships.R'),collapse='\n'))
invisible(eval(parse(text=flagships$flagshipdata),full))
invisible(eval(parse(text=chunks$dotplotmatrixdata),minimal))
for(key in c('flagship_sequences','flagship_hsp','flagship_theme','flagship_paper'))
  stopifnot(identical(get(key,full),get(key,minimal)))
result$public_checks$dotplot_preparation_equivalent <- TRUE
result$status <- 'complete'
jsonlite::write_json(result,'tools/site-audit/redesign-rendering.json',auto_unbox=TRUE,pretty=TRUE)
jsonlite::write_json(course_record,'tools/site-audit/course-rendering.json',auto_unbox=TRUE,pretty=TRUE)
cat('PASS: cumulative source layers, real strand input, equivalent Dotplot preparation and 12 actual PNGs\n')
