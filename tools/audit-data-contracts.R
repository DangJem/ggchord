# Rebuild v0.16.2 data evidence without changing the v0.16.0 audit snapshot.
pkgload::load_all(quiet = TRUE)
functions <- c('read_bed','read_paf','read_blast','read_gff3','read_fasta_lengths','read_genbank',
  'as_ggchord_seq','as_ggchord_feature','as_ggchord_link','validate_ggchord_data',
  'validate_ggchord_hierarchy','clean_ggchord_data','prepare_ggchord_data',
  'find_ggchord_loci','focus_ggchord_data','focus_ggchord_track','segment_ggchord_data',
  'map_ggchord_position','build_ggchord_cds','project_ggchord_features','project_ggchord_links',
  'filter_ggchord_ribbons','deduplicate_ggchord_ribbons','merge_ggchord_ribbons',
  'bundle_ggchord_ribbons','optimize_ggchord_layout','stat_ribbon_bundle','stat_ribbon_density',
  'link_ggchord_clusters','find_common_features','primer_catalog','find_primer_bindings',
  'find_restriction_sites','filter_restriction_sites')
inputs <- list.files('design-roadmap/audit-reviewed', '[.]csv$', full.names = TRUE)
files <- sort(unique(c(inputs, list.files('R', '[.]R$', full.names = TRUE),
  list.files('man', '[.]Rd$', full.names = TRUE), 'NAMESPACE','DESCRIPTION',
  'tools/audit-data-contracts.R','tools/validate-data-contracts.R','tools/validate-genbank.R',
  'NEWS.md','design-roadmap/API-DATA-CONTRACTS-v0.16.2.md')))
hashes <- tools::md5sum(files)
reviewed <- do.call(rbind, lapply(inputs, function(path) {
  x <- read.csv(path, na.strings = character(), check.names = FALSE)
  x <- x[x[['function']] %in% functions, , drop = FALSE]
  x$prior_review_path <- rep(path, nrow(x))
  x
}))
key <- paste(reviewed[['function']], reviewed$parameter)
stopifnot(!anyDuplicated(key))
parameters <- do.call(rbind, lapply(functions, function(fn) {
  args <- formals(getExportedValue('ggchord', fn))
  do.call(rbind, lapply(seq_along(args), function(i) {
    parameter <- names(args)[i]
    index <- match(paste(fn, parameter), key)
    if (is.na(index)) stop('Missing semantic review: ', fn, '/', parameter)
    row <- reviewed[index, , drop = FALSE]
    row$current_default <- if (parameter == '...') '<ellipsis>' else
      paste(deparse(args[[i]], width.cutoff = 400), collapse = ' ')
    if (!nzchar(row$current_default)) row$current_default <- '<required>'
    row$prior_review_revision <- row$source_revision
    row$prior_rationale <- row$rationale
    row$source_revision <- 'v0.16.2 working sources; data-source-hashes-v0.16.2.csv'
    row$decision <- 'KEEP_PUBLIC_ENTRY; CORRECTNESS_FIX_OR_DOCUMENTED_BOUNDARY'
    row$rationale <- 'Current contract retained with documented fixes/boundaries; see API-DATA-CONTRACTS-v0.16.2.md. Historical rationale is separate.'
    row$current_review <- 'SOURCE_REVIEWED_v0.16.2; see API-DATA-CONTRACTS-v0.16.2.md'
    row$runtime_evidence <- 'LIMITED: completion record; not a per-parameter execution claim'
    if (parameter %in% c('file','files')) {
      row$type_contract <- 'NULL or character path vector / glob patterns; nonmissing and nonempty'
      row$actual_usage <- 'Unified expansion; file then files; normalized-path deduplication; first spelling/order retained; UTF-8 input; source lines or record/feature IDs and import report retained'
      row$rationale <- 'F16 KEEP dual entries; FIX asymmetric expansion, identifier parsing and source provenance'
    }
    if (fn == 'validate_ggchord_data' && parameter == 'strict')
      row$actual_usage <- 'Collect structural/reference/bounds issues; strict TRUE stops after collection if errors exist; FALSE returns full report'
    if (fn == 'deduplicate_ggchord_ribbons' && parameter %in% c('by','min_reciprocal_overlap'))
      row$actual_usage <- 'Reciprocal intersection covers each interval; normalize descending bounds for overlap; retain distinct relative endpoint directions; report removed representative mapping'
    if (fn == 'merge_ggchord_ribbons')
      row$actual_usage <- paste(row$actual_usage, 'v0.16.2: pair IDs always group; singleton length unchanged; merged .length_estimated; report source_length/source_pident lists')
    if (fn == 'find_common_features')
      row$actual_usage <- paste(row$actual_usage, 'v0.16.2 explicit boundary: standard-code linear six-frame only; source geneticCode not interpreted; confidence heuristic not probability')
    if (fn %in% c('clean_ggchord_data','prepare_ggchord_data','deduplicate_ggchord_ribbons',
        'merge_ggchord_ribbons','bundle_ggchord_ribbons','stat_ribbon_bundle','stat_ribbon_density'))
      row$actual_usage <- paste(row$actual_usage, 'F47: zero endpoint order is unknown unless supplied direction evidence exists; clean retains pre-clip evidence; single-base endpoints stay independent in merge/bundle')
    row
  }))
}))
stopifnot(identical(hashes, tools::md5sum(files)))
out <- 'design-roadmap/audit-generated'
write.csv(parameters, file.path(out,'data-parameters-v0.16.2.csv'), row.names = FALSE)
contract_lines <- readLines('design-roadmap/API-DATA-CONTRACTS-v0.16.2.md', encoding = 'UTF-8')
contract_lines <- grep('^\\| `[a-z0-9_]+` \\|', contract_lines, value = TRUE)
contracts <- do.call(rbind, lapply(contract_lines, function(line) {
  parts <- trimws(strsplit(line, '|', fixed = TRUE)[[1L]])
  data.frame(function_name = gsub('`', '', parts[2L], fixed = TRUE),
    input_null_na = parts[3L], return_provenance_errors = parts[4L],
    evidence = 'SOURCE_CONTRACT; finite runtime scope in completion record',
    stringsAsFactors = FALSE)
}))
stopifnot(setequal(contracts$function_name, functions), !anyDuplicated(contracts$function_name))
write.csv(contracts, file.path(out,'data-return-contracts-v0.16.2.csv'), row.names = FALSE)
write.csv(data.frame(path = names(hashes), md5 = unname(hashes)),
  file.path(out,'data-source-hashes-v0.16.2.csv'), row.names = FALSE)
cat(length(functions), 'public data entries;', nrow(parameters), 'current parameters; source hashes saved.\n')
