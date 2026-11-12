# Existing AA annotations and explicit format controls; no annotation engine.
pkgload::load_all(quiet = TRUE)
library(dplyr)
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(pattern, conditionMessage(error)))
}
for (accession in c('P00517', 'P51587')) {
  json_path <- paste0('inst/extdata/', accession, '-pfam.json')
  tsv_path <- paste0('inst/extdata/', accession, '-signatures.tsv')
  json <- read_signature_matches(json_path, format = 'entry_json', source_file = TRUE)
  tsv <- read_signature_matches(tsv_path, source_file = TRUE)
  a <- json$features$signature; b <- tsv$features$signature
  # These particular real caches contain one fragment per native location.
  stopifnot(identical(a[c('accver', 'start', 'end', 'signature_id', 'source_length')],
                      b[c('accver', 'start', 'end', 'signature_id', 'source_length')]),
    all(is.na(a$feature_id)), all(is.na(b$feature_id)),
    all(vapply(a$source_score, is.numeric, logical(1))),
    all(vapply(b$source_score, is.character, logical(1))),
    all(a$.source_file == json_path), all(b$.source_file == tsv_path),
    identical(json$report$raw_files[[1L]], readLines(json_path, warn = FALSE)),
    !any(c('pident', 'evalue', 'confidence', 'read_support', 'sequence') %in% names(a)))
  source <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  row <- 0L
  for (entry in source$results) for (protein in entry$proteins) for (location in protein$entry_protein_locations) {
    row <- row+1L
    stopifnot(a$accver[row] == protein$accession,
      a$source_length[row] == protein$protein_length,
      identical(a$source_score[[row]], location$score),
      identical(a$representative[row], location$representative),
      nrow(a$segments[[row]]) == length(location$fragments))
  }
  stopifnot(row == nrow(a), all(json$seq_data$sequence_unit == 'aa'))
}
# Explicit synthetic discontinuity/pagination controls, no biological claim.
fixture <- tempfile(fileext = '.json')
control <- list(count = 2L, `next` = 'https://example.invalid/not-followed', previous = NULL,
  results = list(list(metadata = list(accession = 'FORMAT_MODEL', name = 'Format control',
    source_database = 'FORMAT_LIBRARY', type = 'domain'),
    proteins = list(list(accession = 'FORMAT_ONLY', protein_length = 100L,
      entry_protein_locations = list(list(model = 'FORMAT_MODEL', score = NULL,
        representative = FALSE, fragments = list(
          list(start = 10L, end = 20L, `dc-status` = 'N_TERMINAL_DISC'),
          list(start = 40L, end = 50L, `dc-status` = 'C_TERMINAL_DISC')))))))))
jsonlite::write_json(control, fixture, auto_unbox = TRUE, null = "null")
segmented <- read_signature_matches(fixture, format = 'entry_json')
feature <- segmented$features$signature
stopifnot(nrow(feature) == 1L, nrow(feature$segments[[1L]]) == 2L,
  !feature$representative, is.null(feature$source_score[[1L]]),
  segmented$report$records[[1L]]$pagination$next_page == 'https://example.invalid/not-followed')
# Existing geometry creates separate visible runs, retaining one source match.
normalized <- ggchord_feature_data(feature, fixed_shape = 'block')
stopifnot(length(unique(normalized$.biological_source_row)) == 1L,
  nrow(normalized$.feature_runs[[1L]]) == 2L)
unlink(fixture)
# Real repeat matches stay distinct; model IDs are not feature-instance IDs.
observed <- read_signature_matches('inst/extdata/P51587-pfam.json', format = 'entry_json')
features <- observed$features$signature
stopifnot(nrow(features) == 13L, anyDuplicated(features$signature_id) > 0L,
  all(is.na(features$feature_id)), observed$seq_data$length == 3418)
plot <- ggchord(observed$seq_data) + geom_seq(show.legend = FALSE) + geom_seq_label() +
  geom_feature(aes(feature_fill = signature_id), data = features, shape = 'block',
    position = position_feature_stack()) + coord_collinear() +
  labs(title = 'Existing protein signature locations', subtitle = 'Declared AA coordinates; repeat matches remain distinct')
layout <- export_ggchord_layout(plot, original_data = TRUE)
ggplot2::ggsave('tools/long-term-audit/v0.27.0-signature-features.png', plot,
  width = 10, height = 4.5, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', sources = 'inst/extdata/signature-sources.json',
  original_accession_case = TRUE, original_full_aa_lengths = TRUE,
  real_locations = nrow(features), native_scores_and_flags_preserved = TRUE,
  json_and_explicit_tsv_reencoding = TRUE, one_match_with_native_segments = TRUE,
  no_gap_fill_or_pagination_fetch = TRUE, no_annotation_engine_run = TRUE,
  no_inferred_score_metric_identity_or_feature_instance_id = TRUE),
  'tools/long-term-audit/v0.27.0-signature-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('Signature adapters: passed\n')
