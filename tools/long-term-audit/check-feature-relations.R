# Offline real plant observations plus explicitly synthetic schema controls.
pkgload::load_all(quiet = TRUE)
library(dplyr)
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(pattern, conditionMessage(error)))
}
# This companion four-column tool table is not GFF3. The caller explicitly
# supplies its positions as 1-based inclusive; the reader does not guess this.
positions <- function(path) {
  read.table(path, header = FALSE, col.names = c('accver', 'feature_id', 'start', 'end'),
    colClasses = c('character', 'character', 'numeric', 'numeric'), na.strings = character())
}
pair_path <- 'inst/extdata/plant-feature-pairs.homology'
group_path <- 'inst/extdata/plant-feature-groups.txt'
features <- positions('inst/extdata/plant-pair-features.gff')
mapped <- read_feature_relations(pair_path, feature_data = features, source_file = TRUE)
native <- read.table(pair_path, colClasses = 'character')
links <- mapped$ribbon_data
q <- features[match(native$V1, features$feature_id), ]
s <- features[match(native$V2, features$feature_id), ]
stopifnot(nrow(links) == 20L, nrow(mapped$membership) == 0L,
  identical(links$qfeature_id, native$V1), identical(links$sfeature_id, native$V2),
  identical(links$source_score, native$V3), identical(links$qaccver, q$accver),
  identical(links$saccver, s$accver), identical(links$qstart, q$start),
  identical(links$qend, q$end), identical(links$sstart, s$start), identical(links$send, s$end),
  all(is.na(links$direction)), !any(links$orientation_known),
  identical(links$.source_line, seq_len(20L)), all(links$.source_file == pair_path),
  identical(mapped$report$raw_files[[1L]], readLines(pair_path)),
  identical(mapped$report$feature_data, features), is.null(mapped$seq_data),
  !any(c('pident', 'evalue', 'support', 'length') %in% names(links)))
group_features <- positions('inst/extdata/plant-group-features.gff')
groups <- read_feature_relations(group_path, format = 'named_groups', feature_data = group_features)
tokens <- strsplit(readLines(group_path), '\t')
stopifnot(nrow(groups$ribbon_data) == 0L, nrow(groups$pairs) == 0L,
  nrow(groups$membership) == 12L, all(groups$membership$mapped),
  identical(groups$membership$group_id, rep(c('1166', '1167', '1168'), each = 4L)),
  identical(groups$membership$feature_id, unlist(lapply(tokens, function(x) x[-1L]), use.names = FALSE)))
unplaced <- read_feature_relations(pair_path)
stopifnot(nrow(unplaced$ribbon_data) == 0L, nrow(unplaced$pairs) == 20L,
  nrow(unplaced$report$unplaced_pairs) == 20L)
fail(read_feature_relations(pair_path, feature_data = group_features), 'unmatched')
partial <- read_feature_relations(pair_path, feature_data = features[-1L, ], unmatched = 'report')
stopifnot(nrow(partial$pairs) == 20L, nrow(partial$report$unplaced_pairs) == 1L,
  identical(partial$report$unmapped_ids, features$feature_id[1L]))
fail(read_feature_relations(pair_path, feature_data = rbind(features, features[1L, ])), 'duplicate')
# Format controls only: literal NA/numeric-looking IDs, physical comments,
# multiple files, negative annotation strand, duplicate members and singleton.
control <- tempfile(); second <- tempfile()
writeLines(c('# format only', '', '001 NA raw_score', 'NA 001'), control)
anchors <- data.frame(feature_id = c('001', 'NA'), accver = c('format_A', 'format_B'),
  start = c(2, 10), end = c(4, 20), strand = c('-', '+'))
stopifnot(file.copy(control, second))
result <- read_feature_relations(control, files = second, feature_data = anchors)
stopifnot(nrow(result$ribbon_data) == 4L, all(is.na(result$ribbon_data$direction)),
  identical(result$ribbon_data$qfeature_strand, rep(c('-', '+'), 2L)),
  identical(result$ribbon_data$qstart, rep(c(2, 10), 2L)),
  identical(result$pairs$.source_input, rep(1:2, each = 2L)),
  identical(result$pairs$.source_line, rep(3:4, 2L)), is.na(result$pairs$source_score[2L]))
writeLines(c('# format only', '001 001 NA', 'NA'), control)
result <- read_feature_relations(control, format = 'groups', feature_data = anchors)
stopifnot(nrow(result$membership) == 4L, all(is.na(result$membership$group_id)),
  identical(result$membership$.source_group, c(1L, 1L, 1L, 2L)), nrow(result$ribbon_data) == 0L)
writeLines('001 NA score surplus', control)
fail(read_feature_relations(control), '2 or 3 fields')
writeLines(c('# empty schema control', ''), control)
stopifnot(nrow(read_feature_relations(control)$pairs) == 0L)
unlink(c(control, second))
jsonlite::write_json(list(status = 'passed', sources = 'inst/extdata/feature-relation-sources.json',
  original_pairs = nrow(links), original_groups = length(tokens), original_members = 12L,
  exact_id_and_interval_mapping = TRUE, original_scores_and_records = TRUE,
  groups_never_expanded = TRUE, feature_strand_not_alignment_orientation = TRUE,
  no_source_length_identity_support_inference = TRUE,
  missing_ambiguous_and_namespace_controls = TRUE,
  coordinate_interpretation = 'Caller explicitly declares companion positions as 1-based-inclusive anchors; file suffix is not a schema detector.',
  scope = 'Offline input adaptation only; no homology inference or sequence analysis; no plot without independent source lengths.'),
  'tools/long-term-audit/v0.27.0-feature-relation-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('Feature relations: passed\n')
