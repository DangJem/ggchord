# Direct adapters/real-data check. Requires optional rtracklayer for binary paths.
# No package tests, build/check, vignette compilation or release device matrix.
pkgload::load_all(quiet = TRUE)
library(dplyr)
if (!requireNamespace('rtracklayer', quietly = TRUE)) stop('This check requires optional rtracklayer')
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
text_path <- system.file('extdata', 'pUC19-gc-100.bedGraph', package = 'ggchord')
binary_path <- system.file('extdata', 'pUC19-gc-100.bw', package = 'ggchord')
text <- read_bedgraph(text_path, source_file = TRUE)
binary <- read_bigwig(binary_path, source_file = TRUE)
data(plasmid_example_pUC19)
sequence <- plasmid_example_pUC19$sequence[[1L]]
actual <- vapply(seq_len(nrow(text)), function(i) {
  bases <- strsplit(substr(sequence, text$start[i], text$end[i]), '', fixed = TRUE)[[1L]]
  100*mean(bases %in% c('G', 'C'))
}, numeric(1))
stopifnot(nrow(text) == 27, all(abs(text$value-actual) < 1e-10),
  text$start[1] == 1, tail(text$end, 1) == nchar(sequence),
  identical(text$start, text$chromStart+1), identical(text$end, text$chromEnd),
  all(text$.source_line == seq.int(4L, nrow(text)+3L)),
  all(text$.source_file == text_path))
# Compare the represented position-value signal, not original-window counts:
# a binary format may encode runs rather than the source table's observations.
expand_signal <- function(d) {
  out <- rep(NA_real_, nchar(sequence))
  for (i in seq_len(nrow(d))) out[seq.int(d$start[i], d$end[i])] <- d$value[i]
  out
}
stopifnot(max(abs(expand_signal(text)-expand_signal(binary))) < 1e-4,
  all(binary$value == binary$score))
metadata <- attr(binary, 'ggchord_import')
stopifnot(metadata$sequence_info[[1L]]$length == nchar(sequence),
  metadata$sequence_info[[1L]]$accver == 'pUC19', metadata$decoder$package == 'rtracklayer')
# Real partial windows retain decoder clipping and union provenance.
query <- tibble(accver = 'pUC19', start = c(50, 100), end = c(150, 250))
selected <- read_bigwig(binary_path, which = query)
stopifnot(min(selected$start) == 50, max(selected$end) == 250,
  sum(selected$end-selected$start+1) == 201,
  identical(attr(selected, 'ggchord_import')$query_input, query),
  nrow(attr(selected, 'ggchord_import')$selected_intervals[[1L]]) == 1)
empty <- read_bigwig(binary_path, which = query[0, ])
stopifnot(nrow(empty) == 0, identical(names(empty), c('accver', 'start', 'end', 'value', 'score')))
fail(read_bigwig(binary_path, which = data.frame(accver = 'wrong', start = 1, end = 10)), 'absent')
fail(read_bigwig(binary_path, which = data.frame(accver = 'pUC19', start = 1, end = 3000)), 'exceeds')
# Explicitly synthetic format-only boundaries. No biological measurement claim.
synthetic <- GenomicRanges::GRanges('track', IRanges::IRanges(c(1, 21, 41), c(10, 30, 50)),
  score = c(-2, 4, 0), seqlengths = c(track = 100))
fixture <- tempfile(fileext = '.bw')
rtracklayer::export(synthetic, fixture, format = 'BigWig')
gap <- read_bigwig(fixture, which = data.frame(accver = 'track', start = 11, end = 20))
zero <- read_bigwig(fixture, which = data.frame(accver = 'track', start = 41, end = 50))
stopifnot(nrow(gap) == 0, nrow(zero) == 1, zero$value == 0)
unlink(fixture)
# Header keywords can be literal IDs; actual directives remain provenance.
path <- tempfile(fileext = '.bedGraph')
writeLines(c('track type=bedGraph name="format-only fixture"',
  'browser position track:1-30', '# metadata comment',
  'track\t0\t1\t-2.5', 'browser\t10\t15\t0', 'track name=literal\t20\t30\t1e-3'), path)
literal <- read_bedgraph(path)
stopifnot(identical(literal$accver, c('track', 'browser', 'track name=literal')),
  identical(literal$value, c(-2.5, 0, .001)),
  nrow(attr(literal, 'ggchord_import')$headers[[1L]]) == 3)
writeLines(c('track\t0\t1', 'browser\t2\t3', 'track name=literal\t4\t5'), path)
bed <- read_bed(path)
stopifnot(identical(bed$accver, c('track', 'browser', 'track name=literal')))
writeLines('browser\t0\t1\ttrack\t2\t4\tpair\t.\t+\t-', path)
paired <- read_bedpe(path)
stopifnot(paired$qaccver == 'browser', paired$saccver == 'track')
writeLines(c('track type=bedGraph name="one"', 'pUC19\t0\t10\t1',
  'track type=bedGraph name="two"', 'pUC19\t10\t20\t2'), path)
fail(read_bedgraph(path), 'one signal')
writeLines('pUC19\t0\t1\tNA', path)
fail(read_bedgraph(path), 'finite real')
writeLines('pUC19\t1\t1\t1', path)
fail(read_bedgraph(path), 'chromEnd > chromStart')
unlink(path)
# No line interpolation across a missing interval: use interval rectangles.
band <- position_track('GC', side = 'inside', width = unit(7, 'mm'))
plot <- ggchord(plasmid_example_pUC19) + geom_seq(show.legend = FALSE) +
  geom_track(aes(xmin = start-1, xmax = end, ymin = 0, ymax = value, fill = value),
    data = text, geom = 'rect', track = band, show.legend = TRUE) +
  scale_track_value_continuous('GC', limits = c(0, 100), unit = 'GC (%)') +
  scale_track_fill_gradientn('GC', colours = c('#F1E51D', '#31849B'), limits = c(0, 100), name = 'GC (%)') +
  coord_circular() + theme_ggchord_publication()
layout <- export_ggchord_layout(plot, original_data = TRUE)
stopifnot(length(layout$track_source_data) == 1L)
retained <- layout$track_source_data[[1L]]
stopifnot(all(vapply(names(text), function(field) identical(retained[[field]], text[[field]]), logical(1))))
ggplot2::ggsave('tools/long-term-audit/v0.27.0-signal-track.png', plot,
  width = 7, height = 7, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', real_source = 'examples/plasmid/pUC19.fna',
  source_windows = nrow(text), full_length = nchar(sequence), interpretation = 'actual GC percent, not read coverage',
  binary_signal_preserved = TRUE, binary_native_precision = TRUE, selection_union_and_clipping = TRUE,
  absent_is_not_zero = TRUE, header_like_identifiers_preserved = TRUE,
  physical_lines_and_headers_retained = TRUE, no_fake_support_or_identity = TRUE),
  'tools/long-term-audit/v0.27.0-signal-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Signal adapters: passed\n')
