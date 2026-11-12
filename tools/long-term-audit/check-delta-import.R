# Active real example is Drosophila PROMER. Historical nucleotide example stays
# archived with its original provenance; never load it for active checks.
# Only show-coords parses existing data.
# No aligner, package test suite or release device matrix.
pkgload::load_all(quiet = TRUE)
program <- Sys.getenv('DELTA_DECODER', 'show-coords')
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(pattern, conditionMessage(error)))
}
results <- list(); imported <- list()
for (mode in 'promer') {
  path <- paste0('examples/genome-alignments/', mode, '.delta')
  value <- read_delta(path, executable = program, source_file = TRUE)
  data <- value$ribbon_data
  raw <- readLines(path)
  native <- value$report$records[[1L]]$native_output
  stopifnot(nrow(data) > 0, all(data$.source_file == path),
    identical(value$report$raw_files[[1L]], raw), length(native) == nrow(data),
    all(value$seq_data$sequence_unit == 'bp'),
    all(vapply(data$delta_script, function(x) length(x) > 0 && tail(x, 1) == 0, logical(1))),
    !any(c('length', 'read_support', 'sequence', 'source_hsp_id') %in% names(data)))
  # Compare every returned row to its actual physical source alignment header.
  for (i in seq_len(nrow(data))) {
    header <- strsplit(trimws(raw[data$source_header_line[i]]), '[[:space:]]+')[[1L]]
    fields <- as.numeric(strsplit(trimws(raw[data$.source_line[i]]), '[[:space:]]+')[[1L]])
    stopifnot(data$qaccver[i] == header[2L], data$saccver[i] == substring(header[1L], 2L),
      identical(unname(unlist(data[i, c('sstart', 'send', 'qstart', 'qend')])), fields[1:4]),
      data$delta_errors[i] == fields[5L], data$delta_similarity_errors[i] == fields[6L],
      data$delta_stops[i] == fields[7L],
      data$qlen[i] == as.numeric(header[4L]), data$slen[i] == as.numeric(header[3L]))
    decoded <- strsplit(native[i], '\t', fixed = TRUE)[[1L]]
    stopifnot(data$pident[i] == as.numeric(decoded[7L]))
  }
  stopifnot(any(data$qstart > data$qend), any(data$qstart < data$qend),
    all(data$qspan == abs(data$qend-data$qstart)+1),
    all(data$sspan == abs(data$send-data$sstart)+1))
  if (mode == 'promer') stopifnot(all(data$alignment_unit == 'aa'),
    all(c('similarity_percent', 'stop_percent', 'subject_frame', 'query_frame') %in% names(data)),
    any(abs(data$query_frame) > 1)) else stopifnot(all(data$alignment_unit == 'bp'))
  imported[[mode]] <- value
  results[[mode]] <- list(source_sha256 = digest::digest(file = path, algo = 'sha256'),
    original_alignments = nrow(data), native_identity_preserved = TRUE,
    identity_unit = unique(data$alignment_unit), coordinate_unit = 'bp',
    full_source_lengths_and_direction = TRUE, exact_physical_source_rows = TRUE,
    signed_scripts_retained = TRUE)
}
# Synthetic nucleotide format controls, not source-backed observations.
fixture <- tempfile(fileext = '.delta')
writeLines(c('format-reference.fna format-query.fna', 'NUCMER', '>format_reference format_query 100 200',
  '1 10 20 29 0 0 0', '0', '21 30 100 91 1 1 0', '0'), fixture)
format_control <- read_delta(fixture, executable = program)
stopifnot(nrow(format_control$ribbon_data) == 2L,
  identical(format_control$ribbon_data$qstart, c(20, 100)),
  identical(format_control$ribbon_data$qend, c(29, 91)),
  all(format_control$ribbon_data$alignment_unit == 'bp'))
# Synthetic format-only identity collision, not new biological evidence.
writeLines(c('reference.fna query.fna', 'NUCMER', '>same same 100 200',
  '1 10 200 191 0 0 0', '0'), fixture)
fail(read_delta(fixture, executable = program), 'IDs collide')
aliased <- read_delta(fixture, executable = program,
  query_ids = c(same = 'query.same'), subject_ids = c(same = 'reference.same'))
stopifnot(aliased$ribbon_data$qaccver == 'query.same',
  aliased$ribbon_data$saccver == 'reference.same', aliased$ribbon_data$q_source_accver == 'same',
  aliased$ribbon_data$qstart == 200, aliased$ribbon_data$qend == 191)
unlink(fixture)
# Plot a declared subset of existing Drosophila relations; no synthetic links.
value <- imported$promer
selected <- value$ribbon_data[value$ribbon_data$qaccver == '3210101', ]
ids <- unique(c(selected$qaccver, selected$saccver))
sequences <- value$seq_data[value$seq_data$accver %in% ids, ]
plot <- ggchord(sequences) + geom_seq(show.legend = FALSE) + geom_seq_label() +
  geom_link_ribbon(data = selected) + coord_collinear() +
  labs(title = 'Existing translated alignment intervals',
    subtitle = 'Drosophila source coordinates in bp; native identity in AA')
exported <- export_ggchord_layout(plot, original_data = TRUE)
ggplot2::ggsave('tools/long-term-audit/v0.27.0-delta-active-links.png', plot,
  width = 9, height = 4.5, dpi = 120, bg = 'white')
jsonlite::write_json(list(status = 'passed', modes = results,
  explicit_identity_maps = TRUE, no_aligner_run = TRUE,
  active_real_sample = 'unchanged upstream Drosophila PROMER; no historical nucleotide sample loaded',
  nucleotide_control = 'explicitly synthetic two-record NUCMER schema/direction control; not biological observations',
  no_reconstructed_characters_or_support = TRUE,
  decoder = list(executable = program, source_tag_for_this_check = 'v4.0.0rc1',
    executable_sha256 = digest::digest(file = Sys.which(program), algo = 'sha256'))),
  'tools/long-term-audit/v0.27.0-delta-active-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('Delta adapter: passed\n')
