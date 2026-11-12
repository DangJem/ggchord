# Existing observed plant alignments only; the decoder reads an unchanged prefix.
# No aligner, assembly analysis, variant caller or package/release checks.
pkgload::load_all(quiet = TRUE)
folder <- 'examples/genome-alignments'
metadata <- jsonlite::read_json(file.path(folder, 'arabidopsis-source.json'))
path <- file.path(folder, metadata$file)
stopifnot(digest::digest(file = path, algo = 'sha256') == metadata$product_sha256)
program <- Sys.getenv('DELTA_DECODER', 'show-coords')
value <- read_delta(path, executable = program, source_file = TRUE)
links <- value$ribbon_data
raw <- readLines(path)
native <- value$report$records[[1L]]$native_output
stopifnot(nrow(links) == metadata$selected_alignment_records,
  identical(raw[3L], metadata$selected_header), length(native) == nrow(links),
  identical(value$report$raw_files[[1L]], raw),
  all(links$alignment_unit == 'bp'), all(value$seq_data$sequence_unit == 'bp'),
  all(links$.source_file == path), any(links$qstart > links$qend),
  any(links$qstart < links$qend), all(links$sstart <= links$send),
  !any(c('length', 'read_support', 'sequence', 'source_hsp_id') %in% names(links)))
for (i in seq_len(nrow(links))) {
  header <- strsplit(raw[links$source_header_line[i]], '[[:space:]]+')[[1L]]
  fields <- as.numeric(strsplit(raw[links$.source_line[i]], '[[:space:]]+')[[1L]])
  stopifnot(identical(unname(unlist(links[i, c('sstart', 'send', 'qstart', 'qend')])), fields[1:4]),
    links$qaccver[i] == header[2L], links$saccver[i] == substring(header[1L], 2L),
    links$qlen[i] == as.numeric(header[4L]), links$slen[i] == as.numeric(header[3L]),
    links$delta_errors[i] == fields[5L], links$delta_similarity_errors[i] == fields[6L],
    links$delta_stops[i] == fields[7L],
    links$pident[i] == as.numeric(strsplit(native[i], '\t', fixed = TRUE)[[1L]][7L]))
  next_line <- if (i < nrow(links)) links$.source_line[i+1L]-1L else length(raw)
  script <- as.numeric(raw[seq.int(links$.source_line[i]+1L, next_line)])
  stopifnot(identical(links$delta_script[[i]], script), tail(script, 1L) == 0)
}
original <- value$ribbon_data
plot <- ggchord(value$seq_data) + geom_seq(show.legend = FALSE) + geom_seq_label() +
  geom_link_ribbon(data = links, aes(ribbon_fill = pident), alpha = .35) +
  scale_ribbon_fill_gradientn(colours = c('#CAE0E8', '#548CA0', '#254B5B')) +
  coord_collinear() + labs(title = '39 existing Arabidopsis alignment records',
    subtitle = 'Complete source pair; original lengths and directional bp endpoints')
ggplot2::ggsave('tools/long-term-audit/v0.27.0-plant-delta.png', plot,
  width = 10, height = 5, dpi = 110, bg = 'white')
stopifnot(identical(value$ribbon_data, original))
jsonlite::write_json(list(status = 'passed for the unchanged first complete source-pair block',
  source = metadata, alignments = nrow(links), sequence_ids = value$seq_data$accver,
  complete_source_lengths = value$seq_data$length, native_identity_preserved = TRUE,
  every_physical_header_endpoint_error_and_script_checked = TRUE,
  both_query_directions_preserved = TRUE, source_preamble_retained_not_opened = TRUE,
  decoder_sha256 = digest::digest(file = Sys.which(program), algo = 'sha256'),
  no_aligner_or_variant_caller_run = TRUE,
  unverified = c('remaining source sequence-pair blocks', 'character-level consistency without source FASTA',
    'independent reconstruction of exact unrounded identity', 'release device/platform matrix')),
  'tools/long-term-audit/v0.27.0-plant-delta-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('Plant delta source/decoder preservation: passed\n')
