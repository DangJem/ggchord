# Direct checks of the audited public contracts; no package suite or exports.
pkgload::load_all(quiet = TRUE)
source_models <- jsonlite::read_json('examples/cds/expected.json', simplifyVector = TRUE)
model <- source_models$models[source_models$models$gene == 'rps16', ]
pieces <- model$segments[[1]]
cds_input <- data.frame(feature_id = model$feature_id, accver = model$accver,
  start = min(pieces$start), end = max(pieces$end), strand = model$strand,
  codon_start = model$codon_start, genetic_code = model$genetic_code)
cds_input$segments <- I(list(pieces))
cds <- build_ggchord_cds(cds_input)
# Selected residues cross the observed exon join. This is a return-shape control,
# not a detected domain or a new biological alignment.
selected <- data.frame(feature_id = 'selected_residues', parent_id = model$feature_id,
  start = 10, end = 20)
feature <- project_ggchord_features(selected, cds)
stopifnot(nrow(feature) == 1L, nrow(feature$segments[[1]]) == 2L,
  feature$feature_id == selected$feature_id, feature$parent_id == selected$parent_id)
# Explicit same-residue correspondence for the utility shape check, not a BLAST HSP.
relation <- data.frame(qfeature_id = model$feature_id, sfeature_id = model$feature_id,
  qstart = 10, qend = 20, sstart = 10, send = 20)
links <- project_ggchord_links(relation, cds, ungapped = TRUE)
stopifnot(nrow(links) > 1L)
window <- data.frame(accver = seq_data_example$accver[1], start = 1, end = 100, locus_id = 'selected_window')
segmented <- segment_ggchord_data(seq_data_example[1, ], windows = window)
stopifnot(nrow(segmented$seq_data) == 1L)
bar_error <- tryCatch(geom_track_bar(baseline = 1), error = identity)
stopifnot(inherits(bar_error, 'error'), grepl('geom_track_col', conditionMessage(bar_error)))
alignment <- read_aligned_fasta('examples/alignment/opuntia-aligned.fasta', 'DNA')
lookup <- map_alignment_position(alignment, alignment$accver[1], c(1, NA_real_))
stopifnot(nrow(lookup) == 2L)
plot <- ggchord(seq_data_example[1, ]) + geom_seq() + coord_collinear()
location <- map_ggchord_position(plot, seq_data_example$accver[1], 100)
stopifnot(nrow(location) == 1L, is.numeric(location$x))
result <- list(feature_output_rows = nrow(feature), feature_segment_rows = nrow(feature$segments[[1]]),
  paired_link_rows = nrow(links), windows_argument = 'accepted', nonzero_bar_baseline = 'rejected with col alternative',
  alignment_numeric_position = 'accepted', genomic_numeric_position = 'accepted',
  scope = 'Selected residues on source-backed rps16 exons and minimal documented API calls; no inferred biological alignment, package suite or device matrix')
jsonlite::write_json(result, 'tools/site-audit/semantic-check.json', auto_unbox = TRUE, pretty = TRUE)
cat('PASS: feature/list-column versus link row splitting, windows, bar baseline, numeric position contracts\n')
