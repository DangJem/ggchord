# Focused audit probes, with explicitly synthetic inputs; not a package test suite.
# Records observations rather than accepting the current defects as a contract.
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver = c('A', 'B'), length = 1000)
results <- list(baseline = 'f3d9ef4', R = R.version.string,
  ggplot2 = as.character(utils::packageVersion('ggplot2')))

observations <- data.frame(accver = 'A', pos = c(100, 100), category = c('one', 'two'))
p <- ggchord(seqs) + geom_seq() + geom_track(
  aes(x = pos, fill = category, group = 1), data = observations,
  geom = 'bar', track = position_track('counts')) + coord_collinear()
track <- get_chord_layout(p)$track_statistics[[1]]
reference <- ggplot2::ggplot(observations, aes(x = pos, fill = category, group = 1)) +
  ggplot2::geom_bar()
native <- suppressWarnings(ggplot2::ggplot_build(reference)$data[[1]])
results$explicit_group <- list(track_count = track$count, native_count = native$count,
  defect_observed = !identical(track$count, native$count))

points <- data.frame(accver = 'A', pos = c(100, 150, 200))
p <- ggchord(seqs) + geom_seq() + geom_track(aes(x = pos), data = points,
  geom = 'col', stat = 'bin', breaks = c(100, 200),
  track = position_track('bins')) + coord_collinear()
binned <- get_chord_layout(p)$track_statistics[[1]]
results$bin_source_membership <- list(count = binned$count,
  source_rows = unname(binned$source_rows), expected_source_rows = list(1:3),
  defect_observed = !identical(unname(binned$source_rows[[1]]), 1:3))

hits <- data.frame(qaccver = 'A', saccver = 'B', qstart = 100, qend = 200,
  sstart = 300, send = 400, pident = 90, source_hsp_id = 'observed-1')
p <- ggchord(seqs, hits) + geom_seq() + geom_link_ribbon() + geom_link_ribbon()
acceptance <- lapply(list(chord = coord_chord(), collinear = coord_collinear(),
  radial = coord_radial_synteny()), function(coord) {
  tryCatch({ get_chord_layout(p + coord); list(accepted = TRUE) },
    error = function(e) list(accepted = FALSE, message = conditionMessage(e)))
})
results$declared_duplicate_hsp <- list(views = acceptance,
  defect_observed = acceptance$chord$accepted && acceptance$collinear$accepted && !acceptance$radial$accepted)

source <- ggchord(seqs) + geom_seq(seq_order = c(A1 = 'A', B = 'B', A2 = 'A')) + coord_chord()
detail <- ggchord(seqs) + geom_seq() + coord_collinear()
region <- data.frame(accver = 'A', start = 100, end = 200, instance_note = 'A2')
selected <- export_ggchord_fusion(magnify_ggchord(source, region, detail))$focus_regions
plain <- region[c('accver', 'start', 'end')]
without_note <- tryCatch({export_ggchord_fusion(magnify_ggchord(source, plain, detail)); 'accepted'},
  error = conditionMessage)
results$metadata_selector <- list(actual_source_instance = selected$source_instance,
  without_metadata = without_note,
  defect_observed = identical(selected$source_instance, 'A2') && grepl('ambiguous', without_note))

jsonlite::write_json(results, 'tools/long-term-audit/review-v10-2026-11-12-results.json',
  pretty = TRUE, auto_unbox = TRUE)
cat('Audit probes completed; observed defects:',
  sum(vapply(results[4:7], function(x) isTRUE(x$defect_observed), logical(1))), '\n')
