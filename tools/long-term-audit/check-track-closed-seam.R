# Directed regression for a confirmed closed-backbone geometry defect.
pkgload::load_all(quiet = TRUE)
sequence <- plasmid_example_pUC19$sequence[[1L]]
starts <- seq(1, nchar(sequence), by = 100)
windows <- data.frame(accver = 'pUC19', start = starts,
  end = pmin(starts+99, nchar(sequence)))
windows$value <- vapply(seq_len(nrow(windows)), function(i) {
  bases <- strsplit(substr(sequence, windows$start[i], windows$end[i]), '', fixed = TRUE)[[1L]]
  100*mean(bases %in% c('G', 'C'))
}, numeric(1))
base <- ggchord(plasmid_example_pUC19) + geom_seq(show.legend = FALSE) + coord_circular()
layout <- get_chord_layout(base)
positions <- seq(0, nchar(sequence), length.out = 65)
frame <- ggchord_track_frame(positions, 'pUC19', layout)
stopifnot(!ggchord_front_invalid(ggchord_track_point(positions, 'pUC19', -.05, 0, layout), frame$base),
  ggchord_front_invalid(ggchord_track_point(positions, 'pUC19', -2, 0, layout), frame$base))
plot <- base + geom_track(aes(xmin = start-1, xmax = end, ymin = 0, ymax = value, fill = value),
  data = windows, geom = 'rect', track = position_track('GC', side = 'inside', width = unit(7, 'mm'))) +
  scale_track_value_continuous('GC', limits = c(0, 100), unit = 'GC (%)') +
  scale_track_fill_gradientn('GC', colours = c('#F1E51D', '#31849B'), limits = c(0, 100), name = 'GC (%)')
retained <- export_ggchord_layout(plot, original_data = TRUE)$track_source_data[[1L]]
stopifnot(all(vapply(names(windows), function(field) identical(retained[[field]], windows[[field]]), logical(1))))
ggplot2::ggsave('tools/long-term-audit/v0.27.0-track-closed-seam.png', plot,
  width = 7, height = 7, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', narrow_inward_band = TRUE,
  actual_fold_still_detected = TRUE, original_gc_windows_preserved = TRUE,
  source = 'plasmid_example_pUC19', no_biological_circularity_inference = TRUE),
  'tools/long-term-audit/v0.27.0-track-closed-seam-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('Closed Track seam: passed\n')
