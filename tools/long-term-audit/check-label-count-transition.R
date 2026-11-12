# Real unmodified site records; only the explicitly selected row count changes.
helper <- readLines('tools/long-term-audit/check-label-quality.R', warn = FALSE)
end <- grep('^# Explicit geometry controls', helper)[1L]-1L
stopifnot(is.finite(end))
eval(parse(text = helper[seq_len(end)]))
audit_assert_source_snapshot(audit_session)
audit_session <- audit_source_snapshot(additional = c(audit_session$additional, 'tools/long-term-audit/check-label-count-transition.R'))
prefix <- 'v0.27.0-label-transition-'
data(plasmid_example_pUC19)
sequence <- plasmid_example_pUC19
features <- find_common_features(sequence)
sites <- find_restriction_sites(sequence)
counts <- c(0L, 1L, 3L, 10L, 29L, 30L, 31L, 58L)
stopifnot(nrow(sites) >= max(counts))
selected_counts <- sub('^--counts=', '', grep('^--counts=', commandArgs(trailingOnly = TRUE), value = TRUE))
if (length(selected_counts)) {
  stopifnot(length(selected_counts) == 1L)
  selected_counts <- as.integer(strsplit(selected_counts, ',', fixed = TRUE)[[1L]])
  stopifnot(length(selected_counts) > 0L, all(selected_counts %in% counts))
  counts <- selected_counts
}
canvas <- c(width = 10, height = 8)
for (axis in c('width', 'height')) {
  option <- grep(paste0('^--', axis, '='), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(option)) {
    stopifnot(length(option) == 1L)
    value <- as.numeric(sub(paste0('^--', axis, '='), '', option))
    stopifnot(is.finite(value), value > 0)
    canvas[axis] <- value
  }
}
canvas_suffix <- if (identical(unname(canvas), c(10, 8))) '' else
  paste0('-', canvas[['width']], 'x', canvas[['height']])
prefix <- paste0('v0.27.0-label-transition-device', canvas_suffix, '-')
position <- position_feature_stack(spacing = .085, base_position = position_plasmid())
results <- list()
for (count in counts) {
  cat('Site count transition:', count, '\n'); flush.console()
  selected <- sites[seq_len(count), , drop = FALSE]
  plot <- ggchord(sequence[c('accver', 'length')], feature_data = features) +
    geom_seq(seq_style = 'double') + geom_feature_plasmid(position = position) +
    geom_feature_label_repel(position = position, external = TRUE,
      label_fit = 'shrink', label_min_size = 2.5) +
    geom_restriction_site(data = selected, label_size = 2.5) +
    coord_circular() + theme_ggchord_plasmid()
  results[[as.character(count)]] <- c(measure(plot, paste0('sites-', count), canvas[['width']], canvas[['height']]),
    list(selected_original_rows = seq_len(count), selection = 'unchanged prefix of find_restriction_sites output'))
  stopifnot(identical(selected, sites[seq_len(count), , drop = FALSE]),
    results[[as.character(count)]]$n_source_labels == nrow(features)+count,
    results[[as.character(count)]]$n_visible_source_labels == nrow(features)+count)
}
audit_assert_source_snapshot(audit_session)
sources <- list.files('R', '[.]R$', full.names = TRUE)
jsonlite::write_json(list(scope = paste('pUC19 feature labels with', length(counts), 'selected real-site row count(s) on the declared PNG canvas; not arbitrary counts or devices'),
  canvas_inches = as.list(canvas), source_features = nrow(features), source_sites = nrow(sites), cases = results,
  runtime_sha256 = as.list(audit_session$engine),
  source_tables_preserved = TRUE, per_case_manual_tuning = 0L,
  caveats = 'Contacts are observations; changing source count can legitimately change placement. No pixel-equality or unmeasured quality guarantee.'),
  paste0('tools/long-term-audit/v0.27.0-label-transition-device', canvas_suffix, '-',
    if (length(selected_counts)) paste(counts, collapse = '-') else 'all', '-results.json'), pretty = TRUE, auto_unbox = TRUE)
cat('Site count transition observations recorded\n')
