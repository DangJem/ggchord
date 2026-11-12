# Native source-space estimation comparison, not a package/device matrix.
pkgload::load_all(quiet = TRUE)
library(dplyr)
data(plasmid_example_pUC19)
sequence <- plasmid_example_pUC19
bins <- read_bedgraph('inst/extdata/pUC19-gc-100.bedGraph') |>
  mutate(pos = (start + end) / 2, metric = 'GC (%)')
events <- strsplit(sequence$sequence, '', fixed = TRUE)[[1L]]
events <- data.frame(accver = sequence$accver, pos = seq_along(events), metric = events) |>
  filter(metric %in% c('G', 'C'))
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(pattern, conditionMessage(error)))
}
get_statistics <- function(data, mapping, stat, parameters = list(), geom = 'line') {
  layer <- do.call(geom_track, c(list(mapping = mapping, data = data, geom = geom, stat = stat,
    track = position_track(stat), position = 'identity'), parameters))
  plot <- ggchord(sequence) + geom_seq() + layer + coord_collinear()
  original <- plot$ggchord$data
  result <- export_ggchord_layout(plot, include = 'track')
  stopifnot(identical(original, plot$ggchord$data))
  # The same mature Stat, parameters and source-space inputs, with no projected
  # canvas coordinates or Track numerical transformation entering estimation.
  native_layer <- do.call(getExportedValue('ggplot2', paste0('geom_', geom)),
    c(list(mapping = mapping, data = data, stat = stat, position = 'identity'), parameters))
  native <- ggplot2::ggplot_build(ggplot2::ggplot() + native_layer)$data[[1L]]
  computed <- result$track_statistics[[1L]]
  stopifnot(nrow(computed) == nrow(native), nrow(computed) > 0L)
  for (key in intersect(c('x', 'y', 'ymin', 'ymax', 'density', 'ecdf', 'group'), names(native)))
    stopifnot(isTRUE(all.equal(computed[[key]], native[[key]], check.attributes = FALSE)))
  for (i in seq_len(nrow(computed))) {
    rows <- computed$source_rows[[i]]
    stopifnot(length(rows) > 0L, length(unique(data$metric[rows])) == 1L,
      identical(sort(rows), which(data$metric == data$metric[rows[1L]])))
  }
  list(plot = plot, report = result, data = computed)
}
smooth <- get_statistics(bins, aes(x = pos, y = value, group = metric), 'smooth',
  list(method = 'lm', formula = y ~ x, se = TRUE, n = 40L))
density <- get_statistics(events, aes(x = pos, y = after_stat(density), colour = metric, group = metric),
  'density', list(bw = 100, trim = TRUE, n = 64L))
ecdf <- get_statistics(events, aes(x = pos, y = after_stat(ecdf), colour = metric, group = metric),
  'ecdf', list(pad = FALSE), geom = 'step')
stopifnot(all(ecdf$data$y >= 0 & ecdf$data$y <= 1), all(is.finite(ecdf$data$x)))
fail(geom_track(aes(x = pos), data = events, stat = 'ecdf', track = position_track('ecdf')), 'pad=FALSE')
# Same raw estimates under an explicit transformed Track value scale; the
# transform is a display mapping, never a silent log-space model fit.
transformed <- smooth$plot + scale_track_value_continuous('smooth', transform = 'log10')
stopifnot(identical(export_ggchord_layout(transformed, include = 'track')$track_statistics,
  smooth$report$track_statistics))
# Same group observations through curved and vertical projection, not a new fit.
for (coord in list(coord_circular(), coord_collinear(seq_axis = 'vertical'))) {
  out <- export_ggchord_layout(smooth$plot + coord, include = 'track')
  stopifnot(identical(out$track_statistics, smooth$report$track_statistics))
}
# NA-only formatting control: native removes missing y before group membership.
control <- rbind(bins, transform(bins[1L, ], value = NA_real_))
plot <- ggchord(sequence) + geom_seq() + geom_track(aes(x = pos, y = value, group = metric),
  data = control, geom = 'line', stat = 'smooth', method = 'lm', formula = y ~ x,
  se = FALSE, na.rm = TRUE, track = position_track('eligible')) + coord_collinear()
eligible <- export_ggchord_layout(plot, include = 'track')$track_statistics[[1L]]
stopifnot(all(vapply(eligible$source_rows, function(x) identical(x, seq_len(nrow(bins))), logical(1))))
# Independent CI ribbon and mean line share the existing physical band/scale.
plot <- ggchord(sequence) + geom_seq(show.legend = FALSE) +
  geom_track(aes(x = pos, y = value), data = bins, geom = 'ribbon', stat = 'smooth',
    method = 'lm', formula = y ~ x, se = TRUE, fill = '#387D88', alpha = .2,
    track = position_track('GC', width = unit(15, 'mm'))) +
  geom_track(aes(x = pos, y = value), data = bins, geom = 'line', stat = 'smooth',
    method = 'lm', formula = y ~ x, se = TRUE, colour = '#285B65',
    track = position_track('GC', width = unit(15, 'mm'))) +
  scale_track_value_continuous('GC', name = 'GC', unit = '%') + coord_collinear() +
  labs(title = 'GC windows: native linear fit and confidence band',
    subtitle = 'Existing pUC19 GC windows; both layers share source observations') +
  theme(plot.title.position = 'plot', plot.subtitle = element_text(hjust = .5))
ggplot2::ggsave('tools/long-term-audit/v0.27.0-track-estimates.png', plot,
  width = 10, height = 4.5, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', ggplot2 = as.character(packageVersion('ggplot2')),
  sources = 'unchanged bundled pUC19 sequence and GC bedGraph; inst/extdata/signal-sources.json',
  stats = c('smooth', 'density', 'ecdf'), native_source_estimates_equal = TRUE,
  group_source_membership_before_position = TRUE, missing_rows_excluded = TRUE,
  raw_estimates_before_track_transforms = TRUE, curved_vertical_projections_preserve_statistics = TRUE,
  finite_ecdf_explicit_pad = FALSE, model_output_not_new_observations = TRUE),
  'tools/long-term-audit/v0.27.0-track-estimate-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Track source estimates: passed\n')
