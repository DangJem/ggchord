# Directed audit of omitted and explicitly specified fixed line styles.
# A finite real plant relation table; no package checks or device matrix.
pkgload::load_all(quiet = TRUE)
plant <- read_delta('examples/genome-alignments/arabidopsis-first-pair.delta',
  executable = Sys.getenv('DELTA_DECODER', 'show-coords'))
base <- ggchord(plant$seq_data, plant$ribbon_data)
# Explicit observed alignment-start pairs, not inferred biological breakpoints.
points <- plant$ribbon_data |>
  dplyr::mutate(qpos = qstart, spos = sstart)
cases <- list(sequence_default = list(layer = geom_seq(), field = 'linewidth', omitted = TRUE),
  sequence_explicit = list(layer = geom_seq(linewidth = .9), field = 'linewidth', omitted = FALSE),
  link_default = list(layer = geom_link_line(data = points), field = 'link_linewidth', omitted = TRUE),
  link_explicit = list(layer = geom_link_line(data = points, linewidth = .6), field = 'link_linewidth', omitted = FALSE),
  ribbon_default = list(layer = geom_link_ribbon(), field = 'linewidth', omitted = TRUE),
  ribbon_explicit = list(layer = geom_link_ribbon(linewidth = .25), field = 'linewidth', omitted = FALSE),
  alignment_default = list(layer = geom_alignment_segment(), field = 'link_linewidth', omitted = TRUE),
  alignment_explicit = list(layer = geom_alignment_segment(linewidth = .4), field = 'link_linewidth', omitted = FALSE))
if ('--after' %in% commandArgs(trailingOnly = TRUE)) {
  cases$sequence_mapped <- list(layer = geom_seq(aes(linewidth = length)),
    field = 'linewidth', omitted = TRUE, mapped = TRUE)
  cases$link_mapped <- list(layer = geom_link_line(aes(link_linewidth = pident), data = points),
    field = 'link_linewidth', omitted = TRUE, mapped = TRUE)
  cases$sequence_inherited <- list(layer = geom_seq(inherit.aes = TRUE),
    base_mapping = aes(linewidth = length), field = 'linewidth', omitted = TRUE, mapped = TRUE)
  cases$sequence_fixed_over_mapping <- list(layer = geom_seq(aes(linewidth = length), linewidth = .9),
    field = 'linewidth', omitted = FALSE)
}
after <- '--after' %in% commandArgs(trailingOnly = TRUE)
observations <- list()
file <- tempfile(fileext = '.png')
grDevices::png(file, width = 12, height = 12, units = 'in', res = 72)
for (name in names(cases)) {
  case <- cases[[name]]
  plot <- base
  if (!is.null(case$base_mapping)) plot <- plot + case$base_mapping
  plot <- plot + case$layer + adapt_ggchord('scale')
  original_data <- plot$ggchord$data
  requested <- plot$layers[[1]]$aes_params[[case$field]]
  adjusted <- ggchord_apply_adaptive_styles(plot)
  effective <- adjusted$layers[[1]]$aes_params[[case$field]]
  stopifnot(identical(plot$ggchord$data, original_data),
    identical(plot$layers[[1]]$aes_params[[case$field]], requested))
  if (after) {
    stopifnot(is.numeric(requested), is.numeric(effective))
    if (case$omitted && !isTRUE(case$mapped)) stopifnot(effective > requested) else stopifnot(identical(effective, requested))
  }
  observations[[name]] <- list(aesthetic = case$field, omitted = case$omitted,
    requested = requested, effective = effective,
    mapped = isTRUE(case$mapped),
    line_factor = adjusted$ggchord$adaptive_report$line_factor,
    source_data_and_original_layer_preserved = TRUE)
}
grDevices::dev.off(); unlink(file)
rendered <- list()
if (after) {
  measure <- function(plot, id, keep = FALSE) {
    file <- if (keep) paste0('tools/long-term-audit/v0.27.0-adaptive-default-lines-', id, '.png') else tempfile(fileext = '.png')
    ggchord_diagnostic_device(file, 'png', 12, 12, 90)
    on.exit({grDevices::dev.off(); if (!keep) unlink(file)})
    source <- plot$ggchord$data
    table <- ggplot2::ggplotGrob(plot)
    grid::grid.newpage(); grid::grid.draw(table)
    stopifnot(identical(source, plot$ggchord$data))
    built <- attr(table, 'ggchord_render_built')
    list(linewidth = built$data[[1]]$linewidth,
      link_linewidth = if (length(built$data) > 1L) built$data[[2]]$link_linewidth else NULL,
      source_preserved = TRUE, file = if (keep) file else NULL)
  }
  mapped <- base + geom_seq(aes(linewidth = length))
  original <- measure(mapped + adapt_ggchord('fixed'), 'mapped-fixed')
  scaled <- measure(mapped + adapt_ggchord('scale'), 'mapped-scale')
  stopifnot(length(original$linewidth) > 0L,
    isTRUE(all.equal(original$linewidth, scaled$linewidth)))
  fixed_mapped <- measure(base + geom_seq(aes(linewidth = length), linewidth = .9) +
    adapt_ggchord('scale'), 'fixed-over-mapping')
  stopifnot(length(fixed_mapped$linewidth) > 0L, all(fixed_mapped$linewidth == .9))
  pair <- base + geom_seq() + geom_seq(linewidth = .9) + adapt_ggchord('scale')
  pair_device <- tempfile(fileext = '.png')
  grDevices::png(pair_device, width = 12, height = 12, units = 'in', res = 72)
  adjusted_pair <- ggchord_apply_adaptive_styles(pair)
  stopifnot(adjusted_pair$layers[[1]]$aes_params$linewidth > .9,
    adjusted_pair$layers[[2]]$aes_params$linewidth == .9,
    pair$layers[[1]]$aes_params$linewidth == .9,
    pair$layers[[2]]$aes_params$linewidth == .9)
  fixed <- ggchord_apply_adaptive_styles(base + geom_seq() + adapt_ggchord('fixed'))
  stopifnot(fixed$layers[[1]]$aes_params$linewidth == .9)
  grDevices::dev.off(); unlink(pair_device)
  combined <- base + geom_seq() + geom_link_line(data = points) + adapt_ggchord('scale')
  rendered <- list(mapped_final_values_preserved = TRUE,
    explicit_fixed_width_wins_over_mapping = TRUE,
    two_same_type_layers_independent = TRUE, fixed_policy_unchanged = TRUE,
    combined = measure(combined, 'combined', TRUE))
}
sources <- list.files('R', '[.]R$', full.names = TRUE)
jsonlite::write_json(list(scope = 'Default versus explicit fixed line aesthetics at a larger physical canvas; eight native constructor controls and selected mapped/inherited counterparts',
  observations = observations, rendered = rendered, runtime_sha256 = as.list(stats::setNames(vapply(sources,
    function(path) digest::digest(file = path, algo = 'sha256'), character(1)), sources)),
  caveats = 'Not every geometry constant, native from_theme expression, exact ink or platform certification'),
  paste0('tools/long-term-audit/v0.27.0-adaptive-default-lines', if (after) '-after' else '-before', '-results.json'),
  pretty = TRUE, auto_unbox = TRUE)
cat('Selected adaptive default line observations recorded\n')
