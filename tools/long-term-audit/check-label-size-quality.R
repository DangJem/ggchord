# Controlled canvas experiment, one PNG device. It records failures as evidence;
# no case-specific tuning, package tests or publication device matrix.
pkgload::load_all(quiet = TRUE)
folder <- 'tools/long-term-audit'
adaptive <- '--adaptive' %in% commandArgs(trailingOnly = TRUE)
prefix <- if (adaptive) 'v0.27.0-label-size-adaptive-' else 'v0.27.0-label-size-'
if ('--after' %in% commandArgs(trailingOnly = TRUE)) prefix <- paste0(prefix, 'after-')
if ('--refined' %in% commandArgs(trailingOnly = TRUE)) prefix <- paste0(prefix, 'refined-')
if ('--current' %in% commandArgs(trailingOnly = TRUE)) prefix <- paste0(prefix, 'current-')
cases <- c('pUC19', 'pTRE_Tight_BI', 'pETDuet_1')
profiles <- if (adaptive) list(compact = c(6, 6)) else list(compact = c(6, 6), wide = c(12, 8))
if ('--compact-only' %in% commandArgs(trailingOnly = TRUE)) profiles <- profiles['compact']
observations <- list()
for (id in cases) {
  inputs <- new.env()
  data(list = paste0('plasmid_example_', id), package = 'ggchord', envir = inputs)
  sequence <- inputs[[paste0('plasmid_example_', id)]]
  features <- find_common_features(sequence)
  sites <- find_restriction_sites(sequence)
  nearest <- if (nrow(sites)) unique(vapply(sequence$length*c(1, 3, 5)/6,
    function(x) which.min(abs(sites$position-x)), integer(1))) else integer()
  position <- position_feature_stack(spacing = .085, base_position = position_plasmid())
  plot <- ggchord(sequence[c('accver', 'length')], feature_data = features) +
    geom_seq(seq_style = 'double') + geom_feature_plasmid(position = position) +
    geom_feature_label_repel(position = position, external = TRUE,
      label_fit = 'shrink', label_min_size = 2.5) +
    geom_restriction_site(data = sites[nearest, ], label_size = 2.5) +
    coord_circular() + theme_ggchord_plasmid()
  if (adaptive) plot <- plot + adapt_ggchord('adaptive')
  original <- plot$ggchord$data
  for (profile in names(profiles)) {
    cat('Controlled canvas:', id, profile, '\n'); flush.console()
    size <- profiles[[profile]]
    path <- file.path(folder, paste0(prefix, id, '-', profile, '.png'))
    ggchord_diagnostic_device(path, 'png', size[1], size[2], 110)
    table <- ggplot2::ggplotGrob(plot)
    grid::grid.newpage(); grid::grid.draw(table)
    built <- attr(table, 'ggchord_render_built')
    layout <- attr(table, 'ggchord_render_layout')
    context <- attr(table, 'ggchord_render_context')
    text <- ggchord_rbind_fill(lapply(seq_along(built$data), function(i) {
      d <- built$data[[i]]; layer <- built$plot$layers[[i]]
      d <- ggchord_diagnostic_text_styles(d, layer)
      txt <- ggchord_diagnostic_text(d, built$plot$coordinates,
        built$layout$panel_params[[1L]], context,
        parse = layer$geom_params$parse %||% FALSE,
        size_unit = layer$geom_params$size.unit %||% 'mm')
      if (!nrow(txt)) return(NULL)
      txt$layer_id <- layer$ggchord_layer_id %||% paste0('native-', i)
      txt
    })) %||% data.frame()
    quality <- ggchord_label_quality(layout, text, built, context)
    visibility <- ggchord_label_visibility(layout)
    axis <- ggchord_axis_geometry(layout)
    clearance <- if ('.axis_clearance_status' %in% names(axis)) unique(axis[
      axis$.component %in% 'text' & !is.na(axis$.axis_clearance_status),
      intersect(c('accver', '.axis_label', '.axis_automatic_breaks', '.axis_offset_explicit',
        '.axis_clearance_status', '.axis_clearance_mm'), names(axis)), drop = FALSE]) else data.frame()
    grDevices::dev.off()
    stopifnot(identical(original, plot$ggchord$data), context$status == 'measured')
    observations[[paste(id, profile, sep = '/')]] <- list(
      width_in = size[1], height_in = size[2], quality = quality,
      visibility = visibility, font_range_mm = range(text$size_mm, na.rm = TRUE),
      axis_clearance = clearance,
      source_preserved = TRUE, case_specific_manual_changes = 0L,
      file = path, png_sha256 = digest::digest(file = path, algo = 'sha256'))
  }
}
sources <- list.files('R', '[.]R$', full.names = TRUE)
jsonlite::write_json(list(scope = 'Three unchanged catalogue sources, common layer settings and explicit physical sizes; one PNG device. Observations are not assumed passes.',
  layout_policy = if (adaptive) 'adaptive' else 'existing fixed physical styles',
  observations = observations, runtime_sha256 = as.list(stats::setNames(vapply(sources,
    function(path) digest::digest(file = path, algo = 'sha256'), character(1)), sources)),
  unverified = c('exact ink, non-polygon strokes/markers and font portability',
    'all catalogue sources at arbitrary sizes', 'publication device/platform matrix')),
  file.path(folder, paste0(prefix, 'results.json')),
  pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Controlled canvas observations recorded\n')
