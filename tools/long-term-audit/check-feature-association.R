# Original catalogue settings; one PNG device, not release checks.
source('tools/audit-source-snapshot.R')
source('tools/long-term-audit/feature-association-metrics.R')
arguments <- commandArgs(trailingOnly = TRUE)
stage <- if (length(arguments)) arguments[1L] else 'before'
stopifnot(stage %in% c('before', 'after', 'prototype', 'control', 'reverse', 'lanes', 'balanced', 'cluster', 'ordered', 'coordinated', 'mirror'))
audit_session <- audit_source_snapshot(additional = c(
  'tools/long-term-audit/check-feature-association.R',
  'tools/long-term-audit/feature-association-metrics.R'))
pkgload::load_all(quiet = TRUE)
library(dplyr)
cases <- list()
ids <- if (length(arguments)>1L) arguments[-1L] else c('pTRE_Tight_BI', 'pCAMBIA1300')
stopifnot(length(ids)>0L,all(ids %in% c('pUC19','pTRE_Tight_BI','pCAMBIA1300',
  'pSB1C3','pBR322','pBluescript_II_SK_plus','pET_28a_plus','pETDuet_1','pDONR221',
  'pcDNA3_1_plus','pTRIPZ','pEarleyGate_201','pSpCas9_BB_2A_GFP_PX458')))
for (id in ids) {
  cat('Association:', stage, id, '\n'); flush.console()
  env <- new.env()
  data(list = paste0('plasmid_example_', id), package = 'ggchord', envir = env)
  sequence <- env[[paste0('plasmid_example_', id)]]
  features <- find_common_features(sequence)
  sites <- find_restriction_sites(sequence)
  nearest <- unique(vapply(sequence$length*c(1,3,5)/6,
    function(x) which.min(abs(sites$position-x)), integer(1)))
  position <- position_feature_stack(spacing = .085, base_position = position_plasmid())
  plot <- ggchord(sequence[c('accver', 'length')], feature_data = features) +
    geom_seq(seq_style = 'double') + geom_feature_plasmid(position = position) +
    geom_feature_label_repel(position = position, external = TRUE, label_fit = 'shrink',
      label_min_size = 2.5) + geom_restriction_site(data = sites[nearest, ], label_size = 2.5) +
    coord_circular(rotation = if (stage %in% c('reverse','mirror')) 48 else 0,
      direction = if (stage %in% c('reverse','mirror')) 'counterclockwise' else 'clockwise') + theme_ggchord_plasmid()
  path <- paste0('tools/long-term-audit/v0.27.0-association-', stage, '-', id, '.png')
  ggchord_diagnostic_device(path, 'png', 10, 8, 110)
  conditions <- character()
  table <- withCallingHandlers(ggplot2::ggplotGrob(plot), warning = function(w) {
    conditions <<- c(conditions, conditionMessage(w))
    invokeRestart('muffleWarning')
  })
  grid::grid.newpage(); grid::grid.draw(table)
  layout <- attr(table, 'ggchord_render_layout')
  context <- attr(table, 'ggchord_render_context')
  saveRDS(list(layout = layout, context = context, built = attr(table, 'ggchord_render_built')), paste0('/tmp/ggchord-association-', stage, '-', id, '.rds'))
  cases[[id]] <- c(association_metrics(layout), list(file = path,
    input_features = features, selected_sites = sites[nearest, ],
    label_visibility = ggchord_label_visibility(layout), converged = context$converged,
    warnings = conditions, text_quality = association_text_quality(table)$summary))
  if (length(conditions)) cat('Observed warnings:', paste(unique(conditions), collapse=' | '), '\n')
  grDevices::dev.off()
}
audit_assert_source_snapshot(audit_session)
jsonlite::write_json(list(scope = 'Same source tables and original 10x8 inch / 110dpi PNG settings',
  stage = stage, cases = cases, runtime_sha256 = as.list(audit_session$engine),
  input_sha256 = as.list(audit_session$inputs),
  not_checked = c('exact ink/stroke thickness/occlusion', 'arbitrary inputs or release matrix')),
  paste0('tools/long-term-audit/v0.27.0-association-', stage, '-results.json'),
  auto_unbox = TRUE, pretty = TRUE, na = 'null')
cat('Feature association observations completed\n')
