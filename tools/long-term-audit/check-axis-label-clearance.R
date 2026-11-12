# Directed real-source controls, on one PNG device. No package/release checks.
source('tools/audit-source-snapshot.R')
audit_session <- audit_source_snapshot(additional = c('tools/long-term-audit/check-axis-label-clearance.R'))
pkgload::load_all(quiet = TRUE)
audit_assert_source_snapshot(audit_session)
data(plasmid_example_pTRE_Tight_BI)
source <- plasmid_example_pTRE_Tight_BI
features <- find_common_features(source)
position <- position_feature_stack(spacing = .085, base_position = position_plasmid())
base <- ggchord(source[c('accver', 'length')], feature_data = features) +
  geom_seq(seq_style = 'double') + geom_feature_plasmid(position = position) +
  geom_feature_label_repel(position = position, external = TRUE,
    label_fit = 'shrink', label_min_size = 2.5) + coord_circular()
measure <- function(p) {
  file <- tempfile(fileext = '.png')
  ggchord_diagnostic_device(file, 'png', 10, 8, 110)
  on.exit({grDevices::dev.off(); unlink(file)})
  original <- p$ggchord$data
  table <- ggplot2::ggplotGrob(p)
  grid::grid.newpage(); grid::grid.draw(table)
  layout <- attr(table, 'ggchord_render_layout')
  axis <- ggchord_axis_geometry(layout)
  stopifnot(identical(original, p$ggchord$data))
  list(axis = unique(axis[axis$.component %in% 'text' & axis$.axis_seam_combined %in% TRUE,
    intersect(c('.axis_label', '.axis_offset_explicit', '.axis_clearance_status', '.axis_clearance_mm'), names(axis)), drop = FALSE]),
    axis_all = unique(axis[axis$.component %in% 'text',
      intersect(c('.axis_label', '.axis_automatic_breaks', '.axis_clearance_status', '.axis_clearance_mm'), names(axis)), drop = FALSE]),
    ticks = layout$axis_ticks[c('accver', 'x0', 'y0', 'x1', 'y1')])
}
automatic <- measure(base + theme_ggchord_plasmid())
shortcut <- measure(base + theme_ggchord_plasmid(axis.text.offset = grid::unit(1.15, 'mm')))
native <- measure(base + theme_ggchord_plasmid() + theme(ggchord.axis.text.offset = grid::unit(1.15, 'mm')))
combined <- measure(base + (theme_ggchord_plasmid() +
  theme(ggchord.axis.text.offset = grid::unit(1.15, 'mm'))))
unrelated <- measure(base + (theme_ggchord_plasmid() + theme(legend.position = 'none')))
disabled <- measure(base + theme_ggchord_plasmid() + theme(ggchord.axis.text.check.overlap = FALSE))
manual_breaks <- measure(base + theme_ggchord_plasmid() +
  scale_seq_position_continuous(breaks = c(1, 1000, 2000), minor_breaks = NULL))
stopifnot(automatic$axis$.axis_clearance_status == 'feature_clearance',
  abs(automatic$axis$.axis_clearance_mm) > 0,
  shortcut$axis$.axis_clearance_status == 'explicit_offset',
  native$axis$.axis_clearance_status == 'explicit_offset',
  combined$axis$.axis_clearance_status == 'explicit_offset',
  all(combined$axis$.axis_clearance_mm == 0),
  unrelated$axis$.axis_clearance_status == 'feature_clearance',
  all(shortcut$axis$.axis_clearance_mm == 0), all(native$axis$.axis_clearance_mm == 0),
  disabled$axis$.axis_clearance_status == 'disabled',
  all(!manual_breaks$axis_all$.axis_automatic_breaks),
  setequal(manual_breaks$axis_all$.axis_label, c('1', '1000', '2000')),
  isTRUE(all.equal(automatic$ticks, shortcut$ticks)), isTRUE(all.equal(automatic$ticks, native$ticks)))
rotated <- base + theme_ggchord_plasmid() + coord_circular(rotation = 35, direction = 'counterclockwise')
diagnostics <- diagnose_ggchord(rotated, width = 10, height = 8, device = 'png')
stopifnot(diagnostics$label_quality$summary$n_text_overlaps == 0L,
  diagnostics$label_quality$summary$n_text_feature_unresolved == 0L,
  diagnostics$label_quality$summary$n_label_outside_canvas == 0L)
audit_assert_source_snapshot(audit_session)
jsonlite::write_json(list(runtime_sha256 = as.list(audit_session$engine), status = 'passed for recorded source and selected device',
  automatic = automatic$axis, explicit_shortcut = shortcut$axis, explicit_native = native$axis,
  explicit_combined_theme = combined$axis, unrelated_combined_theme = unrelated$axis,
  explicit_breaks = manual_breaks$axis_all,
  overlap_policy_disabled = disabled$axis, tick_anchors_preserved = TRUE,
  equal_valued_explicit_offsets_are_not_inferred_automatic = TRUE,
  rotated_counterclockwise = diagnostics$label_quality$summary,
  source_inputs_preserved = TRUE, scope = 'bounded automatic circular parallel text; no general axis solver or platform certification'),
  if ('--joint' %in% commandArgs(trailingOnly = TRUE))
    'tools/long-term-audit/v0.27.0-axis-clearance-precise-results.json' else
  if ('--refined' %in% commandArgs(trailingOnly = TRUE))
    'tools/long-term-audit/v0.27.0-axis-clearance-refined-results.json' else
    'tools/long-term-audit/v0.27.0-axis-clearance-results.json', pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Automatic seam clearance and explicit offset ownership: passed\n')
