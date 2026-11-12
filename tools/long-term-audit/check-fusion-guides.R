# Direct outer-guide budget checks; no package tests or release device matrix.
pkgload::load_all(quiet = TRUE)
library(dplyr)
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
data(seq_data_example); data(ribbon_data_example)
base <- ggchord(seq_data_example) + geom_seq(show.legend = FALSE)
hit <- ribbon_data_example |> slice_head(n = 1)
identity_fusion <- function(layout = 'fixed') {
  ggchord_fusion(base, layout = layout) +
    fusion_panel(base + coord_collinear(), 'detail') +
    fusion_ribbon(hit, qpanel = 'center', spanel = 'detail')
}
fixed <- export_ggchord_fusion(identity_fusion(), 12, 8)
compact <- export_ggchord_fusion(identity_fusion('adaptive'), 12, 3)
g <- compact$metadata$guide_layout$guides
stopifnot(fixed$metadata$guide_side == 'right',
  compact$metadata$guide_layout$status == 'feasible',
  g$direction == 'horizontal', g$final_position %in% c('top', 'bottom'),
  g$n_items == 6, !g$position_locked, !g$direction_locked,
  identical(compact$relationship_data[[1]], hit))
# Multiple native guides own their individual positions; scales stay native.
make <- function(layout = 'fixed') {
  ggchord_fusion(base, layout = layout) +
    fusion_panel(base + coord_collinear(), 'detail') +
    fusion_ribbon(ribbon_data_example, aes(fill = qaccver, alpha = pident),
      qpanel = 'center', spanel = 'detail') +
    ggplot2::guides(fill = ggplot2::guide_legend(position = 'bottom', nrow = 1),
      alpha = ggplot2::guide_legend(position = 'right'))
}
a <- export_ggchord_fusion(make(), 12, 7)
b <- export_ggchord_fusion(make('adaptive'), 12, 7)
stopifnot(a$metadata$guide_side == 'multiple', b$metadata$guide_side == 'multiple',
  setequal(b$metadata$guide_bounds_mm$side, c('bottom', 'right')),
  all(b$metadata$guide_layout$guides$position_locked),
  b$metadata$guide_layout$guides$rows[b$metadata$guide_layout$guides$aesthetic == 'fill'] == 1,
  identical(b$relationship_data[[1]], ribbon_data_example),
  nrow(b$connections) == nrow(ribbon_data_example))
# Source-derived colour and opacity mappings remain the same after placement.
styles <- function(d) d$paths |> distinct(layer_id, source_row, fill, alpha) |> arrange(layer_id, source_row)
stopifnot(identical(styles(a), styles(b)))
boxes <- b$metadata$guide_bounds_mm
for (i in seq_len(nrow(boxes))) {
  guide <- boxes[i, ]
  panels <- b$panels
  stopifnot(!any(pmin(guide$xmax, panels$container_xmax) > pmax(guide$xmin, panels$container_xmin) &
    pmin(guide$ymax, panels$container_ymax) > pmax(guide$ymin, panels$container_ymin)))
}
# Explicit user layout wins, even if it makes a narrow canvas impossible.
locked <- identity_fusion('adaptive') + ggplot2::theme(legend.position = 'right',
  legend.direction = 'vertical', legend.key.height = unit(90, 'mm'),
  legend.text = element_text(size = 12))
fail(export_ggchord_fusion(locked, 12, 3), 'cannot fit')
# An absolute panel box participates in the candidate's remaining-space budget.
physical <- ggchord_fusion(base, layout = 'adaptive', center_box = unit(c(5, 245, 5, 65), 'mm')) +
  fusion_ribbon(hit |> mutate(qpanel = 'center', spanel = 'detail')) +
  fusion_panel(base + coord_collinear(), 'detail', box = unit(c(250, 295, 5, 65), 'mm'))
# At this height, the long default ruler must be placed horizontally; absolute
# boxes leave no legal horizontal reserve, so silently moving/shrinking them is forbidden.
fail(export_ggchord_fusion(physical, 12, 3), 'cannot fit')
# Inside composition guides remain rejected, rather than silently relocated.
inside <- identity_fusion('adaptive') +
  ggplot2::guides(fill = ggplot2::guide_coloursteps(position = 'inside'))
fail(export_ggchord_fusion(inside, 12, 8), 'cannot fit')
none <- identity_fusion('adaptive') + ggplot2::theme(legend.position = 'none')
n <- export_ggchord_fusion(none, 8, 6)
stopifnot(n$metadata$guide_side == 'none', n$metadata$guide_layout$status == 'no_guides')
# The shared single-panel guide engine keeps its independent native contract.
single <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_link_ribbon() + adapt_ggchord()
s <- diagnose_ggchord(single, width = 8, height = 6, device = 'png')
stopifnot(s$guide_layout$status == 'feasible')
ggplot2::ggsave('tools/long-term-audit/v0.27.0-fusion-guides.png', make('adaptive'),
  width = 12, height = 7, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', compact_identity = g, multiple_guides = b$metadata$guide_layout$guides,
  full_source_records = TRUE, native_mappings_preserved = TRUE,
  guides_outside_panels = TRUE, explicit_constraints_fail_honestly = TRUE,
  absolute_boxes_preserved = TRUE, inside_rejected = TRUE, shared_single_panel_engine = TRUE),
  'tools/long-term-audit/v0.27.0-fusion-guide-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Fusion outer Guide budget: passed\n')
