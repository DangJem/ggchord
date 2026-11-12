# Measurement collector only; original renderer and source tables are retained.
source('tools/audit-source-snapshot.R')
snapshot <- audit_source_snapshot(c('tools/long-term-audit/observe-linear-refinement.R',
  'tools/long-term-audit/linear-lesson-source.R','tools/documentation-examples.R'))
pkgload::load_all(quiet=TRUE)
source('tools/long-term-audit/linear-lesson-source.R')
lesson <- load_linear_lesson()
original <- lesson$show_features
namespace <- asNamespace('ggchord')
original_context <- get('ggchord_panel_context',namespace)
events <- list()
collector <- function(table,built) {
  context <- original_context(table,built)
  events[[length(events)+1L]] <<- list(status=context$status,
    units_per_mm=context$units_per_mm,panel_mm=context$panel_mm)
  context
}
unlockBinding('ggchord_panel_context',namespace)
assign('ggchord_panel_context',collector,namespace)
lockBinding('ggchord_panel_context',namespace)
result <- diagnose_ggchord(lesson$p_linear_labels,width=10,height=6,dpi=144,device='png')
stopifnot(identical(original,lesson$show_features))
audit_assert_source_snapshot(snapshot)
jsonlite::write_json(list(scope='Native linear-label lesson at original 10x6 inches/144 dpi; measurement collector returns original contexts unchanged',
  events=events,measurement=result$measurement,checks=result$checks,warnings=result$warnings,
  label_quality=result$label_quality,source_features_preserved=TRUE,
  runtime_sha256=as.list(snapshot$engine)),
  'tools/long-term-audit/v0.27.0-linear-refinement-observations.json',pretty=TRUE,auto_unbox=TRUE,na='null',digits=12)
cat('Native linear refinement observations recorded\n')
