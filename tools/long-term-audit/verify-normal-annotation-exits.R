# Verify actual normal-size renders and the explicit-style ownership guard.
# No package checks, redraw matrix or changes to scientific inputs.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot()
pkgload::load_all(quiet = TRUE)
folder <- 'tools/long-term-audit'
parts <- lapply(c('a', 'b'), function(part) jsonlite::read_json(
  file.path(folder, paste0('v0.27.0-all-exits-', part, '-results.json'))))
stopifnot(identical(parts[[1L]]$runtime_sha256, parts[[2L]]$runtime_sha256))
render_engine <- parts[[1L]]$runtime_sha256
current_engine <- as.list(audit$engine)
stopifnot(setequal(names(render_engine), names(current_engine)))
changed <- names(current_engine)[!vapply(names(current_engine), function(path)
  identical(current_engine[[path]], render_engine[[path]]), logical(1))]
stopifnot(identical(changed, 'R/annotation-exits.R'))
# The only post-render runtime edit excludes straight/elbow styles from the
# automatic route search. Prove the previous source hash instead of upgrading it.
path <- 'R/annotation-exits.R'
current <- readChar(path, file.info(path)$size, useBytes = TRUE)
guard <- paste0("      # Explicit straight/elbow styles retain their already chosen routes.\n",
  "      if ('.leader_style' %in% names(labels)) labels <- labels[\n",
  "        !labels$.leader_style %in% c('straight','elbow'),,drop=FALSE]\n")
stopifnot(grepl(guard, current, fixed = TRUE))
previous <- sub(guard, '', current, fixed = TRUE)
stopifnot(identical(digest::digest(charToRaw(previous), algo = 'sha256', serialize = FALSE),
  render_engine[[path]]))
environment <- new.env(parent = asNamespace('ggchord'))
eval(parse(text = previous), environment)
previous_solver <- environment$ggchord_clear_simple_annotation_exits
cases <- do.call(c, lapply(parts, `[[`, 'cases'))
baseline <- jsonlite::read_json(file.path(folder, 'v0.27.0-compact-final-results.json'))
stopifnot(!anyDuplicated(names(cases)), setequal(names(cases), names(baseline$cases)))
rows <- list()
retention_path <- file.path(folder, 'v0.27.0-normal-image-retention.json')
retention <- if (file.exists(retention_path)) jsonlite::read_json(retention_path) else NULL
retained_image <- function(path) {
  reused <- retention$reused_images[[basename(path)]]
  if (is.null(reused)) return(path)
  retained <- file.path(folder, reused$file)
  stopifnot(identical(digest::digest(file = retained, algo = 'sha256'), reused$sha256))
  retained
}
for (id in names(baseline$cases)) {
  part <- if (id %in% names(parts[[1L]]$cases)) 'a' else 'b'
  original <- readRDS(paste0('/tmp/v0.27.0-compact-final-', id, '.rds'))
  layout <- readRDS(paste0('/tmp/v0.27.0-all-exits-', part, '-', id, '.rds'))
  stopifnot(identical(original$layer_inputs, layout$layer_inputs))
  for (owner in layout$layer_geometry) {
    data <- owner$restriction_site
    if (is.data.frame(data) && nrow(data)) stopifnot(all(data$.leader_style %in% 'radial'))
  }
  stopifnot(identical(previous_solver(layout$layer_geometry, layout),
    ggchord_clear_simple_annotation_exits(layout$layer_geometry, layout)))
  case <- cases[[id]]; q <- case$quality
  stopifnot(isTRUE(case$converged), !length(case$warnings),
    all(vapply(case$visibility, function(row) isTRUE(row$visible), logical(1))),
    q$n_text_overlaps == 0, q$n_leader_crossings == 0,
    q$n_leaders_through_text == 0, q$n_label_outside_canvas == 0,
    case$association$anchors_off_own_feature == 0,
    case$association$leaders_over_one_bend == 0)
  rows[[id]] <- list(source_labels = length(case$visibility), complete_inputs_identical = TRUE,
    default_radial_style = TRUE, post_render_guard_geometry_equivalent = TRUE,
    file = case$file, retained_file = retained_image(case$file),
    png_sha256 = digest::digest(file = retained_image(case$file), algo = 'sha256'))
}
normal <- jsonlite::read_json(file.path(folder, 'v0.27.0-normal-flagship-matched-results.json'))
stopifnot(identical(normal$runtime_sha256, current_engine))
q <- normal$profiles$matched
stopifnot(q$n_visible_source_labels == q$n_source_labels, q$n_text_overlaps == 0,
  q$n_leader_crossings == 0, q$n_leaders_through_text == 0,
  q$n_label_outside_canvas == 0, isTRUE(q$geometry_converged))
explicit <- jsonlite::read_json(file.path(folder, 'v0.27.0-explicit-exit-style-results.json'))
stopifnot(identical(explicit$runtime_sha256, current_engine),
  all(unlist(explicit$explicit_styles_keep_routes)))
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(cases = cases, runtime_sha256 = render_engine,
  input_sha256 = parts[[1L]]$input_sha256,
  scope = 'Actual 13 complete-source 10x8in/110dpi renders before the explicit-style guard; not relabelled as new-engine device runs'),
  file.path(folder, 'v0.27.0-all-exits-results.json'), pretty = TRUE, auto_unbox = TRUE)
jsonlite::write_json(list(status = 'passed normal-size source and measured-layout controls',
  cases = rows, source_labels = sum(vapply(rows, `[[`, integer(1), 'source_labels')),
  rendered_runtime_sha256 = render_engine, current_runtime_sha256 = current_engine,
  guarded_styles = c('straight', 'elbow'), only_runtime_delta = changed,
  current_native_flagship = q, current_explicit_style_control = explicit,
  normal_sizes = '13 maps: 10x8in; complete flagship: 10x8in and 16.5x8in; no new extreme-size gate',
  not_checked = c('arbitrary canvases or custom grobs', 'exact glyph ink, occlusion and font portability',
    'publication package/device/platform checks')),
  file.path(folder, 'v0.27.0-normal-layout-verification.json'), pretty = TRUE, auto_unbox = TRUE)
cat('13 actual source-preserving renders + exact guard delta + current normal flagship: passed\n')
