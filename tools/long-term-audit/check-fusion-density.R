# Source-backed dense composition: 39 distinct observed plant records, no
# fabricated scientific duplicates. One PNG, plus one alternate measured size.
pkgload::load_all(quiet = TRUE)
library(dplyr)
direct <- '--direct' %in% commandArgs(trailingOnly = TRUE)
route_mode <- if (direct) 'direct' else 'avoid'
after <- '--after' %in% commandArgs(trailingOnly = TRUE)
prefix <- paste0('tools/long-term-audit/v0.27.0-fusion-density', if (direct) '-direct' else '', if (after) '-after' else '')
input <- read_delta('examples/genome-alignments/arabidopsis-first-pair.delta',
  executable = Sys.getenv('DELTA_DECODER', 'show-coords'), source_file = TRUE)
records <- input$ribbon_data |>
  mutate(source_hsp_id = paste0('plant-delta-physical-line-', .source_line))
stopifnot(nrow(records) == 39L, !anyDuplicated(records$source_hsp_id))
subject <- input$seq_data |> filter(accver == first(records$saccver))
base <- ggchord(subject) + geom_seq(show.legend = FALSE) +
  theme(ggchord.axis = element_blank())
center <- ggchord(input$seq_data) + geom_seq(show.legend = FALSE) +
  geom_seq_label(size = 2.5) + coord_chord() +
  theme(ggchord.axis = element_blank())
targets <- list(horizontal = base + coord_collinear(),
  vertical = base + coord_collinear(seq_axis = 'vertical'),
  mirror = base + coord_sankey())
boxes <- list(horizontal = c(.65, .98, .70, .98), vertical = c(.65, .98, .36, .64),
  mirror = c(.65, .98, .02, .30))
fusion <- ggchord_fusion(center, center_box = c(.02, .42, .15, .85))
declared <- list()
for (i in seq_along(targets)) {
  id <- names(targets)[i]
  target <- targets[[i]] + labs(title = id)
  part <- records |> slice(seq.int(i, nrow(records), by = length(targets)))
  if (id == 'mirror') {
    instances <- get_chord_layout(target)$display_instances
    part$s_instance <- instances$instance_id[instances$accver == subject$accver & instances$side == 'right']
  }
  declared[[i]] <- part
  fusion <- fusion + fusion_panel(target, id, box = boxes[[id]]) +
    fusion_ribbon(part, qpanel = 'center', spanel = id, route = route_mode,
      fill = '#548CA0', alpha = .35, show.legend = FALSE)
}
measure <- function(width, height, file = NULL) {
  if (is.null(file)) {
    close <- ggchord_measurement_device(width, height)
  } else {
    ggchord_diagnostic_device(file, 'png', width, height, 95)
    close <- function() grDevices::dev.off()
  }
  on.exit(close())
  grid::grid.newpage()
  value <- tryCatch(ggchord_build_fusion(fusion), error = identity)
  if (inherits(value, 'error')) return(list(status = 'explicit_failure',
    width = width, height = height, message = conditionMessage(value)))
  grid::grid.draw(value$grob)
  exported <- value$export
  stopifnot(nrow(exported$connections) == 39L, all(exported$connections$retained),
    nrow(exported$source_observations) == 39L,
    all(exported$source_observations$displayed_connections == 1L))
  for (i in seq_along(declared)) stopifnot(identical(exported$relationship_data[[i]], declared[[i]]))
  list(status = 'measured', width = width, height = height,
    connections = exported$connections, observations = exported$source_observations,
    panel_boxes = exported$panels, route_counts = as.list(table(exported$connections$route_status)),
    obstacle_contacts = sum(exported$connections$obstacle_collisions),
    source_rows_retained = TRUE, explicit_mirror_instance = TRUE,
    routing_obstacles = exported$routing_obstacles,
    panel_measurements = exported$diagnostics)
}
results <- list(landscape = measure(16, 10, paste0(prefix, '.png')),
  alternate = measure(10, 8))
jsonlite::write_json(list(scope = '39 distinct original records split across three independent views; no cross-panel aggregation', route = route_mode,
  runtime_sha256 = as.list(stats::setNames(vapply(list.files('R', '[.]R$', full.names = TRUE),
    function(path) digest::digest(file = path, algo = 'sha256'), character(1)), list.files('R', '[.]R$', full.names = TRUE))),
  results = results, unverified = c('connection versus connection crossings/readability',
    'exact ink, transparent occlusion and all panel label-feature contacts', 'release device/platform matrix')),
  paste0(prefix, '-results.json'), pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Dense Fusion outcomes:', paste(vapply(results, `[[`, '', 'status'), collapse = ', '), '\n')
