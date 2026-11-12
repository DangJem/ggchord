# Native plot furniture surrounds the biological composition; no dummy axes.
ggchord_fusion_has_decoration <- function(x) {
  any(vapply(x$styles, function(style) {
    if (inherits(style, 'theme')) return(any(grepl('^(plot[.]|panel[.]background$|aspect[.]ratio$)', names(style))))
    inherits(style, 'ggplot2::labels') && any(names(style) %in% c('title', 'subtitle', 'caption', 'tag'))
  }, logical(1)))
}

ggchord_fusion_decoration <- function(x) {
  p <- ggplot2::ggplot() + ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  for (style in x$styles) if (inherits(style, 'theme') || inherits(style, 'ggplot2::labels')) p <- p + style
  p <- p + ggplot2::theme(legend.position = 'none')
  # Complete user themes must not reveal the carrier's meaningless 0--1 axes.
  for (name in grep('^axis[.]', names(p$theme), value = TRUE))
    if (!grepl('length', name)) p$theme[[name]] <- ggplot2::element_blank()
  p$theme$panel.grid <- p$theme$panel.grid.major <- p$theme$panel.grid.minor <- ggplot2::element_blank()
  built <- ggplot2::ggplot_build(p)
  table <- ggplot2::ggplot_gtable(built)
  context <- ggchord_panel_context(table, built)
  if (!identical(context$status, 'measured') || min(context$panel_mm) < 15)
    ggchord_stop('Fusion: native plot decorations and margins leave insufficient composition space; enlarge the canvas or revise explicit plot theme')
  background <- ggchord_fusion_partition(table, TRUE)
  foreground <- ggchord_fusion_partition(table, FALSE)
  at <- which(grepl('^panel($|[-.])', table$layout$name))
  for (i in at) foreground$grobs[[i]] <- grid::nullGrob()
  list(background = background, foreground = foreground, context = context,
    labels = built$plot$labels, table = table)
}

ggchord_build_fusion <- function(x) {
  if (!ggchord_fusion_has_decoration(x)) return(ggchord_build_fusion_core(x))
  previous <- ggplot2::last_plot()
  on.exit(ggplot2::set_last_plot(previous))
  canvas <- c(width = grid::convertWidth(grid::unit(1, 'npc'), 'mm', TRUE),
    height = grid::convertHeight(grid::unit(1, 'npc'), 'mm', TRUE))
  origin <- grid::deviceLoc(grid::unit(0, 'npc'), grid::unit(0, 'npc'), valueOnly = TRUE)
  origin <- c(origin$x, origin$y)*25.4
  decoration <- ggchord_fusion_decoration(x)
  bounds <- decoration$context$panel_bounds_mm - rep(origin, each = 2L)
  dimensions <- stats::setNames(c(diff(bounds[1:2]), diff(bounds[3:4])), c("width", "height"))
  offset <- bounds[c(1, 3)]
  viewport <- grid::viewport(x = mean(bounds[1:2]), y = mean(bounds[3:4]),
    width = dimensions[1], height = dimensions[2], default.units = 'mm', clip = 'off')
  inner <- x
  inner$styles <- c(inner$styles, list(ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL, tag = NULL),
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0), aspect.ratio = NULL)))
  result <- local({
    grid::pushViewport(viewport)
    on.exit(grid::popViewport())
    ggchord_with_render_context(list(available_canvas_mm = dimensions), ggchord_build_fusion_core(inner))
  })
  result$export <- ggchord_fusion_nested_export(result$export, NULL, offset)
  for (id in names(result$panel_registry)) result$panel_registry[[id]]$context$panel_bounds_mm <-
    result$panel_registry[[id]]$context$panel_bounds_mm + rep(offset, each = 2L)
  for (id in names(result$export$diagnostics)) {
    d <- result$export$diagnostics[[id]]
    if (!is.null(d$render_context$panel_bounds_mm)) d$render_context$panel_bounds_mm <-
      d$render_context$panel_bounds_mm + rep(offset, each = 2L)
    if (isTRUE(d$nested)) {
      if (length(d$content_bounds_mm) == 4L) d$content_bounds_mm <- d$content_bounds_mm + rep(offset, each = 2L)
      if (!is.null(d$plot_decoration)) d$plot_decoration$bounds_mm <-
        ggchord_fusion_translate_table(d$plot_decoration$bounds_mm, offset, c('xmin', 'xmax'), c('ymin', 'ymax'))
      d$composition_guide_bounds_mm <- ggchord_fusion_translate_table(d$composition_guide_bounds_mm,
        offset, c('xmin', 'xmax'), c('ymin', 'ymax'))
      d$guide_layout$bounds <- ggchord_fusion_translate_table(d$guide_layout$bounds,
        offset, c('xmin_mm', 'xmax_mm'), c('ymin_mm', 'ymax_mm'))
    }
    result$export$diagnostics[[id]] <- d
  }
  events <- lapply(result$render_events, function(event) {
    event$background <- grid::grobTree(event$background, vp = viewport)
    event$foreground <- grid::grobTree(event$foreground, vp = viewport)
    if (!is.null(event$bounds)) event$bounds <- event$bounds + rep(offset, each = 2L)
    if (!is.null(event$container_bounds)) event$container_bounds <- event$container_bounds + rep(offset, each = 2L)
    event
  })
  # Unique wrappers prevent two native gtables named "layout" from aliasing
  # their children in a shared grid drawing list.
  events <- c(list(list(panel_id = NULL, background = grid::grobTree(decoration$background),
    foreground = grid::nullGrob(), bounds = bounds, container_bounds = c(0, canvas[1], 0, canvas[2]),
    composition_background = TRUE)), events, list(list(panel_id = NULL,
      background = grid::nullGrob(), foreground = grid::grobTree(decoration$foreground), bounds = NULL)))
  connection_viewport <- viewport; connection_viewport$clip <- grid::viewport(clip = "inherit")$clip
  connections <- lapply(result$connection_grobs, function(g) grid::grobTree(g, vp = connection_viewport))
  overlays <- lapply(result$overlay_grobs, function(g) grid::grobTree(g, vp = viewport))
  obstacles <- ggchord_fusion_guide_obstacles(decoration$table, '', origin)
  result$export$routing_obstacles <- ggchord_rbind_fill(list(result$export$routing_obstacles, obstacles))
  result$export <- ggchord_fusion_recheck_routes(result$export, result$export$routing_obstacles)
  result$export$metadata$content_canvas_mm <- dimensions
  result$export$metadata$content_bounds_mm <- bounds
  result$export$metadata$canvas_mm <- canvas
  result$export$metadata$plot_decoration <- list(status = 'native_gtable',
    labels = decoration$labels[intersect(names(decoration$labels), c('title', 'subtitle', 'caption', 'tag'))],
    bounds_mm = obstacles, coordinate_axes = 'none; child axes belong to each child plot')
  result$render_events <- events; result$connection_grobs <- connections; result$overlay_grobs <- overlays
  result$grob <- ggchord_fusion_render_events(events, connections, overlays, canvas,
    isTRUE(x$allow_overlap) || nrow(result$export$overlaps) > 0)
  result
}
