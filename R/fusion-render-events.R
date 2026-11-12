# One rendering ledger for ordinary, overlapping and decorated compositions.
ggchord_fusion_render_events <- function(render_events, connection_grobs, overlay_grobs, canvas, overlap_mode) {
  if (overlap_mode) {
    containers <- Filter(Negate(is.null), lapply(render_events, `[[`, "container_bounds"))
    clip_regions <- function(regions) lapply(regions, function(bounds)
      ggchord_fusion_clip_connections(connection_grobs, bounds, canvas))
    children <- clip_regions(ggchord_fusion_rect_subtract(c(0,canvas[1],0,canvas[2]), containers))
    for (i in seq_along(render_events)) {
      event <- render_events[[i]]
      children <- c(children, list(event$background))
      if (!is.null(event$bounds) && (!is.null(event$panel_id) || isTRUE(event$composition_background))) {
        later <- if (i < length(render_events)) Filter(Negate(is.null),
          lapply(render_events[seq.int(i+1L,length(render_events))], `[[`, "container_bounds")) else list()
        children <- c(children, clip_regions(ggchord_fusion_rect_subtract(event$bounds, later)))
      }
      children <- c(children, list(event$foreground),
        unname(overlay_grobs[names(overlay_grobs) == event$panel_id]))
    }
  } else children <- c(lapply(render_events, `[[`, "background"), connection_grobs,
    lapply(render_events, `[[`, "foreground"), unname(overlay_grobs))
  grob <- do.call(grid::grobTree, c(unname(children), list(vp = grid::viewport(xscale=c(0,1),yscale=c(0,1)))))
  grob
}
