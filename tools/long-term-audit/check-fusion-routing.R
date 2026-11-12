# Direct measured-routing checks. Source observations are never aggregated.
pkgload::load_all(quiet = TRUE)
library(dplyr)
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
data(seq_data_example); data(ribbon_data_example)
hit <- ribbon_data_example |> slice_head(n = 1)
base <- ggchord(seq_data_example) + geom_seq(show.legend = FALSE) +
  coord_collinear(seq_axis = 'vertical') + theme(ggchord.axis = element_blank())
frame <- ggchord_fusion(base, center_box = c(.02, .28, .12, .88)) +
  fusion_panel(base, 'target', box = c(.72, .98, .12, .88)) +
  fusion_panel(base, 'obstacle', box = c(.4, .6, .22, .78))
direct <- frame + fusion_ribbon(hit, qpanel = 'center', spanel = 'target', fill = '#548CA0')
routed <- frame + fusion_ribbon(hit, qpanel = 'center', spanel = 'target', fill = '#548CA0', route = 'avoid')
a <- export_ggchord_fusion(direct, 12, 7)
b <- export_ggchord_fusion(routed, 12, 7)
stopifnot(a$connections$route_status == 'direct_collision',
  b$connections$route_status == 'routed_clear', b$connections$obstacle_collisions == 0,
  identical(a$fronts, b$fronts), identical(b$relationship_data[[1]], hit))
obstacles <- ggchord_route_obstacles(b$routing_obstacles, c('center', 'target'))
stopifnot(!any(ggchord_route_polygon_hits(b$paths[c('x', 'y')], obstacles)))
# Endpoint departures retain the physical interval normal, not a centre anchor.
qn <- nrow(filter(b$fronts, endpoint == 'q')); sn <- nrow(filter(b$fronts, endpoint == 's'))
edge_n <- (nrow(b$paths)-qn-sn)/2
edge <- b$paths[qn + seq_len(edge_n), c('x', 'y')]
front <- filter(b$fronts, endpoint == 'q')
delta <- unlist(edge[2, ]-edge[1, ])
normal <- c(front$control_normal_x[1], front$control_normal_y[1])
stopifnot(sum(delta*normal) > 0, abs(delta[1]*normal[2]-delta[2]*normal[1]) < .001)
# Group three actual HSPs in a shared DISPLAY corridor, preserving every record.
hits <- ribbon_data_example |> filter(qaccver == hit$qaccver, saccver == hit$saccver)
hits$route_group <- 'explicit-shared-channel'
hits$source_hsp_id <- paste0('bundled-real-HSP-', seq_len(nrow(hits)))
bundled <- frame + fusion_ribbon(hits, qpanel = 'center', spanel = 'target',
  route = 'bundle', fill = '#548CA0')
c <- export_ggchord_fusion(bundled, 12, 7)
stopifnot(nrow(c$connections) == nrow(hits), all(c$connections$route_status == 'bundled_clear'),
  identical(c$relationship_data[[1]], hits), nrow(c$source_observations) == nrow(hits),
  all(c$source_observations$displayed_connections == 1))
corridors <- split(c$routing_corridors, c$routing_corridors$source_row)
stopifnot(all(vapply(corridors, function(d) identical(unname(as.matrix(d[c('vertex', 'x', 'y')])),
  unname(as.matrix(corridors[[1]][c('vertex', 'x', 'y')]))), logical(1))))
for (row in seq_len(nrow(hits))) stopifnot(!any(ggchord_route_polygon_hits(
  filter(c$paths, source_row == row)[c('x', 'y')], obstacles)))
# Subject decoration reuses the renderer's exact tip polygon for collision checks.
tipped <- frame + fusion_ribbon(hits |> slice_tail(n = 1), qpanel = 'center', spanel = 'target',
  route = 'avoid', gradient = TRUE, tip = 'subject', colour = '#325766', linewidth = .3)
t <- export_ggchord_fusion(tipped, 12, 7)
stopifnot(t$connections$obstacle_collisions == 0, t$connections$tip == 'subject')
ggplot2::ggsave('/tmp/ggchord-routing-tip.png', tipped, width = 12, height = 7, dpi = 90, bg = 'white')
# Decorative leaders use the same physical obstacle contract, no HSP records.
region <- tibble(accver = seq_data_example$accver[1], start = 1000, end = 1500)
focus <- magnify_ggchord(base, region, base, route = 'avoid',
  center_box = c(.02, .28, .12, .88), detail_box = c(.72, .98, .12, .88)) +
  fusion_panel(base, 'obstacle', box = c(.4, .6, .22, .78))
f <- export_ggchord_fusion(focus, 12, 7)
stopifnot(nrow(f$leaders) == 2, nrow(f$connections) == 0,
  all(f$leaders$obstacle_collisions == 0), nrow(f$leader_paths) > 4,
  identical(f$focus_region_data, region))
for (boundary in c('start', 'end')) {
  line <- filter(f$leader_paths, .data$boundary == .env$boundary)
  leader <- filter(f$leaders, .data$boundary == .env$boundary)
  stopifnot(identical(line$x[c(1, nrow(line))], c(leader$x, leader$xend)),
    identical(line$y[c(1, nrow(line))], c(leader$y, leader$yend)))
}
# A wall spans the full canvas; an explicit failure must replace a false clear result.
wall <- ggchord_fusion(base, center_box = c(.02, .28, .12, .88)) +
  fusion_panel(base, 'target', box = c(.72, .98, .12, .88)) +
  fusion_panel(base, 'wall', box = c(.4, .6, 0, 1)) +
  fusion_ribbon(hit, qpanel = 'center', spanel = 'target', route = 'avoid', fill = '#548CA0')
fail(export_ggchord_fusion(wall, 12, 7), 'no collision-free')
fail(fusion_ribbon(hit, route = 'bundle'), 'route_group')
# Parent additions must not leave a stale clear claim on a child route.
parent <- ggchord_fusion(routed, allow_overlap = TRUE, center_box = c(.02, .98, .02, .98)) +
  fusion_panel(base, 'cover', box = c(.32, .68, 0, 1))
fail(export_ggchord_fusion(parent, 12, 7), 'enclosing layout obstructs')
# Polygon enclosure counts even when no sampled centreline or edge cuts the box.
box <- data.frame(xmin = 4, xmax = 6, ymin = 4, ymax = 6)
stopifnot(ggchord_route_polygon_hits(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10)), box))
# Final native text/guide envelopes are remeasured with actual fonts.
labelled <- base + geom_seq_label(size = 3) + ggplot2::labs(title = 'Measured panel title')
measured <- export_ggchord_fusion(ggchord_fusion(labelled) + fusion_panel(labelled, 'detail'), 12, 7)
stopifnot(any(measured$routing_obstacles$kind == 'text'),
  any(measured$routing_obstacles$kind == 'decoration'))
# Rebuild at a different physical aspect ratio; no previous pixels are reused.
wider <- export_ggchord_fusion(routed, 16, 6)
stopifnot(wider$connections$obstacle_collisions == 0,
  identical(wider$relationship_data[[1]], hit))
example <- new.env()
source('examples/fusion-routing.R', local = example)
ggplot2::ggsave('tools/long-term-audit/v0.27.0-fusion-routing.png', example$routed,
  width = 12, height = 7, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', real_HSPs = nrow(hits),
  every_source_retained = TRUE, complete_fronts_unchanged = TRUE,
  direct_obstacle_reported = TRUE, routed_polygon_clear = TRUE,
  shared_corridor_without_aggregation = TRUE, focus_leaders_clear = TRUE,
  no_route_explicit_failure = TRUE, physical_rebuild = TRUE,
  scope = 'measured other containers, native layer text and gtable decoration/guide envelopes; no global optimum claim'),
  'tools/long-term-audit/v0.27.0-routing-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Fusion physical routing: passed\n')
