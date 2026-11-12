# Direct automatic-placement checks. Not a device/release matrix.
pkgload::load_all(quiet = TRUE)
source('examples/fusion-auto.R')
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
a <- export_ggchord_fusion(auto_fusion, 16, 12)
stopifnot(nrow(a$panels) == 13, nrow(a$connections) == 12,
  all(a$connections$retained), nrow(a$source_observations) == 1,
  a$source_observations$displayed_connections == 12,
  a$metadata$panel_layout$status == 'feasible',
  !a$metadata$panel_layout$placement_rescales_fonts)
for (i in seq_along(a$relationship_data)) stopifnot(identical(a$relationship_data[[i]], observed))
check_boxes <- function(d) {
  for (i in seq_len(nrow(d))) for (j in seq_len(i-1L)) stopifnot(
    min(d$container_xmax[i], d$container_xmax[j]) <= max(d$container_xmin[i], d$container_xmin[j]) ||
    min(d$container_ymax[i], d$container_ymax[j]) <= max(d$container_ymin[i], d$container_ymin[j]))
  stopifnot(all(d$container_xmin >= 0), all(d$container_ymin >= 0))
}
check_boxes(a$panels)
# Physical angle is relative to the actual central allocation on a rectangular canvas.
center <- filter(a$panels, panel_id == 'center')
children <- filter(a$panels, side == 'auto')
dx <- (children$container_xmin+children$container_xmax-center$container_xmin-center$container_xmax)/2
dy <- (children$container_ymin+children$container_ymax-center$container_ymin-center$container_ymax)/2
angle <- (atan2(dy, dx)*180/pi) %% 360
stopifnot(all(abs(((angle-children$auto_angle+180) %% 360)-180) < .01))
# An explicit fixed panel/absolute center and explicit auto size remain constraints.
plain <- ggchord(target_seqs) + geom_seq(show.legend = FALSE) + coord_collinear() +
  theme(ggchord.axis = element_blank())
mixed <- ggchord_fusion(plain, center_box = c(.4, .6, .4, .6)) +
  fusion_panel(plain, 'fixed', box = c(.78, .98, .35, .65)) +
  fusion_panel(plain, 'sized', side = 'auto', size = c(.12, .16)) +
  fusion_panel(plain, 'auto-2', side = 'auto') +
  fusion_panel(plain, 'auto-3', side = 'auto')
m <- export_ggchord_fusion(mixed, 12, 8)
check_boxes(m$panels)
fixed <- filter(m$panels, panel_id == 'fixed')
sized <- filter(m$panels, panel_id == 'sized')
stopifnot(abs(fixed$container_xmin/304.8-.78) < .001,
  abs((sized$container_xmax-sized$container_xmin)/304.8-.12) < .001,
  abs((sized$container_ymax-sized$container_ymin)/203.2-.16) < .001)
# Same object rebuilds on another aspect ratio, keeping every source row/front.
b <- export_ggchord_fusion(auto_fusion, 12, 8)
check_boxes(b$panels)
stopifnot(nrow(b$connections) == 12, identical(a$relationship_data, b$relationship_data))
# Magnification can route multiple independent details to repeated auto sides.
region <- tibble(accver = target_seqs$accver, start = 1000, end = 1500,
  detail_panel = c('one', 'two'))
focus <- magnify_ggchord(plain, region, list(one = plain, two = plain),
  side = c(one = 'auto', two = 'auto'))
f <- export_ggchord_fusion(focus, 10, 7)
stopifnot(nrow(f$focus_regions) == 2, nrow(f$leaders) == 4,
  all(f$panels$side[-1L] == 'auto'), identical(f$focus_region_data, region))
# Explicit sizes are never silently shrunk to force a packing success.
impossible <- ggchord_fusion(plain, center_box = c(.2, .8, .2, .8)) +
  fusion_panel(plain, 'too-large', side = 'auto', size = c(.8, .8))
fail(export_ggchord_fusion(impossible, 4, 4), 'automatic perimeter panels cannot fit')
ggplot2::ggsave('tools/long-term-audit/v0.27.0-fusion-auto.png', auto_fusion,
  width = 16, height = 12, dpi = 95, bg = 'white')
jsonlite::write_json(list(status = 'passed', independent_panels = 13,
  displayed_connections = 12, scientific_observations = 1,
  clockwise_physical_placement = TRUE, boxes_disjoint = TRUE,
  explicit_boxes_sizes_preserved = TRUE, source_rows_preserved = TRUE,
  multiple_auto_details = TRUE, alternate_aspect_ratio = TRUE,
  impossible_size_rejected = TRUE, placement = a$metadata$panel_layout),
  'tools/long-term-audit/v0.27.0-fusion-auto-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Fusion automatic perimeter: passed\n')
