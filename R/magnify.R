# Focus+context uses the same panel/frame projection as biological links, but
# its leaders and highlights never enter the relationship/HSP pipeline.

#' Compose a source region and an explicit local detail plot
#'
#' Marks real source intervals in registered genomic panels and connects them to
#' the same intervals in independent genomic or registered character details.
#' A detail or source can itself be a nested Fusion.
#' Returns a Fusion, not a new Coord or a crop of an existing image. Ordinary
#' ggsave rebuilds both panels, highlights and leaders at the actual target size.
#'
#' Build `detail_plot` explicitly with [focus_ggchord_data()] or other documented
#' source-aware operations. This function does not copy/filter arbitrary layers
#' or silently choose a local/global Stat policy. The supplied detail plot owns
#' its data, mappings, statistics and labels. Original source origins, including
#' renamed/repeated focus windows, are used to resolve the same interval.
#'
#' @param plot A registered ggchord genomic source plot or nested Fusion.
#' @param region Nonempty data frame with `accver`, `start`, `end` in original
#'   one-based inclusive source coordinates. Bounds must be ordered positive
#'   integers. Optional `instance` and `detail_instance` select exact display
#'   instances in the source and detail panels for every row. Ambiguity fails.
#'   Explicit wrap = TRUE permits start > end across a circular source origin.
#'   Both selected leaves must declare matching full source lengths; character
#'   sources require as_ggchord_alignment(circular = TRUE). The interval is
#'   drawn as two pieces with only its two external boundary leaders. With a list of details,
#'   a `detail_panel` column selects the named target for each row. Nested leaves
#'   use full paths, for example detail/center; optional source_panel similarly
#'   selects a source leaf under center. Ambiguous container selection fails.
#' @param detail_plot An explicit registered ggchord genomic plot, nested Fusion,
#'   a character plot from [register_ggchord_alignment()], or a uniquely named
#'   list of such plots. Every panel must have
#'   at least one source region; the name `center` is reserved.
#'   The selected source interval must be fully visible in both panels, in the
#'   same bp/aa unit. Character plots require explicit source/column registration;
#'   ordinary unregistered plots are rejected.
#' @param side A named side from fusion_panel(), including diagonal sides.
#'   For a list of plots, supply a character vector named by panel ID, with
#'   distinct sides. Multiple panels require explicit routing and placement.
#' @param center_size,detail_size Optional panel width/height fractions, following
#'   [ggchord_fusion()] and [fusion_panel()]. Defaults reserve a smaller source
#'   panel and a wider detail container for a single detail plot. With multiple
#'   details, use a named list of sizes or NULL for Fusion defaults; center_size
#'   also follows Fusion defaults. Out-of-frame allocations fail; overlap
#'   requires explicit allow_overlap.
#' @param detail_box Explicit detail allocation instead of side/detail_size,
#'   using the box convention in fusion_panel(). For multiple plots, supply a
#'   named list covering every detail. Physical boxes are rechecked at target size.
#' @param center_box Explicit source allocation instead of center_size.
#' @param allow_overlap Deliberately permit overlapping containers, following
#'   ggchord_fusion(). Later panels cover earlier ones; default FALSE.
#' @param connector `"leader"` draws two decorative boundary leaders per row; `"none"`
#'   retains the highlights only. These are not biological ribbons or alignments.
#' @param route Direct boundary leaders (default), or avoid to route around
#'   measured other containers and native text/guide envelopes. Endpoints stay
#'   fixed; no clear route fails rather than hiding the selection.
#' @param colour Fixed highlight/leader colour.
#' @param alpha Highlight opacity in `[0,1]`. Leaders use the opaque colour.
#' @param highlight_width Positive highlight band width in physical mm.
#' @param linewidth Nonnegative leader width using the usual ggplot2 convention.
#' @param linetype Leader line type, default dashed.
#' @return A [ggchord_fusion()] object. [export_ggchord_fusion()] additionally
#'   returns focus_regions, focus_fronts and leaders in mm, plus original
#'   focus_region_data (without metadata-name collisions), preserving the input
#'   region, exact instances, source origins and sampled physical path-length ratio.
#'   The derived tables use region_row to identify the original input row and
#'   region_piece to identify the origin-split piece. The raw input stays intact.
#'   Stat scope is recorded as explicit_detail_plot. Leaders do not appear in
#'   connections or relationship_data. leader_paths contains the final decorative
#'   polyline vertices. Routing is deterministic in measured physical space, not
#'   a guarantee of globally optimal layout or automatic layer rebuilding. The measured ratio may be below one if the chosen
#'   detail allocation does not enlarge the region.
#' @export
#' @examples
#' library(dplyr)
#' data(plasmid_example_pUC19)
#' whole <- ggchord(plasmid_example_pUC19) + geom_seq() + coord_circular()
#' region <- tibble(accver = "pUC19", start = 101, end = 1900)
#' local <- focus_ggchord_data(plasmid_example_pUC19, loci = region)
#' detail <- ggchord(local$seq_data) +
#'   geom_seq() + coord_collinear(seq_axis = "vertical")
#' magnify_ggchord(whole, region, detail)
#' magnify_ggchord(whole, region, detail,
#'   center_box = c(.02, .35, .2, .8), detail_box = c(.45, .98, .1, .9))
magnify_ggchord <- function(plot, region, detail_plot,
    side = c("right", "top", "left", "bottom"), center_size = NULL, detail_size = NULL,
    connector = c("leader", "none"), colour = "#B45309", alpha = .2,
    highlight_width = 1.5, linewidth = .35, linetype = 2, detail_box = NULL, allow_overlap = FALSE, center_box = NULL, route = c("direct", "avoid")) {
  route <- match.arg(route)
  side_supplied <- !missing(side)
  ggchord_fusion_check_plot(plot, center = TRUE)
  caller <- "magnify_ggchord()"
  ggchord_require_columns(region, c("accver", "start", "end"), caller)
  wrap <- region[["wrap"]] %||% rep(FALSE, nrow(region))
  if (!is.logical(wrap) || anyNA(wrap)) ggchord_stop(caller, ": wrap must contain logical declarations")
  if (nrow(region) < 1L || !is.character(region$accver) || anyNA(region$accver) || any(!nzchar(region$accver)) ||
      !is.numeric(region$start) || !is.numeric(region$end) ||
      any(!is.finite(c(region$start,region$end))) || any(region$start < 1) || any(ifelse(wrap, region$end >= region$start, region$end < region$start)) ||
      any(c(region$start,region$end) != floor(c(region$start,region$end))) || any(c(region$start,region$end) > 2^53-1))
    ggchord_stop(caller, ": provide ordered positive integral intervals, or start > end with explicit wrap = TRUE")
  for (name in intersect(c("instance", "detail_instance", "source_panel"), names(region)))
    if (!is.character(region[[name]]) || anyNA(region[[name]]) || any(!nzchar(region[[name]])))
      ggchord_stop(caller, ": instance selectors must be nonempty strings")
  if (!is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<0 || alpha>1 ||
      !is.numeric(highlight_width) || length(highlight_width)!=1L || !is.finite(highlight_width) || highlight_width<=0 ||
      !is.numeric(linewidth) || length(linewidth)!=1L || !is.finite(linewidth) || linewidth<0)
    ggchord_stop(caller, ": alpha must be in [0,1], highlight_width positive and linewidth nonnegative")
  if (!is.character(colour) || length(colour)!=1L || is.na(colour))
    ggchord_stop(caller, ": colour must be a single colour string")
  tryCatch(grDevices::col2rgb(colour),error=function(e) ggchord_stop(caller, ": invalid colour"))
  if (length(linetype)!=1L || is.na(linetype)) ggchord_stop(caller, ": linetype must be scalar")
  multiple <- !inherits(detail_plot, "ggplot")
  details <- if (multiple) detail_plot else list(detail = detail_plot)
  ids <- names(details)
  if (!is.list(details) || !length(details) ||
      is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || "center" %in% ids)
    ggchord_stop(caller, ": detail_plot must be a plot or a uniquely named list of plots")
  lapply(details, ggchord_fusion_check_plot)
  boxes <- NULL
  if (!is.null(detail_box)) {
    if (side_supplied || !is.null(detail_size)) ggchord_stop(caller, ": use detail_box instead of side/detail_size")
    boxes <- if (multiple) detail_box else stats::setNames(list(detail_box), ids)
    if (!is.list(boxes) || is.null(names(boxes)) || anyDuplicated(names(boxes)) || !setequal(names(boxes), ids))
      ggchord_stop(caller, ": detail_box must match every detail panel")
    lapply(boxes, ggchord_fusion_box)
    sides <- stats::setNames(rep("free", length(ids)), ids)
  } else if (multiple) {
    if (!is.character(side) || is.null(names(side)) || anyDuplicated(names(side)) ||
        !setequal(names(side), ids) || anyNA(side) ||
        any(!side %in% ggchord_fusion_sides()) || anyDuplicated(side))
      ggchord_stop(caller, ": provide one distinct named side for each detail panel")
    sides <- side[ids]
  } else sides <- stats::setNames(if (side_supplied) match.arg(side, ggchord_fusion_sides()) else "right", ids)
  targets <- region[["detail_panel"]]
  if (is.null(targets)) {
    if (length(details) > 1L) ggchord_stop(caller, ": multiple details require a detail_panel column")
    targets <- rep(ids, nrow(region))
  }
  routes <- unlist(lapply(ids, function(id) ggchord_fusion_panel_paths(details[[id]], id)), use.names = FALSE)
  if (!is.character(targets) || anyNA(targets) || any(!targets %in% routes))
    ggchord_stop(caller, ": detail_panel must select an existing detail plot")
  used <- vapply(ids, function(id) any(targets == id | startsWith(targets, paste0(id, "/"))), logical(1))
  if (any(!used)) ggchord_stop(caller, ": every detail panel must have a source region")
  sizes <- if (multiple && !is.null(detail_size)) {
    if (!is.list(detail_size) || is.null(names(detail_size)) || anyDuplicated(names(detail_size)) ||
        !setequal(names(detail_size), ids))
      ggchord_stop(caller, ": detail_size must be a named list matching the detail plots")
    detail_size[ids]
  } else stats::setNames(rep(list(detail_size),length(ids)),ids)
  if (multiple || !is.null(boxes)) {
    result <- ggchord_fusion(plot, center_size = center_size, allow_overlap = allow_overlap, center_box = center_box)
    for (id in ids) result <- result + if (!is.null(boxes)) fusion_panel(details[[id]], id, box = boxes[[id]]) else
      fusion_panel(details[[id]], id, sides[[id]], sizes[[id]])
  } else {
    side <- sides[[1L]]
    lateral <- side %in% c("right", "left")
    if (is.null(center_size) && is.null(center_box)) center_size <- if (lateral) c(.4,.8) else c(.8,.4)
    if (is.null(detail_size)) detail_size <- if (lateral) c(.42,.88) else c(.88,.42)
    result <- ggchord_fusion(plot, center_size = center_size, allow_overlap = allow_overlap, center_box = center_box) +
      fusion_panel(detail_plot, "detail", side = side, size = detail_size)
  }
  resolved_region <- as.data.frame(region)
  resolved_region$detail_panel <- targets
  resolved_region$wrap <- wrap
  result$focus <- list(input = region, region = resolved_region, connector = match.arg(connector),
    colour = colour, alpha = alpha, width = highlight_width, linewidth = linewidth, linetype = linetype, route = route)
  ggplot2::set_last_plot(result)
  result
}

# Build only after the real panel sizes, position scales and Coord are final.
ggchord_build_focus <- function(focus, panels) {
  empty <- list(input = data.frame(), underlay = list(), overlay = list(), regions = data.frame(), fronts = data.frame(), leaders = data.frame(), leader_paths = data.frame())
  if (is.null(focus)) return(empty)
  pieces <- list()
  for (i in seq_len(nrow(focus$region))) {
    region <- focus$region[i, , drop = FALSE]
    if (isTRUE(region$wrap)) {
      # Resolve both endpoints to the SAME exact leaf/instance before splitting.
      resolve <- function(position, panel, instance) ggchord_resolve_panel_interval(
        panels, region$accver, position, position, panel, instance, "magnify_ggchord() wrap: ")
      q <- resolve(region$start, region[["source_panel"]] %||% "center", region[["instance"]])
      qe <- resolve(region$end, q$panel, q$instance)
      d <- resolve(region$start, region$detail_panel, region[["detail_instance"]])
      de <- resolve(region$end, d$panel, d$instance)
      source_length <- ggchord_focus_wrap_length(panels[[q$panel]], q)
      detail_length <- ggchord_focus_wrap_length(panels[[d$panel]], d)
      if (!is.finite(source_length) || !is.finite(detail_length) || source_length != detail_length)
        ggchord_stop("magnify_ggchord(): wrap requires matching declared full source lengths in source and detail")
      if (max(region$start, region$end) > source_length)
        ggchord_stop("magnify_ggchord(): wrapped endpoints exceed declared source length")
      region$source_panel <- q$panel; region$instance <- q$instance
      region$detail_panel <- d$panel; region$detail_instance <- d$instance
      bounds <- list(c(region$start, source_length), c(1, region$end))
    } else bounds <- list(c(region$start, region$end))
    for (j in seq_along(bounds)) {
      part <- region; part$start <- bounds[[j]][1L]; part$end <- bounds[[j]][2L]
      boundaries <- if (length(bounds) == 1L) c("start", "end") else if (j == 1L) "start" else "end"
      piece <- ggchord_build_focus_region(focus, panels, part, boundaries)
      for (name in c("regions", "fronts", "leaders", "leader_paths")) {
        piece[[name]]$region_row <- rep.int(i, nrow(piece[[name]]))
        piece[[name]]$region_piece <- rep.int(j, nrow(piece[[name]]))
      }
      piece$regions$wrap <- region$wrap
      pieces[[length(pieces) + 1L]] <- piece
    }
  }
  bind <- function(name) do.call(rbind, lapply(pieces, `[[`, name))
  grobs <- function(name) unlist(lapply(pieces, `[[`, name), recursive = FALSE)
  list(input = focus$input, underlay = grobs("underlay"), overlay = grobs("overlay"),
    regions = bind("regions"), fronts = bind("fronts"), leaders = bind("leaders"), leader_paths = bind("leader_paths"))
}

# Full source lengths belong to provenance, not the displayed local span.
ggchord_focus_wrap_length <- function(panel, endpoint) {
  if (isTRUE(panel$layout$character)) {
    d <- panel$layout$character_data
    d <- d[as.character(d$accver) == endpoint$instance, , drop = FALSE]
    if (!"source_circular" %in% names(d) || !isTRUE(all(d$source_circular)) ||
        !"source_length" %in% names(d)) return(NA_real_)
    n <- unique(d$source_length)
  } else {
    source <- ggchord_focus_source(panel$seq_data)
    # Repeated focused windows must agree about their complete source length.
    n <- unique(source$length[source$accver == endpoint$entity])
  }
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 1 || n != floor(n) || n > 2^53 - 1)
    return(NA_real_)
  n
}

ggchord_build_focus_region <- function(focus, panels, region, boundaries = c("start", "end")) {
  resolve <- function(panel, instance) ggchord_resolve_panel_interval(panels,
    region$accver, region$start, region$end, panel, instance,
    paste0("magnify_ggchord() ",panel,": "))
  q <- resolve(region[["source_panel"]] %||% "center", region[["instance"]])
  s <- resolve(region[["detail_panel"]] %||% "detail", region[["detail_instance"]])
  if (q$unit != s$unit) ggchord_stop("magnify_ggchord(): source and detail units must agree")
  qf <- ggchord_fusion_front(panels[[q$panel]], q)
  sf <- ggchord_fusion_front(panels[[s$panel]], s)
  grobs <- list()
  ends <- c(1L,nrow(qf))
  leaders <- data.frame(); leader_paths <- data.frame()
  if (focus$connector == "leader") {
    leaders <- data.frame(boundary = c("start", "end"),
      source_position_boundary = qf$position_boundary[ends],
      x = qf$x[ends],y = qf$y[ends],xend = sf$x[c(1L,nrow(sf))],yend = sf$y[c(1L,nrow(sf))],
      coordinate_space = "composition",units = "mm",semantic = "decorative_focus_leader")
    leaders <- leaders[leaders$boundary %in% boundaries, , drop = FALSE]
    obstacles <- ggchord_route_inflate(ggchord_route_obstacles(focus$obstacles, c(q$panel, s$panel)), focus$linewidth/2)
    paths <- lapply(seq_len(nrow(leaders)), function(i) {
      a <- c(leaders$x[i], leaders$y[i]); b <- c(leaders$xend[i], leaders$yend[i])
      hits <- sum(ggchord_route_segment_hits(a, b, obstacles))
      route <- focus$route %||% "direct"
      line <- if (route == "avoid") ggchord_route_shortest(a, b, obstacles, focus$canvas,
        clearance = focus$linewidth/2 + .6) else rbind(a, b)
      if (is.null(line)) ggchord_stop("magnify_ggchord(): no clear leader route; move/enlarge panels or inspect route = 'direct'")
      leaders$route[i] <<- route
      leaders$route_status[i] <<- if (nrow(line) > 2L) "routed_clear" else if (hits) "direct_collision" else "direct_clear"
      leaders$obstacle_collisions[i] <<- if (route == "avoid") 0L else hits
      leaders$length_mm[i] <<- sum(sqrt(rowSums((line[-1L, , drop = FALSE]-line[-nrow(line), , drop = FALSE])^2)))
      data.frame(boundary = leaders$boundary[i], vertex = seq_len(nrow(line)),
        x = line[, 1], y = line[, 2], coordinate_space = "composition", units = "mm")
    })
    leader_paths <- do.call(rbind, paths)
    grobs <- lapply(paths, function(path) grid::polylineGrob(path$x, path$y,
      default.units = "mm", gp = grid::gpar(col = focus$colour,
        lwd = focus$linewidth * 72.27/25.4, lty = focus$linetype)))
  }
  highlight <- function(front) {
    offset <- focus$width/2
    grid::polygonGrob(
      c(front$x+offset*front$normal_x, rev(front$x-offset*front$normal_x)),
      c(front$y+offset*front$normal_y, rev(front$y-offset*front$normal_y)),
      default.units = "mm", gp = grid::gpar(col = NA,fill = scales::alpha(focus$colour,focus$alpha)))
  }
  highlights <- stats::setNames(list(highlight(qf),highlight(sf)), c(q$panel, s$panel))
  length_mm <- function(front) sum(sqrt(diff(front$x)^2+diff(front$y)^2))
  source_mm <- length_mm(qf); detail_mm <- length_mm(sf)
  if (any(!is.finite(c(source_mm,detail_mm))) || min(source_mm,detail_mm)<=0)
    ggchord_stop("magnify_ggchord(): selected interval has no measurable span at the requested scale")
  report <- region[c("accver", "start", "end")]
  report$source_panel <- q$panel; report$detail_panel <- s$panel
  report$source_instance <- q$instance; report$detail_instance <- s$instance
  report$source_origin <- q$origin; report$detail_origin <- s$origin
  report$sequence_unit <- q$unit; report$source_length_mm <- source_mm; report$detail_length_mm <- detail_mm
  report$length_ratio <- detail_mm/source_mm; report$stat_scope <- "explicit_detail_plot"
  report$collision_status <- if (nrow(leaders)) paste(unique(leaders$route_status), collapse = ";") else "no_leaders"
  report$highlight_width_mm <- focus$width; report$highlight_alpha <- focus$alpha
  report$colour <- focus$colour; report$connector <- focus$connector
  report$leader_linewidth <- focus$linewidth; report$leader_linetype <- focus$linetype
  qf$panel_id <- q$panel; qf$instance_id <- q$instance
  sf$panel_id <- s$panel; sf$instance_id <- s$instance
  list(input = focus$input, underlay = grobs, overlay = highlights, regions = report, fronts = rbind(qf,sf), leaders = leaders, leader_paths = leader_paths)
}
