# Focus+context uses the same panel/frame projection as biological links, but
# its leaders and highlights never enter the relationship/HSP pipeline.

#' Compose a source region and an explicit local detail plot
#'
#' Marks one real source interval in a Circular/Chord plot and connects it to
#' the same interval in an independently constructed Collinear detail plot.
#' Returns a Fusion, not a new Coord or a crop of an existing image. Ordinary
#' ggsave rebuilds both panels, highlights and leaders at the actual target size.
#'
#' Build `detail_plot` explicitly with [focus_ggchord_data()] or other documented
#' source-aware operations. This function does not copy/filter arbitrary layers
#' or silently choose a local/global Stat policy. The supplied detail plot owns
#' its data, mappings, statistics and labels. Original source origins, including
#' renamed/repeated focus windows, are used to resolve the same interval.
#'
#' @param plot A ggchord Circular or Chord source plot.
#' @param region One-row data frame with `accver`, `start`, `end` in original
#'   one-based inclusive source coordinates. Bounds must be ordered positive
#'   integers. Optional `instance` and `detail_instance` select exact display
#'   instances in the source and detail panels. Ambiguity fails; origin-wrapping
#'   intervals must be split before using this single-window interface.
#' @param detail_plot An explicit horizontal or vertical Collinear ggchord plot.
#'   The selected source interval must be fully visible in both panels, in the
#'   same bp/aa unit. Ordinary Cartesian character plots are not supported here.
#' @param side One of right, top, left or bottom for the detail panel.
#' @param center_size,detail_size Optional panel width/height fractions, following
#'   [ggchord_fusion()] and [fusion_panel()]. Defaults reserve a smaller source
#'   panel and a wider detail container for the selected side. Overlap fails.
#' @param connector `"leader"` draws two decorative boundary leaders; `"none"`
#'   retains the highlights only. These are not biological ribbons or alignments.
#' @param colour Fixed highlight/leader colour.
#' @param alpha Highlight opacity in `[0,1]`. Leaders use the opaque colour.
#' @param highlight_width Positive highlight band width in physical mm.
#' @param linewidth Nonnegative leader width using the usual ggplot2 convention.
#' @param linetype Leader line type, default dashed.
#' @return A [ggchord_fusion()] object. [export_ggchord_fusion()] additionally
#'   returns focus_regions, focus_fronts and leaders in mm, plus original
#'   focus_region_data (without metadata-name collisions), preserving the input
#'   region, exact instances, source origins and sampled physical path-length ratio.
#'   Stat scope is recorded as explicit_detail_plot. Leaders do not appear in
#'   connections or relationship_data. Global label/leader collision avoidance,
#'   arbitrary inset boxes, nested/multiple insets and automatic layer rebuilding
#'   are not implemented. The measured ratio may be below one if the chosen
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
magnify_ggchord <- function(plot, region, detail_plot,
    side = c("right", "top", "left", "bottom"), center_size = NULL, detail_size = NULL,
    connector = c("leader", "none"), colour = "#B45309", alpha = .2,
    highlight_width = 1.5, linewidth = .35, linetype = 2) {
  ggchord_fusion_check_plot(plot, center = TRUE)
  ggchord_fusion_check_plot(detail_plot)
  caller <- "magnify_ggchord()"
  ggchord_require_columns(region, c("accver", "start", "end"), caller)
  if (nrow(region) != 1L || !is.character(region$accver) || is.na(region$accver) || !nzchar(region$accver) ||
      !is.numeric(region$start) || !is.numeric(region$end) ||
      any(!is.finite(c(region$start,region$end))) || region$start < 1 || region$end < region$start ||
      any(c(region$start,region$end) != floor(c(region$start,region$end))) || region$end > 2^53-1)
    ggchord_stop(caller, ": provide one ordered positive integral source interval; split origin-wrapping regions explicitly")
  for (name in intersect(c("instance", "detail_instance"), names(region)))
    if (!is.character(region[[name]]) || is.na(region[[name]]) || !nzchar(region[[name]]))
      ggchord_stop(caller, ": instance selectors must be nonempty strings")
  if (!is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<0 || alpha>1 ||
      !is.numeric(highlight_width) || length(highlight_width)!=1L || !is.finite(highlight_width) || highlight_width<=0 ||
      !is.numeric(linewidth) || length(linewidth)!=1L || !is.finite(linewidth) || linewidth<0)
    ggchord_stop(caller, ": alpha must be in [0,1], highlight_width positive and linewidth nonnegative")
  if (!is.character(colour) || length(colour)!=1L || is.na(colour))
    ggchord_stop(caller, ": colour must be a single colour string")
  tryCatch(grDevices::col2rgb(colour),error=function(e) ggchord_stop(caller, ": invalid colour"))
  if (length(linetype)!=1L || is.na(linetype)) ggchord_stop(caller, ": linetype must be scalar")
  side <- match.arg(side)
  lateral <- side %in% c("right", "left")
  if (is.null(center_size)) center_size <- if (lateral) c(.4,.8) else c(.8,.4)
  if (is.null(detail_size)) detail_size <- if (lateral) c(.42,.88) else c(.88,.42)
  result <- ggchord_fusion(plot, center_size = center_size) +
    fusion_panel(detail_plot, "detail", side = side, size = detail_size)
  result$focus <- list(input = region, region = as.data.frame(region), connector = match.arg(connector),
    colour = colour, alpha = alpha, width = highlight_width, linewidth = linewidth, linetype = linetype)
  ggplot2::set_last_plot(result)
  result
}

# Build only after the real panel sizes, position scales and Coord are final.
ggchord_build_focus <- function(focus, panels) {
  empty <- list(input = data.frame(), underlay = list(), overlay = list(), regions = data.frame(), fronts = data.frame(), leaders = data.frame())
  if (is.null(focus)) return(empty)
  region <- focus$region
  resolve <- function(panel, instance) ggchord_resolve_panel_interval(panels,
    region$accver, region$start, region$end, panel, instance,
    paste0("magnify_ggchord() ",panel,": "))
  q <- resolve("center", region$instance)
  s <- resolve("detail", region$detail_instance)
  if (q$unit != s$unit) ggchord_stop("magnify_ggchord(): source and detail units must agree")
  qf <- ggchord_fusion_front(panels[[q$panel]], q)
  sf <- ggchord_fusion_front(panels[[s$panel]], s)
  grobs <- list()
  ends <- c(1L,nrow(qf))
  leaders <- data.frame()
  if (focus$connector == "leader") {
    leaders <- data.frame(boundary = c("start", "end"),
      source_position_boundary = qf$position_boundary[ends],
      x = qf$x[ends],y = qf$y[ends],xend = sf$x[ends],yend = sf$y[ends],
      coordinate_space = "composition",units = "mm",semantic = "decorative_focus_leader")
    grobs[[1L]] <- grid::segmentsGrob(leaders$x, leaders$y, leaders$xend, leaders$yend,
      default.units = "mm",gp = grid::gpar(col = focus$colour,
        lwd = focus$linewidth * 72.27/25.4,lty = focus$linetype))
  }
  highlight <- function(front) {
    offset <- focus$width/2
    grid::polygonGrob(
      c(front$x+offset*front$normal_x, rev(front$x-offset*front$normal_x)),
      c(front$y+offset*front$normal_y, rev(front$y-offset*front$normal_y)),
      default.units = "mm", gp = grid::gpar(col = NA,fill = scales::alpha(focus$colour,focus$alpha)))
  }
  highlights <- list(highlight(qf),highlight(sf))
  length_mm <- function(front) sum(sqrt(diff(front$x)^2+diff(front$y)^2))
  source_mm <- length_mm(qf); detail_mm <- length_mm(sf)
  if (any(!is.finite(c(source_mm,detail_mm))) || min(source_mm,detail_mm)<=0)
    ggchord_stop("magnify_ggchord(): selected interval has no measurable span at the requested scale")
  report <- region[c("accver", "start", "end")]
  report$source_instance <- q$instance; report$detail_instance <- s$instance
  report$source_origin <- q$origin; report$detail_origin <- s$origin
  report$sequence_unit <- q$unit; report$source_length_mm <- source_mm; report$detail_length_mm <- detail_mm
  report$length_ratio <- detail_mm/source_mm; report$stat_scope <- "explicit_detail_plot"
  report$collision_status <- "not globally measured"
  report$highlight_width_mm <- focus$width; report$highlight_alpha <- focus$alpha
  report$colour <- focus$colour; report$connector <- focus$connector
  report$leader_linewidth <- focus$linewidth; report$leader_linetype <- focus$linetype
  qf$panel_id <- "center"; qf$instance_id <- q$instance
  sf$panel_id <- "detail"; sf$instance_id <- s$instance
  list(input = focus$input, underlay = grobs, overlay = highlights, regions = report, fronts = rbind(qf,sf),leaders = leaders)
}
