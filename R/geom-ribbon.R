#' Draw interval-to-interval links
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param route Choose curved or straight routes, or a shared corridor where supported.
#'   On a curved collinear route, an explicit positive `bulge` clears the high
#'   coordinate end of intervening rows by that fraction of the longest visible
#'   backbone. Radial curves with positive `bulge` clear intermediate spoke tips.
#' @param fill_mode Use solid fill or map independent query and subject endpoint colours.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#'   The default outline is muted green (#53776A), 0.18 mm wide. Set
#'   `colour = NA` to omit it. When `linetype` is omitted and unmapped, visible outlines connecting
#'   nonadjacent rows or radial spokes are dashed. Radial adjacency is cyclic.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @param gap Nonnegative clearance in mm between a ribbon front and the edge of
#'   the sequence or avoided feature. Accepts flexible per-sequence values.
#' @param avoid Use `"none"` for a constant sequence gap (the default), or `"smooth"`
#'   to follow the facing feature outlines at ribbon endpoints.
#' @param avoid_strength A number from 0 to 1 controlling how far the rounded
#'   transitions extend. With `avoid = "smooth"`, 0 follows contours most closely
#'   and 1 gives wider transitions. The minimum clearance applies at every strength.
#' @return A ggplot2 layer for addition to a plot.
#' @details Source references are `qaccver` and `saccver`. Point links consume `qpos` and `spos`. Interval links and alignment segments consume closed interval endpoints on both references. Supply reliable `qstrand` and `sstrand` or explicit relative direction when the geometry requires it.
#' Ribbons retain same-reference and non-neighbouring relationships. A
#' same-row linear ribbon needs a curved route. A shared route retains
#' independent source relationships and uses `corridor_width` in mm.
#' Smooth avoidance follows the facing feature outlines while preserving the
#' original interval endpoints and the route through intervening references.
#' Tip decoration is recorded separately from the source interval fronts.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_link_ribbon(aes(fill = pident), data = links, avoid = "smooth",
#'     gap = 0.4, avoid_strength = 0.35, id = "matches") +
#'   coord_chord()
geom_link_ribbon <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", route = "curve", fill_mode = "solid", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = NULL,
    gap = 0.4, avoid = "none", avoid_strength = 0.5) {
  call <- match.call(function(...) NULL, sys.call(), envir = parent.frame())
  if ("fill" %in% names(call) && !"fill_mode" %in% names(call)) {
    # Reserve the exact visual argument so R cannot treat it as an
    # abbreviation of fill_mode before ... is evaluated.
    signature <- sys.function()
    formals(signature) <- c(as.list(formals(signature)), alist(fill = ))
    call <- match.call(signature, call, envir = parent.frame())
    if (is.null(call[["fill_mode"]])) call[["fill_mode"]] <- "solid"
    return(eval(call, parent.frame()))
  }
  ggchord_check_argument_names(sys.call(), sys.function(), dots = "fill")
  caller <- "geom_link_ribbon()"
  fill_mode <- match.arg(fill_mode, c("solid", "endpoint"))
  position <- ggchord_as_feature_position(position, caller)
  if (ggchord_local_position(position) && !identical(position$placement, "offset"))
    ggchord_stop(caller, ": use identity or position_offset()")
  dots <- ggchord_colour_dots(list(...), caller)
  gap <- ggchord_physical_mm(gap, "geom_link_ribbon(): gap", TRUE)
  avoid <- match.arg(avoid, c("none", "smooth"))
  if (!is.numeric(avoid_strength) || length(avoid_strength) != 1L ||
      !is.finite(avoid_strength) || avoid_strength < 0 || avoid_strength > 1)
    ggchord_stop(caller, ": avoid_strength must be a number from 0 to 1")
  controls <- ggchord_link_route(route, dots, caller)
  dots <- controls$dots
  controls$dots <- NULL
  tip <- match.arg(dots$tip %||% "none", c("none", "from", "to", "both"))
  if (tip == "none" && "tip_length" %in% names(dots))
    ggchord_stop(caller, ": tip_length requires a tip")
  tip_length <- ggchord_physical_mm(dots$tip_length %||% 2, "tip_length", TRUE)
  if (length(tip_length) != 1L || tip_length <= 0)
    ggchord_stop(caller, ": tip_length must be positive")
  dots[c("tip", "tip_length")] <- NULL
  used <- c(names(mapping), names(dots))
  if ((fill_mode == "solid" && any(c("fill_from", "fill_to") %in% used)) ||
      (fill_mode == "endpoint" && "fill" %in% used))
    ggchord_stop(caller, ": select the fill channels for the requested fill_mode")
  special <- ggchord_ribbon_stat(stat, dots)
  dots <- special$geom_parameters
  geom <- GeomRibbonFlow
  ggchord_check_geom_params(dots, geom, caller,
    stat = if (is.null(special$stat)) stat else ggplot2::StatIdentity)
  parameters <- if (is.null(special$stat)) ggchord_split_source_stat(stat, dots, geom) else list(geom = dots)
  layer <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), group = integer()),
    mapping = aes(x = x, y = y, group = group),
    stat = "identity", geom = geom, position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes,
    check.aes = FALSE, params = c(list(na.rm = na.rm, flow = fill_mode == "endpoint",
      tip = tip, tip_length = tip_length), parameters$geom)
  )
  layer$ggchord_type <- "ribbon"
  layer$ggchord_params <- list(id = id, type = "ribbon", route = controls,
    ribbon_ctrl_point = controls, front_gap_mm = gap,
    avoid_strength = avoid_strength, link_avoid = avoid,
    feature_position = position, tip = tip, fill_mode = fill_mode,
    physical_layout = TRUE,
    ribbon_color_scheme = "single", ribbon_alpha = .6,
    link_branch = if (controls$route == "shared") "query" else "none",
    link_branch_fraction = .15, source_stat = parameters$stat)
  layer <- ggchord_capture_layer_input(layer, data, mapping,
    c("qaccver", "qstart", "qend", "saccver", "sstart", "send",
      "qstrand", "sstrand", "direction"))
  if (!is.null(special$stat)) layer <- ggchord_apply_ribbon_stat(layer, special)
  layer
}

# Rendering hooks used during layer assembly. Ribbons use the standard
# aesthetic dictionary in both endpoint and solid modes.
ggchord_inherit_ribbon_linetype <- function(layer, mapping) layer
make_ribbon_geom <- function(outline = FALSE, linetype = FALSE) GeomRibbonFlow
