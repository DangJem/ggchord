# The source computation happens before reference-window projection. The
# drawing-side identity Stat makes the source results available to after_stat.
StatGgchordLink <- ggplot2::ggproto(
  "StatGgchordLink", ggplot2::Stat, required_aes = c("x", "y"),
  compute_panel = function(data, scales) {
    keys <- grep("^.ggchord_link_stat_", names(data), value = TRUE)
    for (key in keys) data[[sub("^.ggchord_link_stat_", "", key)]] <- data[[key]]
    data
  }
)

ggchord_ribbon_stat <- function(name, dots) {
  if (!is.character(name) || !name %in% c("link_bundle", "link_density"))
    return(list(geom_parameters = dots))
  controls <- if (name == "link_bundle") c("bins", "min_bundle") else c("bins", "normalization")
  arguments <- dots[intersect(names(dots), controls)]
  list(stat = name, arguments = arguments,
    geom_parameters = dots[setdiff(names(dots), names(arguments))])
}

ggchord_apply_ribbon_stat <- function(layer, specification) {
  type <- specification$stat
  fields <- if (type == "link_bundle") c("n_members", "sum_weight", "is_summary") else
    c("bin_weight", "bin_area", "density_raw", "density")
  layer$stat <- StatGgchordLink
  layer$ggchord_params$link_stat <- list(type = type, arguments = specification$arguments)
  for (field in fields)
    layer$mapping[[paste0(".ggchord_link_stat_", field)]] <- as.name(field)
  layer$ggchord_role_aes <- c(layer$ggchord_role_aes, "group", "weight")
  layer
}

ggchord_prepare_link_statistics <- function(plot, references) {
  for (layer in plot$layers) {
    p <- layer$ggchord_params
    specification <- p$link_stat
    if (is.null(specification)) next
    caller <- ggchord_layer_input_context(layer)
    input <- ggchord_resolve_layer_input(layer, ggchord_layer_fallback(layer, plot$ggchord$data))
    if (is.null(input)) ggchord_stop(caller, ": supply interval link data")
    members <- ggchord_facet_stat_input(input, plot, layer)
    compute <- function(raw) {
      data <- as_ggchord_link(raw)
      group <- raw[intersect("group", names(raw))]
      arguments <- specification$arguments
      grid <- ggchord_link_bins(data, references, arguments$bins %||% 80,
        if (is.null(raw$weight)) NULL else as.numeric(raw$weight), group, caller)
      if (identical(specification$type, "link_bundle"))
        return(ggchord_bundle_binned_links(data, grid, arguments$min_bundle %||% 2))
      fields <- c("bin_weight", "bin_area", "density_raw", "density")
      if (any(fields %in% names(data))) ggchord_stop(caller, ": density fields collide with source columns")
      normalization <- match.arg(arguments$normalization %||% "none", c("none", "pair_max"))
      data$bin_weight <- if (nrow(grid)) stats::ave(grid$weight, grid$bin_group, FUN = sum) else numeric()
      data$bin_area <- grid$area
      data$density_raw <- data$bin_weight / grid$area
      data$density <- data$density_raw
      if (normalization == "pair_max" && nrow(data)) {
        pair_max <- stats::ave(data$density_raw, grid$pair_group, FUN = max)
        data$density <- ifelse(pair_max > 0, data$density_raw / pair_max, 0)
      }
      attr(data, "ggchord_density_normalization") <- normalization
      data$source_members <- I(ggchord_link_source_members(data))
      data
    }
    if (is.null(members)) result <- compute(input) else {
      panels <- plot$ggchord$source_facet_setup$panels
      keys <- ggchord_facet_vars(plot)
      result <- ggchord_rbind_fill(lapply(seq_len(nrow(panels)), function(panel) {
        rows <- members$.ggchord_input_index[as.integer(members$PANEL) == panel]
        part <- input[rows, , drop = FALSE]
        if (!".source_row" %in% names(part)) part$.source_row <- rows
        out <- compute(part)
        out$.ggchord_stat_panel <- rep.int(panel, nrow(out))
        for (key in keys) out[[key]] <- rep(panels[[key]][panel], nrow(out))
        out
      }))
    }
    capture <- p$stage_capture
    if (!is.null(capture)) {
      if ("mapped" %in% capture$requested) capture$mapped <- input
      if ("statistic" %in% capture$requested) capture$statistic <- result
    }
    mapping <- layer$ggchord_input_mapping
    for (name in intersect(layer$ggchord_role_aes, names(mapping)))
      mapping[[name]] <- rlang::new_quosure(rlang::sym(name), emptyenv())
    layer$ggchord_input_mapping <- mapping
    layer$ggchord_input_data <- result
    layer$ggchord_params$link_stat_source <- input
    if (!is.null(p$visual_source)) {
      p$visual_source$data <- result
      p$visual_source$rows <- seq_len(nrow(result))
      p$visual_source$values <- new.env(parent = emptyenv())
    }
  }
  plot
}

#' Aggregate source intervals in a complete reference grid
#' @export
#' @param mapping Aesthetic mapping from `ggplot2::aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param geom Compatible geometry to consume this statistic.
#' @param position A compatible Position for placement of the layer geometry.
#' @param bins Number of bins for both axes, or a named vector with q and s counts.
#' @param min_bundle Minimum number of source records required to aggregate a bin.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Source references are `qaccver` and `saccver`. Point links consume `qpos` and `spos`. Interval links and alignment segments consume closed interval endpoints on both references. Supply reliable `qstrand` and `sstrand` or explicit relative direction when the geometry requires it.
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
#' ggchord(data = seqs) + geom_seq() +
#'   stat_link_bundle(data = links, bins = c(q = 20, s = 20),
#'     min_bundle = 2, aes(fill = after_stat(n_members))) +
#'   scale_fill_viridis_c(name = "Members") + coord_chord()
stat_link_bundle <- function(mapping = NULL, data = NULL, geom = "link_ribbon",
    position = "identity", bins = 80, min_bundle = 2, ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!identical(geom, "link_ribbon")) ggchord_stop("stat_link_bundle(): requires an interval link geometry")
  geom_link_ribbon(mapping, data, stat = "link_bundle", position = position,
    ..., bins = bins, min_bundle = min_bundle, na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Compute source interval density
#' @export
#' @param mapping Aesthetic mapping from `ggplot2::aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param geom Compatible geometry to consume this statistic.
#' @param position A compatible Position for placement of the layer geometry.
#' @param bins Number of bins for both axes, or a named vector with q and s counts.
#' @param normalization Return raw bin density or also normalize relative to the pair maximum.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Source references are `qaccver` and `saccver`. Point links consume `qpos` and `spos`. Interval links and alignment segments consume closed interval endpoints on both references. Supply reliable `qstrand` and `sstrand` or explicit relative direction when the geometry requires it.
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
#' ggchord(data = seqs) + geom_seq() +
#'   stat_link_density(data = links,
#'     aes(fill = after_stat(density)), bins = 20,
#'     normalization = "pair_max") +
#'   scale_fill_viridis_c(name = "Relative density", limits = c(0, 1)) +
#'   coord_chord()
stat_link_density <- function(mapping = NULL, data = NULL, geom = "link_ribbon",
    position = "identity", bins = 80, normalization = "none", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!identical(geom, "link_ribbon")) ggchord_stop("stat_link_density(): requires an interval link geometry")
  geom_link_ribbon(mapping, data, stat = "link_density", position = position,
    ..., bins = bins, normalization = normalization, na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}
