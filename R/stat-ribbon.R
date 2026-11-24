# The source computation happens before reference-window projection. The
# drawing-side identity Stat makes the source results available to after_stat.
StatGgchordLink <- ggproto(
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
    c("q_bin", "s_bin", "bin_weight", "bin_area", "density_raw", "density",
      "normalization", "density_unit")
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
    fields <- c("qstart", "qend", "sstart", "send")
    ggchord_require_columns(input, fields, caller)
    fields <- c(fields, intersect("weight", names(input)))
    for (field in fields) {
      value <- input[[field]]
      if (is.logical(value) && all(is.na(value))) input[[field]] <- as.numeric(value)
    }
    if (!is.null(input$weight) && (!is.numeric(input$weight) ||
        any(is.infinite(input$weight) | is.nan(input$weight) | input$weight < 0, na.rm = TRUE)))
      ggchord_stop(caller, ": weight must be numeric, finite and nonnegative")
    missing <- !complete.cases(input[fields])
    if (!".source_row" %in% names(input)) input$.source_row <- seq_len(nrow(input))
    if (any(missing) && !isTRUE(layer$geom_params$na.rm))
      warning(caller, ": omitted ", sum(missing), " incomplete relationship row(s)", call. = FALSE)
    observations <- input[!missing, , drop = FALSE]
    members <- ggchord_facet_stat_input(observations, plot, layer)
    arguments <- specification$arguments
    density <- identical(specification$type, "link_density")
    normalization <- if (density)
      ggchord_match_arg(arguments$normalization %||% "none", c("none", "pair_max"))
    compute <- function(raw) {
      data <- as_ggchord_link(raw)
      group <- raw[intersect("group", names(raw))]
      grid <- ggchord_link_bins(data, references, arguments$bins %||% 80,
        raw$weight, group, caller)
      if (!density)
        return(ggchord_bundle_binned_links(data, grid, arguments$min_bundle %||% 2))
      computed <- c("q_bin", "s_bin", "bin_weight", "bin_area", "density_raw", "density",
        "normalization", "density_unit")
      if (any(computed %in% names(data))) ggchord_stop(caller, ": density fields collide with source columns")
      data$q_bin <- grid$qbin
      data$s_bin <- grid$sbin
      data$bin_weight <- if (nrow(grid)) ave(grid$weight, grid$bin_group, FUN = sum) else numeric()
      data$bin_area <- grid$area
      data$density_raw <- data$bin_weight / grid$area
      data$density <- data$density_raw
      if (normalization == "pair_max" && nrow(data)) {
        pair_max <- ave(data$density_raw, grid$pair_group, FUN = max)
        data$density <- ifelse(pair_max > 0, data$density_raw / pair_max, 0)
      }
      data$normalization <- rep(normalization, nrow(data))
      data$density_unit <- if (nrow(grid)) paste0("weight/(", grid$q_unit, " × ", grid$s_unit, ")") else character()
      data$source_members <- I(ggchord_link_source_members(data))
      data
    }
    cache <- p$link_stat_cache
    if (!is.null(cache$result)) result <- cache$result else if (is.null(members)) result <- compute(observations) else {
      panels <- plot$ggchord$source_facet_setup$panels
      keys <- ggchord_facet_vars(plot)
      result <- ggchord_rbind_fill(lapply(seq_len(nrow(panels)), function(panel) {
        rows <- members$.ggchord_input_index[as.integer(members$PANEL) == panel]
        part <- observations[rows, , drop = FALSE]
        out <- compute(part)
        out$.ggchord_stat_panel <- rep.int(panel, nrow(out))
        for (key in keys) out[[key]] <- rep(panels[[key]][panel], nrow(out))
        out
      }))
    }
    if (density) attr(result, "ggchord_density_normalization") <- normalization
    if (!is.null(cache)) cache$result <- result
    capture <- p$stage_capture
    if (!is.null(capture)) {
      if ("mapped" %in% capture$requested) capture$mapped <- input
      if ("statistic" %in% capture$requested) capture$statistic <- result
    }
    mapping <- layer$ggchord_input_mapping
    for (name in intersect(layer$ggchord_role_aes, names(mapping)))
      mapping[[name]] <- new_quosure(sym(name), emptyenv())
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
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param geom Compatible geometry to consume this statistic.
#' @param position A compatible Position for placement of the layer geometry.
#' @param bins Number of bins for both axes, or a named vector with q and s counts.
#' @param min_bundle Minimum number of source records required to aggregate a bin.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Omit records with missing endpoint coordinates or weights.
#'   Set to `TRUE` to omit them silently. Invalid references, units and
#'   out-of-bounds coordinates remain errors.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer whose statistic output retains the original members
#'   in `source_members`, with `n_members` and `sum_weight` for each result.
#' @details The statistic bins the midpoints of closed query and subject
#'   intervals over their complete references. Bins are evaluated within
#'   reference pairs, orientation, explicit groups and facets. Bins below
#'   `min_bundle` retain separate source records.
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
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!identical(geom, "link_ribbon")) ggchord_stop("stat_link_bundle(): requires an interval link geometry")
  geom_link_ribbon(mapping, data, stat = "link_bundle", position = position,
    ..., bins = bins, min_bundle = min_bundle, na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Compute source interval density
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param geom Compatible geometry to consume this statistic.
#' @param position A compatible Position for placement of the layer geometry.
#' @param bins Number of bins for both axes, or a named vector with q and s counts.
#' @param normalization Return raw bin density or also normalize relative to the pair maximum.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Omit records with missing endpoint coordinates or weights.
#'   Set to `TRUE` to omit them silently. Invalid references, units and
#'   out-of-bounds coordinates remain errors.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer. Its statistic output retains the source rows and
#'   adds `q_bin`, `s_bin`, `bin_weight`, `bin_area`, `density_raw`,
#'   `density`, `normalization`, `density_unit` and `source_members`.
#' @details Each interval midpoint is assigned to a bin along its full
#'   reference. The query and subject bin widths determine `bin_area`.
#'   The sum of nonnegative weights in a cell gives `bin_weight`, with
#'   unit weight per row by default. `density_raw` is `bin_weight / bin_area`,
#'   expressed in weight per product of the two source units.
#'
#'   With `normalization = "pair_max"`, `density` is relative to the maximum
#'   in the same reference pair, orientation, group and facet. An all-zero
#'   group returns zero. The normalization is recorded in the statistic
#'   output, including when the plot uses facets. Use
#'   `aes(fill = after_stat(density))` to map the result to colour.
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
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!identical(geom, "link_ribbon")) ggchord_stop("stat_link_density(): requires an interval link geometry")
  geom_link_ribbon(mapping, data, stat = "link_density", position = position,
    ..., bins = bins, normalization = normalization, na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}
