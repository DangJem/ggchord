#' Set a target-size rendering policy
#'
#' Scope size adaptation to default appearance and guide layout. Explicit
#' geometric dimensions and source-coordinate computations remain locked.
#'
#' @param mode Use `"fixed"` (default) to retain physical styles, `"scale"` to
#'   resize unlocked default styles, or `"adaptive"` to bound the text
#'   adjustment and reflow unconstrained guides.
#' @param reference_size Width and height for the scaling baseline. Required
#'   in scale mode. Adaptive mode uses 8 by 6 inches when omitted.
#' @param units Physical units of reference_size.
#' @param text_scale_range Lower and upper text factors for adaptive mode.
#' @param min_text_size Optional floor in mm for automatic reductions of
#'   theme text. Defaults already below that size are kept at their size.
#' @return A policy that can be added to a plot or Fusion.
#' @details The target is the actual graphics device or the region allocated to
#' a child plot in Fusion. Scale mode uses the smaller of the width and height
#' ratios. Adaptive mode uses the fourth root of the target-to-reference area
#' ratio, limited by `text_scale_range`. Unlocked line widths and spacing use
#' the square root of the adaptive text factor. Explicit or mapped text sizes,
#' feature thickness and offsets keep their declared values. Explicit guide
#' settings take priority. Add the policy to an individual child plot to adapt
#' that child, or to Fusion to adapt its own outer decorations and guides.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' p_basic + render_policy(mode = "adaptive", text_scale_range = c(.9, 1.1))
render_policy <- function(mode = "fixed", reference_size = NULL, units = "in",
    text_scale_range = c(0.8, 1.2), min_text_size = NULL) {
  ggchord_local_console_errors()
  caller <- "render_policy()"
  mode <- ggchord_match_arg(mode, c("fixed", "scale", "adaptive"))
  units <- ggchord_match_arg(units, c("in", "mm", "cm"))
  if (mode == "scale" && is.null(reference_size))
    ggchord_stop(caller, ": scale mode requires reference_size")
  if (!is.null(reference_size) && (!is.numeric(reference_size) ||
      length(reference_size) != 2L || any(!is.finite(reference_size) | reference_size <= 0)))
    ggchord_stop(caller, ": reference_size must contain positive width and height")
  if (!is.numeric(text_scale_range) || length(text_scale_range) != 2L ||
      any(!is.finite(text_scale_range) | text_scale_range <= 0) ||
      text_scale_range[1L] > text_scale_range[2L])
    ggchord_stop(caller, ": text_scale_range must contain positive increasing multipliers")
  if (mode != "adaptive" && !missing(text_scale_range))
    ggchord_stop(caller, ": text_scale_range only applies in adaptive mode")
  if (!is.null(min_text_size) && (!is.numeric(min_text_size) ||
      length(min_text_size) != 1L || !is.finite(min_text_size) || min_text_size <= 0))
    ggchord_stop(caller, ": min_text_size must be positive physical mm")
  if (mode == "fixed" && !is.null(min_text_size))
    ggchord_stop(caller, ": min_text_size applies to adaptive or scale mode")
  size_in <- unname(reference_size) / switch(units, "in" = 1, "cm" = 2.54, "mm" = 25.4)
  structure(list(mode = mode, reference_size = if (length(size_in)) size_in,
    text_scale_range = unname(text_scale_range), min_text_size = min_text_size),
    class = "ggchord_layout_policy")
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.ggchord_layout_policy <- function(object, plot, ...) {
  ggchord_local_console_errors()
  if (!inherits(plot, "ggchord")) ggchord_stop("render_policy(): add the policy to ggchord or Fusion")
  plot$ggchord$layout_policy <- unclass(object)
  plot$ggchord$ref <- new.env(parent = emptyenv())
  plot$ggchord$layout <- NULL
  plot
}

# Build-local policy: data and the source plot stay untouched.
ggchord_adaptive_factors <- function(policy, canvas = ggchord_render_dimensions()) {
  measured <- length(canvas) == 2L && all(is.finite(canvas)) && all(canvas > 0)
  reference <- policy$reference_size %||% c(8, 6)
  scale <- if (measured) canvas / 25.4 / reference else c(1, 1)
  factor <- switch(policy$mode, fixed = 1, scale = min(scale),
    adaptive = if (!measured) 1 else max(policy$text_scale_range[1L],
      min(policy$text_scale_range[2L], prod(scale)^.25)))
  spacing <- if (policy$mode == "scale") factor else sqrt(factor)
  list(mode = policy$mode, canvas_mm = canvas, reference_size_in = policy$reference_size,
    text_factor = factor, line_factor = spacing, spacing_factor = spacing,
    min_text_size_mm = policy$min_text_size,
    measurement = if (measured) "device_dimensions" else "unavailable")
}

ggchord_apply_adaptive_styles <- function(plot) {
  policy <- plot$ggchord$layout_policy
  if (is.null(policy)) return(plot)
  plot <- ggchord_render_plot_copy(plot)
  report <- ggchord_adaptive_factors(policy)
  explicit <- plot$ggchord$theme_explicit %||% list()
  tree <- ggplot2::get_element_tree()
  ancestors <- function(name) unique(c(name, unlist(lapply(tree[[name]]$inherit, ancestors))))
  supplied <- function(name, field) {
    any(vapply(ancestors(name), function(key) {
      value <- explicit[[key]]
      !is.null(value) && !inherits(value, "element_blank") && !is.null(value[[field]])
    }, logical(1)))
  }
  theme <- plot$theme
  sizes <- list()
  for (name in names(theme)) {
    raw <- theme[[name]]
    if (!inherits(raw, c("element_text", "element_line"))) next
    field <- if (inherits(raw, "element_text")) "size" else "linewidth"
    el <- ggchord_theme_element(plot, name)
    value <- el[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) next
    locked <- supplied(name, field)
    target <- value
    if (!locked && policy$mode != "fixed") {
      target <- value * if (field == "size") report$text_factor else report$line_factor
      if (field == "size" && !is.null(policy$min_text_size))
        target <- max(target, min(value, policy$min_text_size * ggplot2::.pt))
    }
    raw[[field]] <- target
    theme[[name]] <- raw
    sizes[[length(sizes) + 1L]] <- data.frame(element = name, property = field,
      requested = value, effective = target, user_locked = locked,
      units = if (field == "size") "pt" else "mm")
  }
  for (name in c("plot.margin", "legend.spacing", "legend.spacing.x", "legend.spacing.y", "legend.box.spacing"))
    if (is.unit(theme[[name]]) && is.null(explicit[[name]]))
      theme[[name]] <- theme[[name]] * report$spacing_factor
  geom_styles <- list()
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    defaults <- layer$ggchord_params$style_defaults
    fields <- intersect(grep('(^|_)linewidth$', names(defaults), value = TRUE), names(layer$aes_params))
    user <- layer$ggchord_input_mapping
    if (isTRUE(layer$inherit.aes)) {
      inherited <- ggchord_normalize_mapping(plot$mapping)
      user <- c(user, inherited[setdiff(names(inherited), names(user))])
    }
    copied <- FALSE
    for (field in fields) {
      value <- layer$aes_params[[field]]
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) next
      mapped <- isTRUE(defaults[[field]]) && field %in% names(user)
      locked <- !isTRUE(defaults[[field]]) || mapped
      target <- if (!locked && policy$mode != 'fixed') value*report$line_factor else value
      if (!identical(target, value)) {
        if (!copied) layer <- reconstruct_layer(layer, layer$data)
        copied <- TRUE
        layer$aes_params[[field]] <- target
      }
      geom_styles[[length(geom_styles)+1L]] <- data.frame(
        layer_id = layer$ggchord_layer_id %||% paste0('layer-', i), aesthetic = field,
        default_or_fixed = value, effective_fixed = if (mapped) NA_real_ else target,
        mapped = mapped, user_locked = locked, units = 'ggplot2 linewidth aesthetic')
    }
    if (copied) plot$layers[[i]] <- layer
  }
  plot$theme <- theme
  report$styles <- ggchord_rbind_fill(sizes) %||% data.frame()
  report$geom_styles <- ggchord_rbind_fill(geom_styles) %||% data.frame()
  plot$ggchord$adaptive_report <- report
  plot
}

# Field-wise merge preserves an explicitly supplied size even when a later
# theme changes only colour. A new complete theme resets the old constraints.
ggchord_merge_theme_constraints <- function(old, incoming, complete = FALSE) {
  if (complete) old <- list()
  for (name in names(incoming)) {
    value <- incoming[[name]]
    if (is.null(value)) next
    if (inherits(value, c("element_text", "element_line", "element_rect")) &&
        inherits(old[[name]], class(value)[1])) {
      value <- ggplot2::merge_element(value, old[[name]])
    }
    old[[name]] <- value
  }
  old
}
