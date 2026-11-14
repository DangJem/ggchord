#' Choose target-size styles and joint guide layout
#'
#' Add this policy to a ggchord plot. The ordinary print, ggsave, preview and
#' diagnostic paths rebuild at the actual device or Fusion panel size.
#' Plots without a policy retain their existing behaviour.
#'
#' `fixed` keeps default theme sizes and guide furniture at their reference
#' physical sizes. `scale` adjusts theme-driven default text, lines and spacing
#' within bounds. `adaptive` additionally selects a joint arrangement of trained
#' native legends and colour bars, then remeasures the complete plot. Existing
#' label layout is rerun with the resulting sizes; this is not a new global
#' label collision solver. Data, source coordinates and scale mappings are unchanged.
#' Omitted fixed line widths in native backbone, link, ribbon and alignment
#' layers follow the bounded line factor. Explicit widths and mapped widths
#' retain their own values and scale rules. Other geometry constants and native
#' theme-derived aesthetic expressions are not uniformly multiplied.
#'
#' Explicit geom styles, mapped sizes, theme fields and guide constraints take
#' priority. Only unspecified guide grids, directions and positions can change.
#' Continuous/stepped bars stay intact and retain their scale breaks and order.
#' Unknown custom guides are measured but not rearranged. Infeasible layouts
#' warn and report their constraints; categories are never silently removed.
#'
#' Apply a policy to each desired Fusion panel before composition. Cross-panel
#' endpoints are rebuilt after panel measurement; the outer Fusion allocations
#' and its separate relationship guide are not controlled by this policy.
#' Facets and arbitrary annotation grobs are outside the automatic layout contract.
#'
#' @param mode One of `adaptive`, `fixed`, or `scale`.
#' @param reference_size Positive width and height in inches. Scale modes use
#'   the fourth root of target/reference area, with a bounded text multiplier.
#' @param text_scale_range Lower and upper positive text multipliers containing
#'   `1`. Lines and default spacing use the square root of the text multiplier.
#' @param min_text_size Lower bound in mm for automatically reduced theme text.
#'   Defaults already smaller than this value are preserved rather than enlarged.
#'   Explicit text sizes and existing label-fit policies keep their own contracts.
#' @return A plot policy added with `+`. Target-size [diagnose_ggchord()] results
#'   include `adaptive_layout` (effective styles and scale factors) and
#'   `guide_layout` (candidates, constraints, chosen layouts and measured bounds).
#'   A feasible guide layout is not certification of every label or arbitrary
#'   grob. Ordinary layout-only exports report the policy without claiming a
#'   target-device guide arrangement.
#' @export
#' @examples
#' library(dplyr)
#' seqs <- seq_data_example |> select(accver, length)
#' p <- ggchord(seqs, ribbon_data_example) +
#'   geom_seq() + geom_link_ribbon() +
#'   adapt_ggchord(reference_size = c(8, 6))
#' p
#' # ggplot2::ggsave("adaptive.png", p, width = 8, height = 6, units = "in")
adapt_ggchord <- function(mode = c("adaptive", "fixed", "scale"),
    reference_size = c(8, 6), text_scale_range = c(.85, 1.35), min_text_size = 2.5) {
  mode <- match.arg(mode)
  if (!is.numeric(reference_size) || length(reference_size) != 2L ||
      any(!is.finite(reference_size)) || any(reference_size <= 0))
    ggchord_stop("adapt_ggchord(): reference_size must be positive width and height in inches")
  if (!is.numeric(text_scale_range) || length(text_scale_range) != 2L ||
      any(!is.finite(text_scale_range)) || any(text_scale_range <= 0) ||
      text_scale_range[1] > 1 || text_scale_range[2] < 1)
    ggchord_stop("adapt_ggchord(): text_scale_range must contain positive lower/upper multipliers around 1")
  if (!is.numeric(min_text_size) || length(min_text_size) != 1L ||
      !is.finite(min_text_size) || min_text_size < 0)
    ggchord_stop("adapt_ggchord(): min_text_size must be one nonnegative size in mm")
  structure(list(mode=mode,reference_size=unname(reference_size),
    text_scale_range=unname(text_scale_range),min_text_size=min_text_size),class="ggchord_layout_policy")
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.ggchord_layout_policy <- function(object, plot, ...) {
  if (!inherits(plot,"ggchord")) ggchord_stop("adapt_ggchord(): add the policy to an individual ggchord plot")
  if (!inherits(plot$facet,"FacetNull")) ggchord_stop("adapt_ggchord(): faceted layouts are not supported")
  plot$ggchord$layout_policy <- unclass(object)
  plot$ggchord$ref <- new.env(parent=emptyenv())
  plot$ggchord$layout <- NULL
  plot
}

# Build-local policy: data and the source plot stay untouched.
ggchord_adaptive_factors <- function(policy, canvas = ggchord_render_dimensions()) {
  measured <- length(canvas) == 2L && all(is.finite(canvas)) && all(canvas > 0)
  factor <- if (policy$mode == "fixed" || !measured) 1 else
    (prod(canvas / 25.4) / prod(policy$reference_size))^.25
  factor <- max(policy$text_scale_range[1], min(policy$text_scale_range[2], factor))
  if (policy$mode == "fixed") factor <- 1
  list(mode = policy$mode, canvas_mm = canvas, reference_size_in = policy$reference_size,
    text_factor = factor, line_factor = sqrt(factor), spacing_factor = sqrt(factor),
    min_text_size_mm = policy$min_text_size,
    measurement = if (measured) "device_dimensions" else "unavailable")
}

ggchord_apply_adaptive_styles <- function(plot) {
  policy <- plot$ggchord$layout_policy
  if (is.null(policy)) return(plot)
  if (!inherits(plot$facet, "FacetNull"))
    ggchord_stop("adapt_ggchord(): faceted layouts are not supported")
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
      if (field == "size") target <- max(target, min(value, policy$min_text_size * ggplot2::.pt))
    }
    raw[[field]] <- target
    theme[[name]] <- raw
    sizes[[length(sizes) + 1L]] <- data.frame(element = name, property = field,
      requested = value, effective = target, user_locked = locked,
      units = if (field == "size") "pt" else "mm")
  }
  for (name in c("plot.margin", "legend.spacing", "legend.spacing.x", "legend.spacing.y", "legend.box.spacing"))
    if (grid::is.unit(theme[[name]]) && is.null(explicit[[name]]))
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
