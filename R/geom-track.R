# Explicit annotation tracks: scale values are independent of physical bands.
ggchord_track_aes <- function(role, track) {
  track <- ggchord_track_id(track)
  paste0("track_", role, "_", paste(format(charToRaw(enc2utf8(track))), collapse = ""))
}

ggchord_track_stage_mapping <- function(mapping, track, replacements = NULL) {
  if (is.null(replacements)) replacements <- c(color = "colour")
  rewrite <- function(x, staged = FALSE) {
    if (is.symbol(x) && staged && as.character(x) %in% names(replacements))
      return(rlang::sym(replacements[[as.character(x)]]))
    if (!is.call(x)) return(x)
    if (rlang::is_call(x, "stage")) {
      x <- match.call(ggplot2::stage, x)
      if (!is.null(x$after_scale)) x$after_scale <- rewrite(x$after_scale, TRUE)
      return(x)
    }
    active <- staged || rlang::is_call(x, "after_scale")
    if (length(x) > 1L) for (i in seq.int(2L, length(x))) {
      if (!rlang::is_missing(x[[i]])) x[i] <- list(rewrite(x[[i]], active))
    }
    x
  }
  for (name in names(mapping)) {
    q <- mapping[[name]]
    if (rlang::is_quosure(q)) mapping[[name]] <- rlang::quo_set_expr(q, rewrite(rlang::get_expr(q)))
  }
  mapping
}

ggchord_track_id <- function(id) {
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(trimws(id)))
    ggchord_stop("Supply one nonempty track ID, for example track = 'depth'")
  id
}

ggchord_measure_name <- function(measure) {
  if (!is.null(measure) && (!is.character(measure) || length(measure) != 1L ||
      is.na(measure) || !nzchar(trimws(measure)))) ggchord_stop("measure must be one nonempty quantity name")
  measure
}

#' Declare a physical track band
#'
#' @param id Track identifier.
#' @param side Positive or negative side of the path.
#' @param order Order from the backbone. Tracks on the same side require distinct orders.
#' @param thickness,gap Band thickness and preceding gap in millimetres or absolute grid units.
#' @return A physical band declaration.
#' @export
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' signal <- data.frame(accver = rep(c("A", "B"), each = 3),
#'   pos = rep(c(100, 400, 700), 2), value = c(2, 5, 3, 1, 4, 2))
#'
#' ggchord(data = seqs) + geom_seq() +
#'   track_band("signal", thickness = 5, gap = 2) +
#'   geom_track_line(data = signal, aes(x = pos, y = value), track = "signal") +
#'   scale_track_value("signal", name = "Signal") + coord_collinear()
track_band <- function(id = "default", side = "positive", order = 1,
    thickness = 4, gap = 1) {
  ggchord_check_argument_names(sys.call(), sys.function())
  id <- ggchord_track_id(id)
  side <- match.arg(side, c("positive", "negative"))
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
      order < 1 || order != trunc(order))
    ggchord_stop("track_band(): order must be a positive integer")
  width <- ggchord_physical_mm(thickness, "track_band(thickness)", nonnegative = TRUE)
  gap <- ggchord_physical_mm(gap, "track_band(gap)", nonnegative = TRUE)
  if (any(unlist(width, use.names = FALSE) <= 0))
    ggchord_stop("track_band(): thickness must be positive")
  structure(list(id = id, side = if (side == "positive") "outside" else "inside",
    track_order = order, width = width, gap = gap), class = "ggchord_track_band")
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.ggchord_track_band <- function(object, plot, ...) {
  if (!inherits(plot, "ggchord")) ggchord_stop("Add track_band() to a ggchord plot")
  if (!is.null(plot$ggchord$track_bands[[object$id]]))
    ggchord_stop("Track band '", object$id, "' was already declared; declare each ID once")
  if (any(vapply(plot$ggchord$track_bands, function(band)
    band$side == object$side && band$track_order == object$track_order, logical(1))))
    ggchord_stop("track_band(): tracks on the same side must have distinct order values")
  plot$ggchord$track_bands[[object$id]] <- object
  plot$ggchord$ref <- new.env(parent = emptyenv())
  plot$ggchord$layout <- NULL
  plot
}

ggchord_uses_inherited_line_xy <- function(layer, mapping) {
  p <- layer$ggchord_params
  if (!isTRUE(layer$inherit.aes) || !identical(p$track_kind, "line") ||
      !is.null(p$native_geom) || isTRUE(p$closed)) return(FALSE)
  local <- ggchord_track_mapping_aliases(layer$ggchord_input_mapping)
  any(!intersect(c("x", "y"), names(mapping)) %in% names(local))
}

ggchord_native_track_line <- function(layer) {
  # The role-based line already holds the render's data callback, scoped
  # colour names and any legend key. Reuse that configuration for source x/y.
  keys <- layer$geom_params$track_keys
  public <- function(aesthetics) {
    for (role in names(keys)) names(aesthetics)[names(aesthetics) == keys[[role]]] <- role
    aesthetics
  }
  p <- layer$ggchord_params
  result <- do.call(geom_track, c(list(mapping = public(layer$ggchord_input_mapping),
    data = layer$ggchord_input_data, geom = "line", track = p$track_id, id = p$id,
    measure = p$measure, bounds = p$bounds, layout = layer$layout, show.legend = layer$show.legend,
    na.rm = layer$geom_params$na.rm), public(layer$aes_params)))
  result$name <- layer$name
  if (exists("draw_key", layer$geom, inherits = FALSE))
    result$geom$draw_key <- layer$geom$draw_key
  result
}
ggchord_track_layer <- function(kind, mapping, data, track, fill, colour,
    alpha, linewidth, linetype, bounds, baseline, closed, show.legend, inherit.aes, fixed, dots, id = NULL) {
  track <- ggchord_track_id(track %||% "default")
  if (!is.numeric(baseline) || length(baseline) != 1L || !is.finite(baseline)) ggchord_stop("Track baseline must be finite")
  mapping <- ggchord_normalize_mapping(mapping) %||% ggplot2::aes()
  mapping <- ggchord_track_stage_mapping(mapping, track)
  if (any(c("x", "y") %in% names(mapping))) ggchord_stop("Track geometry maps accver/pos or start/end; do not override layout x/y")
  fixed[["colour"]] <- fixed[["colour"]] || "color" %in% names(dots)
  alias <- ggchord_colour_alias(colour, dots, "track geom", sys.call(-1L))
  colour <- alias$colour; dots <- alias$dots
  keys <- c(fill = "fill", colour = "colour")
  geom <- ggplot2::ggproto(NULL, GeomGgchordTrack)
  names(geom$default_aes)[names(geom$default_aes) == "fill"] <- keys[["fill"]]
  names(geom$default_aes)[names(geom$default_aes) == "colour"] <- keys[["colour"]]
  defaults <- list(fill = fill, colour = colour, alpha = if (fixed[["alpha"]]) alpha else NA_real_,
    linewidth = linewidth, linetype = linetype)
  for (name in names(defaults)) {
    aesthetic <- if (name %in% names(keys)) keys[[name]] else name
    geom$default_aes[[aesthetic]] <- defaults[[name]]
  }
  styles <- defaults[fixed]
  fixed_names <- intersect(ggplot2::standardise_aes_names(names(dots)), geom$aesthetics())
  for (role in names(keys)) names(styles)[names(styles) == role] <- keys[[role]]
  local_styles <- names(mapping)
  for (role in names(keys)) local_styles[local_styles == keys[[role]]] <- role
  mapping[intersect(names(mapping), c(names(defaults)[fixed], fixed_names))] <- NULL
  ggchord_check_geom_params(dots, geom, paste0("geom_track_", kind, "()"))
  layer <- ggplot2::layer(data = data.frame(x = numeric(), y = numeric(), group = integer()),
    mapping = ggplot2::aes(x = x, y = y, group = group, label = label,
      .component = .component), geom = geom,
    stat = "identity", position = "identity", show.legend = show.legend,
    inherit.aes = inherit.aes, check.aes = FALSE,
    params = c(styles, list(track_keys = keys, track_kind = kind), dots))
  layer$ggchord_type <- "track"
  layer$ggchord_params <- list(type = "track", id = id, track_kind = kind, track_id = track,
    baseline = baseline, closed = closed, bounds = match.arg(bounds, c("error", "clip")), physical_layout = TRUE)
  ggchord_capture_layer_input(layer, data, mapping, c("accver", "start", "end", "pos", "track_value", "group"))
}

GeomGgchordTrack <- ggplot2::ggproto("GeomGgchordTrack", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(fill = "#86AFB9", colour = NA, alpha = 1,
    linewidth = .25, linetype = 1, size = 2.2, shape = 19, stroke = .5,
    label = NA_character_, angle=0, hjust=.5, vjust=.5, family="",fontface=1,lineheight=1.2, .component = "track_poly"),
  extra_params = c("na.rm", "track_keys", "track_kind", "axis_colour", "axis_size", "native_draw_params",
    "axis_line_params", "axis_tick_params", "axis_text_params", "axis_title_params", "axis_check_overlap"),
  draw_panel = function(data, panel_params, coord, track_keys, track_kind,
      axis_colour = NULL, axis_size = NULL, na.rm = FALSE, native_draw_params=list(),
      axis_line_params = list(), axis_tick_params = list(), axis_text_params = list(), axis_title_params = list(), axis_check_overlap = FALSE) {
    for (role in names(track_keys)) data[[role]] <- data[[track_keys[[role]]]]
    polys <- data[data$.component %in% c("track_poly", "track_area_fill"), , drop = FALSE]
    polys$colour[polys$.component == "track_area_fill"] <- NA
    lines <- data[data$.component %in% "track_line", , drop = FALSE]
    points <- data[data$.component %in% "track_point", , drop = FALSE]
    labels <- data[data$.component %in% "track_text", , drop = FALSE]
    axis <- data[data$.component %in% "track_axis_line", , drop = FALSE]
    tick <- data[data$.component %in% "track_axis_tick", , drop = FALSE]
    text <- data[data$.component %in% "track_axis_text", , drop = FALSE]
    title <- data[data$.component %in% "track_axis_title", , drop = FALSE]
    grobs <- list()
    draw <- function(geom,d,keys=character()) {
      args <- list(data=d,panel_params=panel_params,coord=coord)
      if (!identical(geom,ggplot2::GeomPolygon)) args$na.rm <- na.rm
      do.call(geom$draw_panel,c(args,native_draw_params[intersect(names(native_draw_params),keys)]))
    }
    if (nrow(polys)) grobs[[length(grobs) + 1L]] <- draw(ggplot2::GeomPolygon,polys,"rule")
    if (nrow(lines)) grobs[[length(grobs) + 1L]] <- draw(ggplot2::GeomPath,lines,c("arrow","lineend","linejoin","linemitre"))
    if (nrow(points)) grobs[[length(grobs) + 1L]] <- ggplot2::GeomPoint$draw_panel(points, panel_params, coord, na.rm = na.rm)
    if (nrow(labels)) {
      labels$angle <- labels$angle %||% 0
      labels$hjust <- labels$hjust %||% .5; labels$vjust <- labels$vjust %||% .5
      labels$family <- labels$family %||% ""; labels$fontface <- labels$fontface %||% 1
      labels$lineheight <- labels$lineheight %||% 1.2
      grobs[[length(grobs) + 1L]] <- draw(ggplot2::GeomText,labels,c("parse","check_overlap"))
    }
    style_data <- function(d, style) {
      for (name in setdiff(names(style), "hidden")) if (length(style[[name]])) d[[name]] <- style[[name]]
      if (!is.null(axis_colour)) d$colour <- axis_colour
      d$alpha <- 1
      d
    }
    if (nrow(axis) && !isTRUE(axis_line_params$hidden)) {
      axis <- style_data(axis, axis_line_params)
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(axis, panel_params, coord, na.rm = na.rm,
        lineend = axis_line_params$lineend %||% "butt", linejoin = axis_line_params$linejoin %||% "round")
    }
    if (nrow(tick) && !isTRUE(axis_tick_params$hidden)) {
      tick <- style_data(tick, axis_tick_params)
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomSegment$draw_panel(
        tick, panel_params, coord, na.rm = na.rm)
    }
    for (component in c("text", "title")) {
      d <- if (component == "text") text else title
      style <- if (component == "text") axis_text_params else axis_title_params
      if (!nrow(d) || isTRUE(style$hidden)) next
      d <- style_data(d, style)
      if (!is.null(axis_size)) d$size <- axis_size
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(d, panel_params, coord, na.rm = na.rm,
        check_overlap = component == "text" && axis_check_overlap)
    }
    do.call(grid::grobTree, grobs)
  },
  draw_key = function(data, params, size) {
    for (role in names(params$track_keys)) data[[role]] <- data[[params$track_keys[[role]]]]
    if (params$track_kind %in% c("line", "path", "step", "segment", "linerange", "errorbar")) ggplot2::draw_key_path(data, params, size) else
      if (params$track_kind == "point") ggplot2::draw_key_point(data, params, size) else
        if (params$track_kind == "text") ggplot2::draw_key_text(data, params, size) else ggplot2::draw_key_polygon(data, params, size)
  })

# A single render owns the source-domain statistic. Guides and panel-size
# refinement may rebuild geometry, but use the same observations and groups.
ggchord_render_track_source <- function(plot) {
  pending <- vapply(plot$layers, function(layer)
    (!is.null(layer$ggchord_params$native_geom) && !isTRUE(layer$ggchord_params$native_prepared)) ||
      ggchord_uses_inherited_line_xy(layer, plot$mapping), logical(1))
  if (any(pending)) {
    plot <- ggchord_register_layers(plot)
    plot <- ggchord_prepare_facet_sources(plot)
    sequences <- ggchord_sequence_reference_data(plot)
    for (i in which(pending))
      plot$layers[[i]] <- ggchord_track_native_input(plot$layers[[i]], plot, sequences)
  }
  plot
}

ggchord_prepare_track_layers <- function(plot, sequences = ggchord_sequence_reference_data(plot)) {
  plot$layers <- lapply(plot$layers, function(layer) {
    if (!identical(layer$ggchord_type, "track")) return(layer)
    if (!is.null(layer$ggchord_params$native_geom)) {
      if (!isTRUE(layer$ggchord_params$native_prepared)) return(ggchord_track_native_input(layer, plot, sequences))
      return(layer)
    }
    source <- ggchord_layer_source(layer, plot$ggchord$data$seq_data)
    caller <- ggchord_layer_input_context(layer)
    ggchord_check_interval_unit(source, ggchord_sequence_unit(sequences), caller)
    mapping <- layer$ggchord_input_mapping
    inherited <- if (isTRUE(layer$inherit.aes)) ggchord_normalize_mapping(plot$mapping)
    if (any(c("y", "ymin", "ymax") %in% names(inherited)))
      warning("geom_track_interval(): ignored inherited y aesthetics", call. = FALSE)
    for (name in setdiff(names(inherited), c("x", "y", "ymin", "ymax")))
      if (!name %in% names(mapping)) mapping[[name]] <- inherited[[name]]
    for (name in intersect(names(layer$aes_params), names(mapping))) mapping[[name]] <- NULL
    layer <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    layer$ggchord_input_mapping <- ggchord_track_stage_mapping(mapping, layer$ggchord_params$track_id)
    layer$inherit.aes <- FALSE
    input <- ggchord_resolve_layer_input(layer, plot$ggchord$data$seq_data, data = source)
    ggchord_track_reference_index(input, sequences, caller)
    capture <- layer$ggchord_params$stage_capture
    if (!is.null(capture) && "statistic" %in% capture$requested)
      capture$statistic <- input
    layer$ggchord_input_data <- input
    layer$ggchord_params$track_sources <- source
    for (role in intersect(layer$ggchord_role_aes, names(mapping))) layer$ggchord_input_mapping[[role]] <- NULL
    layer
  })
  numeric_layers <- Filter(function(l) identical(l$ggchord_type, "track") &&
    !identical(l$ggchord_params$track_kind, "tile"), plot$layers)
  ids <- vapply(numeric_layers, function(l) l$ggchord_params$track_id, character(1))
  plot$ggchord$track_quantities <- setNames(lapply(unique(ids), function(id) {
    members <- numeric_layers[ids == id]
    known <- function(field) unique(unlist(lapply(members, function(l) l$ggchord_params[[field]])))
    measures <- known("measure")
    scale <- plot$scales$get_scales(ggchord_track_aes("value", id))
    units <- unique(c(known("value_unit"), scale$ggchord_track_unit))
    if (length(units) > 1L)
      ggchord_stop("Track '", id, "' mixes incompatible value units")
    if (length(measures) > 1L) {
      owners <- vapply(members, function(l) l$ggchord_layer_id, character(1))
      ggchord_stop("Track '", id, "' mixes height measures: ",
        paste(owners, collapse = ", "), ". Use separate bands for different quantities.")
    }
    labels <- known("value_label")
    list(unit = units, measure = measures,
      label = if (length(measures)) measures else if (length(labels) == 1L) labels else id)
  }), unique(ids))
  plot
}

#' @noRd
scale_track_value_continuous <- function(track, name = track, limits = NULL,
    breaks = ggplot2::waiver(), labels = ggplot2::waiver(), transform = "identity",
    oob = scales::censor, unit = NULL, guide = c("axis", "none"), ...,
    facet_domain = c("shared", "panel")) {
  ggchord_check_argument_names(sys.call(), sys.function(), ggplot2::continuous_scale)
  if (!is.null(unit) && (!is.character(unit) || length(unit) != 1L || is.na(unit)))
    ggchord_stop("Track scale unit must be one character label or NULL")
  sc <- ggplot2::continuous_scale(ggchord_track_aes("value", track), palette = identity,
    name = name, limits = limits, breaks = breaks, labels = labels,
    transform = transform, oob = oob, guide = "none", ...)
  axis <- ggchord_track_axis_guide(guide)
  sc$ggchord_facet_domain <- match.arg(facet_domain)
  if (sc$ggchord_facet_domain == "panel" && !axis$show)
    ggchord_stop("Panel-local track domains require a visible numeric axis")
  sc$ggchord_domain_policy <- list(transform = transform, oob = oob, limits = limits,
    na.value = sc$na.value, rescaler = get("rescaler", envir = sc, inherits = FALSE))
  sc$ggchord_track_axis <- axis$show
  sc$ggchord_track_axis_guide <- axis$guide
  sc$ggchord_track_unit <- unit
  sc$ggchord_track_id <- track
  sc
}

#' @noRd
shared_track_limits <- function(..., tracks = "default") {
  ggchord_check_argument_names(sys.call(), sys.function())
  plots <- list(...)
  if (!length(plots) || !all(vapply(plots, inherits, logical(1), "ggchord")) ||
      !is.character(tracks) || !length(tracks) %in% c(1L, length(plots)))
    ggchord_stop("Supply ggchord plots and one track ID per plot or one shared ID")
  tracks <- rep_len(tracks, length(plots))
  policies <- domains <- vector("list", length(plots))
  for (i in seq_along(plots)) {
    id <- ggchord_track_id(tracks[i])
    layout <- get_chord_layout(ggchord_render_plot_copy(plots[[i]]))
    registry <- layout$track_registry
    band <- registry[registry$track_id == id, , drop = FALSE]
    scale <- layout$track_scales[[id]]
    if (is.null(scale) || !NROW(band) || is.null(band$measure) || anyNA(band$measure))
      ggchord_stop("Plot ", i, ", track '", id, "' needs a numeric track with a declared measure")
    if (is.null(scale$ggchord_domain_policy))
      ggchord_stop("Reapply scale_track_value_continuous() to track '", id, "' in plot ", i,
        " to record its scale policy before sharing limits")
    policies[[i]] <- list(measure = unique(band$measure), unit = unique(band$unit),
      coordinate_unit = layout$sequence_unit, scale = scale$ggchord_domain_policy)
    if (!identical(policies[[i]], policies[[1L]]))
      ggchord_stop("Incompatible measure, unit or scale policy in plot ", i, ", track '", id, "'")
    domains[[i]] <- unlist(band$limits, use.names = FALSE)
  }
  values <- unlist(domains, use.names = FALSE)
  if (!length(values)) ggchord_stop("The selected tracks have no numerical domain")
  range(values)
}

#' @noRd
scale_track_fill_gradientn <- function(track, ..., colours, guide = NULL) {
  key <- ggchord_track_aes("fill", track)
  ggplot2::scale_fill_gradientn(..., colours = colours, aesthetics = key,
    guide = guide %||% ggplot2::guide_colourbar(available_aes = key))
}
#' @noRd
scale_track_fill_manual <- function(track, ..., values) {
  ggplot2::scale_fill_manual(..., values = values, aesthetics = ggchord_track_aes("fill", track))
}
#' @noRd
scale_track_colour_gradientn <- function(track, ..., colours, guide = NULL) {
  key <- ggchord_track_aes("colour", track)
  ggplot2::scale_colour_gradientn(..., colours = colours, aesthetics = key,
    guide = guide %||% ggplot2::guide_colourbar(available_aes = key))
}
#' @noRd
scale_track_colour_manual <- function(track, ..., values) {
  ggplot2::scale_colour_manual(..., values = values, aesthetics = ggchord_track_aes("colour", track))
}

ggchord_add_track_scales <- function(plot) {
  seen <- character()
  for (layer in plot$layers) {
    # ggplot2 takes its default labels from the first mapped layer. Use the
    # original expression there, before projection introduced .track_* names.
    mapped <- names(layer$mapping)
    if (isTRUE(layer$inherit.aes)) mapped <- union(mapped, names(plot$mapping))
    first <- setdiff(mapped, c(seen, names(layer$aes_params)))
    seen <- union(seen, first)
    labels <- layer$ggchord_params$track_mapping_labels
    for (key in intersect(first, names(labels))) {
      if (!key %in% c("shape", layer$geom_params$track_keys)) next
      name <- labels[[key]]
      if (is.function(plot$labels[[key]])) {
        plot$labels[[key]] <- plot$labels[[key]](name)
      } else if (!key %in% names(plot$labels)) {
        dict <- plot$labels$dictionary
        plot$labels[[key]] <- if (name %in% names(dict)) unname(dict[name]) else name
      }
    }
    if (!identical(layer$ggchord_type, "track")) next
    mapping <- layer$ggchord_input_mapping
    for (role in c("fill", "colour")) {
      track <- layer$ggchord_params$track_id
      key <- ggchord_track_aes(role, track)
      if (!key %in% names(mapping) || plot$scales$has_scale(key)) next
      if (is.null(ggchord_start_mapping(mapping[[key]]))) next
      values <- lapply(plot$layers, function(other) {
        if (!identical(other$ggchord_type, "track") || !key %in% names(other$ggchord_input_mapping)) return(NULL)
        if (!NROW(other$data)) return(NULL)
        columns <- other$ggchord_params$visual_columns
        if (key %in% names(columns)) return(other$data[[columns[[key]]]])
        q <- ggchord_start_mapping(other$ggchord_input_mapping[[key]])
        if (is.null(q)) return(NULL)
        rlang::eval_tidy(q, other$data)
      })
      values <- Filter(Negate(is.null), values)
      continuous <- vapply(values, is.numeric, logical(1))
      if (length(unique(continuous)) > 1L) ggchord_stop("Conflicting numerical and categorical visual mappings for track '", track, "'")
      value <- if (all(continuous)) unlist(values, use.names = FALSE) else
        unlist(lapply(values, as.character), use.names = FALSE)
      levels <- unique(as.character(value[!is.na(value)]))
      scale <- if (is.numeric(value)) ggplot2::scale_colour_gradientn(
        colours = c("#31507D", "#65B4A8", "#F0D56A"), aesthetics = key,
        guide = ggplot2::guide_colourbar(available_aes = key)) else ggplot2::scale_colour_manual(
          values = stats::setNames(chord_default_palette(length(levels)), levels),
          aesthetics = key, guide=ggchord_role_guide(plot,"track",labels=levels))
      plot$scales$add(scale)
    }
  }
  plot
}
