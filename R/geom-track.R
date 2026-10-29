# Explicit annotation tracks: scale values are independent of physical bands.
ggchord_track_aes <- function(role, track) {
  if (!is.character(track) || length(track) != 1L || is.na(track) || !nzchar(track))
    ggchord_stop("Track scales require one nonempty track ID")
  paste0("track_", role, "_", paste(format(charToRaw(enc2utf8(track))), collapse = ""))
}

ggchord_track_stage_mapping <- function(mapping, track) {
  replacements <- c(fill = ggchord_track_aes("fill", track),
    colour = ggchord_track_aes("colour", track), color = ggchord_track_aes("colour", track),
    track_fill = ggchord_track_aes("fill", track), track_colour = ggchord_track_aes("colour", track))
  rewrite <- function(x, staged = FALSE) {
    if (is.symbol(x) && staged && as.character(x) %in% names(replacements))
      return(rlang::sym(replacements[[as.character(x)]]))
    if (!is.call(x)) return(x)
    active <- staged || identical(as.character(x[[1L]]), "after_scale")
    if (length(x) > 1L) for (i in seq.int(2L, length(x))) x[[i]] <- rewrite(x[[i]], active)
    x
  }
  for (name in names(mapping)) {
    q <- mapping[[name]]
    if (rlang::is_quosure(q)) mapping[[name]] <- rlang::quo_set_expr(q, rewrite(rlang::get_expr(q)))
  }
  mapping
}

#' Place a quantitative layer in an explicit sequence-local band
#'
#' Layers with the same track ID share a band and numerical scale across all
#' displayed sequences and windows. Their declarations must agree. Outside
#' means positive local normal (up in linear coordinates), inside negative.
#' The coordinate system resolves physical placement; feature-stack lanes and
#' sequence ring radii are unrelated to track IDs.
#' @param track One stable nonempty character ID.
#' @param side `"outside"` or `"inside"`.
#' @param order Finite ordering priority, nearest band first. Ties use ID order.
#' @param width,gap Positive band width and nonnegative preceding gap. Absolute
#'   grid units (mm, cm, inches or points) are resolved against the render
#'   context. Numeric values retain native layout units. Numeric per-sequence
#'   vectors/lists follow the usual sequence parameter formats.
#' @return A ggplot2 Position object for the three track geoms.
#' @md
#' @export
#' @examples
#' position_track("depth", width = unit(4, "mm"))
position_track <- function(track, side = c("outside", "inside"), order = 1,
    width = grid::unit(4, "mm"), gap = grid::unit(1, "mm")) {
  if (!is.character(track) || length(track) != 1L || is.na(track) || !nzchar(track))
    ggchord_stop("position_track(): track must be one nonempty ID")
  side <- match.arg(side)
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order)) ggchord_stop("position_track(): order must be finite")
  for (name in c("width", "gap")) {
    value <- if (name == "width") width else gap
    if (grid::is.unit(value)) {
      if (length(value) != 1L || !grid::unitType(value) %in% c("mm", "cm", "inches", "points", "bigpts"))
        ggchord_stop("position_track(): use one absolute physical unit, not npc/native units")
      number <- as.numeric(value)
      if (!is.finite(number) || number < 0 || (name == "width" && number == 0))
        ggchord_stop("position_track(): width must be positive and gap nonnegative")
    } else if (!is.numeric(unlist(value)) || any(!is.finite(unlist(value))) ||
        any(unlist(value) < 0) || (name == "width" && any(unlist(value) == 0)))
      ggchord_stop("position_track(): width must be positive and gap nonnegative")
  }
  ggplot2::ggproto(NULL, ggplot2::PositionIdentity, ggchord_track_position = TRUE,
    track = track, side = side, track_order = order, width = width, gap = gap)
}

#' Draw explicit quantitative tracks along sequence backbones
#'
#' Tile and bar layers require `accver`, `start`, `end`; lines require `accver`
#' and `pos`. Map `track_value` for numeric height and `track_fill` or
#' `track_colour` for visual mapping. Values are never inferred from column
#' names. `group` separates line series; positions must increase strictly
#' within each sequence/series. NA height breaks lines or omits bars, never
#' becomes zero. Inf is rejected. Tiles may overlap in input drawing order;
#' no averaging, interpolation or implicit aggregation is performed.
#'
#' All coordinates are window-local, one-based inclusive.
#' Interval bins occupy their full `end - start + 1` bases, including one-base
#' bins; line samples use the specified positions without binning.
#'
#' Use
#' [focus_ggchord_track()] for source windows. Lines never connect distinct
#' windows, even when they share a display row. Circular lines close only with
#' explicit `closed = TRUE`; missing endpoints are not joined across NA.
#'
#' @param mapping,data Standard ggplot2 mapping and layer input.
#' @param position A [position_track()] declaration. Every layer requires one.
#' @param fill,colour,alpha,linewidth,linetype Fixed styles; fill and colour
#'   mappings use the dedicated track roles and scoped scales.
#' @param baseline Explicit numeric bar baseline, default zero. It is included
#'   in scale training and must remain representable under limits/transform.
#' @param closed Explicitly close a complete circular line series.
#' @param bounds `"error"` (default) or `"clip"` for genomic coordinates.
#'   Clipping is reported and never changes source values.
#' @param show.legend,inherit.aes Standard layer arguments.
#' @param ... Extra geom parameters.
#' @return A ggplot2 layer with independent input and a shared named band.
#' @md
#' @export
#' @examples
#' s <- data.frame(accver = c("A", "B"), length = 1000)
#' d <- data.frame(accver = c("A", "B"), start = 100, end = 200, depth = c(10, 20))
#' ggchord(s) + geom_seq() +
#'   geom_track_bar(aes(track_value = depth), data = d,
#'     position = position_track("coverage")) +
#'   scale_track_value_continuous("coverage", limits = c(0, 30)) + coord_collinear()
geom_track_tile <- function(mapping = NULL, data = NULL, position,
    fill = "#86AFB9", colour = NA, alpha = 1, linewidth = .25,
    linetype = 1, bounds = c("error", "clip"), show.legend = NA, inherit.aes = FALSE, ...) {
  fixed <- c(fill = !missing(fill), colour = !missing(colour), alpha = !missing(alpha),
    linewidth = !missing(linewidth), linetype = !missing(linetype))
  ggchord_track_layer("tile", mapping, data, position, fill, colour, alpha,
    linewidth, linetype, bounds, 0, FALSE, show.legend, inherit.aes, fixed, list(...))
}
#' @rdname geom_track_tile
#' @export
geom_track_line <- function(mapping = NULL, data = NULL, position,
    colour = "#387D88", alpha = 1, linewidth = .4, linetype = 1,
    closed = FALSE, bounds = c("error", "clip"), show.legend = NA, inherit.aes = FALSE, ...) {
  fixed <- c(fill = FALSE, colour = !missing(colour), alpha = !missing(alpha),
    linewidth = !missing(linewidth), linetype = !missing(linetype))
  ggchord_track_layer("line", mapping, data, position, NA, colour, alpha,
    linewidth, linetype, bounds, 0, closed, show.legend, inherit.aes, fixed, list(...))
}
#' @rdname geom_track_tile
#' @export
geom_track_bar <- function(mapping = NULL, data = NULL, position,
    baseline = 0, fill = "#86AFB9", colour = NA, alpha = 1,
    linewidth = .25, linetype = 1, bounds = c("error", "clip"), show.legend = NA,
    inherit.aes = FALSE, ...) {
  fixed <- c(fill = !missing(fill), colour = !missing(colour), alpha = !missing(alpha),
    linewidth = !missing(linewidth), linetype = !missing(linetype))
  ggchord_track_layer("bar", mapping, data, position, fill, colour, alpha,
    linewidth, linetype, bounds, baseline, FALSE, show.legend, inherit.aes, fixed, list(...))
}

ggchord_track_layer <- function(kind, mapping, data, position, fill, colour,
    alpha, linewidth, linetype, bounds, baseline, closed, show.legend, inherit.aes, fixed, dots) {
  if (!inherits(position, "Position") || !isTRUE(position$ggchord_track_position)) ggchord_stop("Track geoms require position_track()")
  if (!is.numeric(baseline) || length(baseline) != 1L || !is.finite(baseline)) ggchord_stop("Track baseline must be finite")
  if (!is.logical(closed) || length(closed) != 1L || is.na(closed)) ggchord_stop("closed must be TRUE or FALSE")
  mapping <- ggchord_normalize_mapping(mapping) %||% ggplot2::aes()
  mapping <- ggchord_track_stage_mapping(mapping, position$track)
  if (any(c("x", "y") %in% names(mapping))) ggchord_stop("Track geometry maps accver/pos or start/end; do not override layout x/y")
  fixed[["colour"]] <- fixed[["colour"]] || "color" %in% names(dots)
  alias <- ggchord_colour_alias(colour, dots, "track geom", sys.call(-1L))
  colour <- alias$colour; dots <- alias$dots
  if ("fill" %in% names(mapping)) names(mapping)[names(mapping) == "fill"] <- "track_fill"
  if ("colour" %in% names(mapping)) names(mapping)[names(mapping) == "colour"] <- "track_colour"
  keys <- c(fill = ggchord_track_aes("fill", position$track), colour = ggchord_track_aes("colour", position$track))
  for (role in names(keys)) names(mapping)[names(mapping) == paste0("track_", role)] <- keys[[role]]
  geom <- ggplot2::ggproto(NULL, GeomGgchordTrack)
  names(geom$default_aes)[names(geom$default_aes) == "fill"] <- keys[["fill"]]
  names(geom$default_aes)[names(geom$default_aes) == "colour"] <- keys[["colour"]]
  defaults <- list(fill = fill, colour = colour, alpha = alpha,
    linewidth = linewidth, linetype = linetype)
  for (name in names(defaults)) {
    aesthetic <- if (name %in% names(keys)) keys[[name]] else name
    geom$default_aes[[aesthetic]] <- defaults[[name]]
  }
  styles <- defaults[fixed]
  for (role in names(keys)) names(styles)[names(styles) == role] <- keys[[role]]
  for (name in intersect(names(styles), names(mapping))) mapping[[name]] <- NULL
  layer <- ggplot2::layer(data = data.frame(x = numeric(), y = numeric(), group = integer()),
    mapping = ggplot2::aes(x = x, y = y, group = group, label = label,
      .component = .component), geom = geom,
    stat = "identity", position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, check.aes = FALSE,
    params = c(styles, list(track_keys = keys, track_kind = kind), dots))
  layer$ggchord_type <- "track"
  layer$ggchord_params <- list(type = "track", track_kind = kind, track_position = position,
    baseline = baseline, closed = closed, bounds = match.arg(bounds, c("error", "clip")), physical_layout = TRUE)
  ggchord_capture_layer_input(layer, data, mapping, c("accver", "start", "end", "pos", "track_value", "group"))
}

GeomGgchordTrack <- ggplot2::ggproto("GeomGgchordTrack", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(fill = "#86AFB9", colour = NA, alpha = 1,
    linewidth = .25, linetype = 1, size = 2.2,
    label = NA_character_, .component = "track_poly"),
  extra_params = c("na.rm", "track_keys", "track_kind", "axis_colour", "axis_size"),
  draw_panel = function(data, panel_params, coord, track_keys, track_kind,
      axis_colour = "#4F5962", axis_size = 2.2, na.rm = FALSE) {
    for (role in names(track_keys)) data[[role]] <- data[[track_keys[[role]]]]
    polys <- data[data$.component %in% "track_poly", , drop = FALSE]
    lines <- data[data$.component %in% "track_line", , drop = FALSE]
    axis <- data[data$.component %in% "track_axis_line", , drop = FALSE]
    text <- data[data$.component %in% "track_axis_text", , drop = FALSE]
    grobs <- list()
    if (nrow(polys)) grobs[[length(grobs) + 1L]] <- ggplot2::GeomPolygon$draw_panel(polys, panel_params, coord)
    if (nrow(lines)) grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(lines, panel_params, coord, na.rm = na.rm)
    if (nrow(axis)) {
      axis$colour <- axis_colour; axis$linewidth <- .25; axis$alpha <- 1; axis$linetype <- 1
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(axis, panel_params, coord, na.rm = na.rm)
    }
    if (nrow(text)) {
      text$colour <- axis_colour; text$alpha <- 1; text$size <- axis_size
      text$angle <- 0; text$hjust <- .5; text$vjust <- .5; text$family <- ""
      text$fontface <- 1; text$lineheight <- 1.2
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(text, panel_params, coord, na.rm = na.rm)
    }
    do.call(grid::grobTree, grobs)
  },
  draw_key = function(data, params, size) {
    for (role in names(params$track_keys)) data[[role]] <- data[[params$track_keys[[role]]]]
    if (params$track_kind == "line") ggplot2::draw_key_path(data, params, size) else ggplot2::draw_key_polygon(data, params, size)
  })

ggchord_prepare_track_layers <- function(plot) {
  plot$layers <- lapply(plot$layers, function(layer) {
    if (!identical(layer$ggchord_type, "track") || !isTRUE(layer$inherit.aes)) return(layer)
    mapping <- layer$ggchord_input_mapping
    inherited <- ggchord_normalize_mapping(plot$mapping)
    for (name in setdiff(names(inherited), c("x", "y"))) if (!name %in% names(mapping)) mapping[[name]] <- inherited[[name]]
    for (role in c("fill", "colour")) {
      key <- ggchord_track_aes(role, layer$ggchord_params$track_position$track)
      public <- paste0("track_", role)
      for (name in c(role, public)) if (name %in% names(mapping)) {
        if (!key %in% names(mapping)) mapping[[key]] <- mapping[[name]]
        mapping[[name]] <- NULL
      }
    }
    for (name in intersect(names(layer$aes_params), names(mapping))) mapping[[name]] <- NULL
    layer <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    layer$ggchord_input_mapping <- ggchord_track_stage_mapping(mapping, layer$ggchord_params$track_position$track)
    layer$inherit.aes <- FALSE
    layer
  })
  plot
}

#' Set a named track's numerical height scale
#'
#' Uses an independent continuous scale for each track ID. The same ID shares
#' its domain across sequences, windows and layers. The scale controls height,
#' not colour. A numerical axis is drawn by default; `guide = "none"` hides it.
#' @param track Track ID from [position_track()].
#' @param name Axis title, default the track ID.
#' @param limits,breaks,labels,transform,oob Standard continuous-scale controls.
#' @param unit Optional data-unit label retained in track exports.
#' @param guide `"axis"` or `"none"` for the numerical track axis.
#' @param ... Additional ggplot2 continuous-scale arguments.
#' @return A ggplot2 scale, scoped to this track.
#' @md
#' @export
#' @examples
#' scale_track_value_continuous("coverage", limits = c(0, 100), unit = "reads")
scale_track_value_continuous <- function(track, name = track, limits = NULL,
    breaks = ggplot2::waiver(), labels = ggplot2::waiver(), transform = "identity",
    oob = scales::censor, unit = NULL, guide = c("axis", "none"), ...) {
  position_track(track)
  if (!is.null(unit) && (!is.character(unit) || length(unit) != 1L || is.na(unit)))
    ggchord_stop("Track scale unit must be one character label or NULL")
  sc <- ggplot2::continuous_scale(ggchord_track_aes("value", track), palette = identity,
    name = name, limits = limits, breaks = breaks, labels = labels,
    transform = transform, oob = oob, guide = "none", ...)
  sc$ggchord_track_axis <- match.arg(guide) == "axis"
  sc$ggchord_track_unit <- unit
  sc$ggchord_track_id <- track
  sc
}

#' Set a track's fill or colour mapping
#' @param track Track ID.
#' @param ... Standard ggplot2 scale arguments.
#' @param colours Gradient colours.
#' @param values Manual category colours.
#' @param guide Optional standard guide for the fill/colour scale.
#' @return A standard scale scoped to the named track's visual role.
#' @md
#' @export
#' @examples
#' scale_track_fill_gradientn("gc", colours = c("navy", "white", "orange"))
scale_track_fill_gradientn <- function(track, ..., colours, guide = NULL) {
  key <- ggchord_track_aes("fill", track)
  ggplot2::scale_fill_gradientn(..., colours = colours, aesthetics = key,
    guide = guide %||% ggplot2::guide_colourbar(available_aes = key))
}
#' @rdname scale_track_fill_gradientn
#' @export
scale_track_fill_manual <- function(track, ..., values) {
  ggplot2::scale_fill_manual(..., values = values, aesthetics = ggchord_track_aes("fill", track))
}
#' @rdname scale_track_fill_gradientn
#' @export
scale_track_colour_gradientn <- function(track, ..., colours, guide = NULL) {
  key <- ggchord_track_aes("colour", track)
  ggplot2::scale_colour_gradientn(..., colours = colours, aesthetics = key,
    guide = guide %||% ggplot2::guide_colourbar(available_aes = key))
}
#' @rdname scale_track_fill_gradientn
#' @export
scale_track_colour_manual <- function(track, ..., values) {
  ggplot2::scale_colour_manual(..., values = values, aesthetics = ggchord_track_aes("colour", track))
}

ggchord_add_track_scales <- function(plot) {
  for (layer in plot$layers) {
    if (!identical(layer$ggchord_type, "track")) next
    mapping <- layer$ggchord_input_mapping
    for (role in c("fill", "colour")) {
      track <- layer$ggchord_params$track_position$track
      key <- ggchord_track_aes(role, track)
      if (!key %in% names(mapping) || plot$scales$has_scale(key)) next
      input <- layer$ggchord_resolved_input
      if (grepl("after_stat|after_scale", rlang::as_label(mapping[[key]]))) next
      values <- lapply(plot$layers, function(other) {
        if (!identical(other$ggchord_type, "track") || !key %in% names(other$ggchord_input_mapping)) return(NULL)
        q <- other$ggchord_input_mapping[[key]]
        if (grepl("after_stat|after_scale", rlang::as_label(q))) return(NULL)
        rlang::eval_tidy(q, other$ggchord_resolved_input)
      })
      values <- Filter(Negate(is.null), values)
      continuous <- vapply(values, is.numeric, logical(1))
      if (length(unique(continuous)) > 1L) ggchord_stop("Conflicting numerical and categorical visual mappings for track '", track, "'")
      value <- if (all(continuous)) unlist(values, use.names = FALSE) else
        unlist(lapply(values, as.character), use.names = FALSE)
      levels <- unique(as.character(value[!is.na(value)]))
      scale <- if (is.numeric(value)) ggplot2::scale_colour_gradientn(
        colours = c("#31507D", "#65B4A8", "#F0D56A"), aesthetics = key,
        guide = ggplot2::guide_colourbar(available_aes = key),
        name = rlang::as_label(mapping[[key]])) else ggplot2::scale_colour_manual(
          values = stats::setNames(chord_default_palette(length(levels)), levels),
          aesthetics = key, name = rlang::as_label(mapping[[key]]))
      plot$scales$add(scale)
    }
  }
  plot
}
