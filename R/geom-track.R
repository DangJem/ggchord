# Explicit annotation tracks: scale values are independent of physical bands.
ggchord_track_aes <- function(role, track) {
  track <- ggchord_track_id(track)
  paste0("track_", role, "_", paste(format(charToRaw(enc2utf8(track))), collapse = ""))
}

ggchord_track_stage_mapping <- function(mapping, track, replacements = NULL) {
  if (is.null(replacements)) replacements <- c(fill = ggchord_track_aes("fill", track),
    colour = ggchord_track_aes("colour", track), color = ggchord_track_aes("colour", track),
    track_fill = ggchord_track_aes("fill", track), track_colour = ggchord_track_aes("colour", track))
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

#' Declare a quantitative band along sequence backbones
#'
#' Add a band to a ggchord plot and reference its ID from any number of
#' quantitative layers, before or after its declaration. Omit `track` on a
#' layer to use the shared `"default"` band. Declare that ID explicitly to
#' change its placement. Layers on the same band share a numerical scale
#' across displayed sequences and windows. Use separate IDs for different units.
#' @param id One nonempty character ID.
#' @param side `"outside"` (positive local normal, up in linear coordinates)
#'   or `"inside"`.
#' @param order Finite ordering priority, nearest band first. Ties use ID order.
#' @param width,gap Positive band width and nonnegative preceding gap. Supply
#'   an absolute grid unit (mm, cm, inches or points), or native layout units
#'   as a scalar or a per-sequence numeric vector/list.
#' @return A band declaration to add with `+`.
#' @md
#' @export
#' @examples
#' observations <- data.frame(accver = "A", pos = c(25, 75), depth = c(3, 5))
#' ggchord(data.frame(accver = "A", length = 100)) +
#'   track_band("depth", width = grid::unit(8, "mm")) +
#'   geom_track_col(aes(pos, depth), data = observations, track = "depth", width = 10) +
#'   geom_track(aes(pos, depth), data = observations, geom = "point", track = "depth") +
#'   coord_collinear()
track_band <- function(id, side = c("outside", "inside"), order = 1,
    width = grid::unit(4, "mm"), gap = grid::unit(1, "mm")) {
  ggchord_check_argument_names(sys.call(), sys.function())
  id <- ggchord_track_id(id)
  side <- match.arg(side)
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order)) ggchord_stop("track_band(): order must be finite")
  for (name in c("width", "gap")) {
    value <- if (name == "width") width else gap
    if (grid::is.unit(value)) {
      if (length(value) != 1L || !grid::unitType(value) %in% c("mm", "cm", "inches", "points", "bigpts"))
        ggchord_stop("track_band(): use one absolute physical unit for ", name)
      numbers <- as.numeric(value)
    } else {
      numbers <- unlist(value, use.names = FALSE)
    }
    if (!length(numbers) || !is.numeric(numbers) || any(!is.finite(numbers)) ||
        any(numbers < 0) || (name == "width" && any(numbers == 0)))
      ggchord_stop("track_band(): width must be positive and gap nonnegative")
  }
  structure(list(id = id, side = side, track_order = order, width = width, gap = gap),
    class = "ggchord_track_band")
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.ggchord_track_band <- function(object, plot, ...) {
  if (!inherits(plot, "ggchord")) ggchord_stop("Add track_band() to a ggchord plot")
  if (!is.null(plot$ggchord$track_bands[[object$id]]))
    ggchord_stop("Track band '", object$id, "' was already declared; declare each ID once")
  plot$ggchord$track_bands[[object$id]] <- object
  plot$ggchord$ref <- new.env(parent = emptyenv())
  plot$ggchord$layout <- NULL
  plot
}

#' Draw explicit quantitative tracks along sequence backbones
#'
#' Tiles require `accver`, `start`, `end`; lines require `accver`
#' and `pos`. Map `track_value` for numeric height and `track_fill` or
#' `track_colour` for visual mapping. Values are never inferred from column
#' names. `group` separates line series; positions must increase strictly
#' within each sequence/series. NA height breaks lines and never
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
#' @param track Character ID declared with [track_band()]. Omit for the default band.
#' @param id Optional unique layer ID used in layout and source exports.
#' @param position ggplot2 data position. Full-band interval tiles and closed
#'   lines use identity. For tiles with local numerical heights and data
#'   adjustments, use `geom_track(geom = "tile")`. Open lines with local or
#'   inherited x/y mappings, an explicit stat, or a data adjustment use the native line
#'   grammar. Lines with sequence-local roles and identity position preserve
#'   input order, missing-value breaks and strict increasing positions.
#' @param fill,colour,alpha,linewidth,linetype Fixed styles; fill and colour
#'   mappings use the dedicated track roles and scoped scales.
#'   A fixed style overrides its inherited mapping. Choose either a mapping
#'   or a fixed value for that aesthetic within a layer.
#' @param stat Optional source statistic for open lines, or `"coverage"` for full-band tiles.
#' @inheritParams geom_track
#' @param closed Explicitly close a complete circular line series.
#' @param bounds `"error"` (default) or `"clip"` for genomic coordinates.
#'   Clipping is reported and never changes source values.
#' @param ... Extra geom parameters.
#' @return A ggplot2 layer with independent input and a shared named band.
#' @md
#' @export
#' @examples
#' library(ggchord)
#' library(dplyr)
#' data(plasmid_example_pUC19)
#' track_seq <- plasmid_example_pUC19 |>
#'   select(accver, length)
#' starts <- seq.int(1L, track_seq$length, by = 100L)
#' gc_fraction <- function(sequence) {
#'   bases <- strsplit(sequence, "")[[1]]
#'   mean(bases %in% c("G", "C"))
#' }
#' bins <- tibble(
#'   accver = track_seq$accver,
#'   start = starts,
#'   end = pmin(starts + 99L, track_seq$length)
#' ) |>
#'   rowwise() |>
#'   mutate(gc = gc_fraction(substr(plasmid_example_pUC19$sequence, start, end))) |>
#'   ungroup()
#' doc_theme <- theme_ggchord_publication(
#'   axis.text = element_text(size = 9), axis.text.check.overlap = TRUE
#' ) +
#'   theme(legend.position = "bottom")
#' p_tracks <- ggchord(data = track_seq) +
#'   geom_seq() +
#'   geom_seq_center_label() +
#'   track_band("gc", width = unit(6, "mm")) +
#'   geom_track_col(aes(track_value = gc, track_fill = gc),
#'     data = bins,
#'     track = "gc"
#'   ) +
#'   scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
#'   scale_track_fill_gradientn("gc", colours = c("#EDF1EA", "#8AAE9C", "#345F58")) +
#'   coord_circular() +
#'   doc_theme
#' ggchord(track_seq) + geom_track_line(aes(x = (start + end) / 2, y = gc), data = bins) + coord_collinear()
#' ggchord(track_seq) +
#'   geom_track_tile(aes(fill = after_stat(coverage)), data = bins,
#'     stat = "coverage", binwidth = 250) +
#'   coord_collinear()
geom_track_tile <- function(mapping = NULL, data = NULL, position = NULL,
    fill = "#86AFB9", colour = NA, alpha = 1, linewidth = .25,
    linetype = 1, bounds = c("error", "clip"), show.legend = NA, inherit.aes = TRUE, ..., track = NULL, id = NULL, stat = "identity") {
  ggchord_check_argument_names(sys.call(), sys.function())
  fixed <- c(fill = !missing(fill), colour = !missing(colour), alpha = !missing(alpha),
    linewidth = !missing(linewidth), linetype = !missing(linetype))
  if (!ggchord_is_data_identity(position %||% "identity"))
    ggchord_stop("Full-band tiles require identity position; use geom_track(geom = 'tile') for data adjustments")
  if (!identical(stat, "identity")) {
    if (!identical(stat, "coverage")) ggchord_stop("Full-band tiles accept identity or coverage statistics")
    dots <- list(...)
    for (key in names(fixed)[fixed]) dots[[key]] <- get(key)
    return(do.call(geom_track, c(list(mapping = mapping, data = data, geom = "tile", stat = stat,
      track = track, position = position, bounds = match.arg(bounds), show.legend = show.legend,
      inherit.aes = inherit.aes, id = id), dots)))
  }
  ggchord_track_layer("tile", mapping, data, track, fill, colour, alpha,
    linewidth, linetype, bounds, 0, FALSE, show.legend, inherit.aes, fixed, list(...), id = id)
}
#' @rdname geom_track_tile
#' @export
geom_track_line <- function(mapping = NULL, data = NULL, position = NULL,
    colour = "#387D88", alpha = 1, linewidth = .4, linetype = 1,
    closed = FALSE, bounds = c("error", "clip"), show.legend = NA, inherit.aes = TRUE, ..., track = NULL, id = NULL, stat = NULL, measure = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  fixed <- c(fill = FALSE, colour = !missing(colour), alpha = !missing(alpha),
    linewidth = !missing(linewidth), linetype = !missing(linetype))
  if (!is.logical(closed) || length(closed) != 1L || is.na(closed)) ggchord_stop("closed must be TRUE or FALSE")
  dots <- list(...)
  identity <- ggchord_is_data_identity(position %||% "identity")
  if (closed && (!identity || !is.null(stat)))
    ggchord_stop("Closed lines require identity position and no computed statistic")
  native <- any(c("x", "y") %in% names(mapping)) || !identity || !is.null(stat)
  if (native && !closed) {
    for (key in names(fixed)[fixed]) dots[[key]] <- get(key)
    return(do.call(geom_track, c(list(mapping = mapping, data = data, geom = "line", stat = stat,
      track = track, position = position, bounds = match.arg(bounds), show.legend = show.legend,
      inherit.aes = inherit.aes, id = id, measure = measure), dots)))
  }
  if (closed && !is.null(mapping)) {
    for (key in intersect(c("x", "y"), names(mapping))) {
      role <- if (key == "x") "pos" else "track_value"
      if (role %in% names(mapping)) ggchord_stop("Specify only one of ", key, " and ", role)
      names(mapping)[names(mapping) == key] <- role
    }
  }
  layer <- ggchord_track_layer("line", mapping, data, track, NA, colour, alpha,
    linewidth, linetype, bounds, 0, closed, show.legend, inherit.aes, fixed, list(...), id = id)
  layer$ggchord_params$measure <- ggchord_measure_name(measure)
  layer
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
#' @rdname geom_track
#' @export
geom_track_bar <- function(mapping = NULL, data = NULL, position = NULL,
    baseline = 0, fill = "#86AFB9", colour = NA, alpha = 1,
    linewidth = .25, linetype = 1, bounds = c("error", "clip"), show.legend = NA,
    inherit.aes = TRUE, ..., track = NULL, stat = "count", id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!is.numeric(baseline) || length(baseline) != 1L || !is.finite(baseline))
    ggchord_stop("geom_track_bar(): baseline must be one finite number")
  if (baseline != 0) ggchord_stop("geom_track_bar(): count uses a zero baseline; use geom_track_col() for supplied values")
  fixed <- c(fill = !missing(fill), colour = !missing(colour), alpha = !missing(alpha),
    linewidth = !missing(linewidth), linetype = !missing(linetype))
  dots <- list(...)
  for (key in names(fixed)[fixed]) dots[[key]] <- get(key)
  do.call(geom_track, c(list(mapping = mapping, data = data, geom = "bar", track = track, stat = stat, id = id,
    position = position, bounds = bounds, show.legend = show.legend, inherit.aes = inherit.aes), dots))
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
  if ("fill" %in% names(mapping)) names(mapping)[names(mapping) == "fill"] <- "track_fill"
  if ("colour" %in% names(mapping)) names(mapping)[names(mapping) == "colour"] <- "track_colour"
  keys <- c(fill = ggchord_track_aes("fill", track), colour = ggchord_track_aes("colour", track))
  for (role in names(keys)) names(mapping)[names(mapping) == paste0("track_", role)] <- keys[[role]]
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
  ggchord_check_track_styles(local_styles, c(names(defaults)[fixed], fixed_names),
    id %||% paste0("geom_track_", kind, "()"))
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
    "axis_line_params", "axis_text_params", "axis_title_params", "axis_check_overlap"),
  draw_panel = function(data, panel_params, coord, track_keys, track_kind,
      axis_colour = NULL, axis_size = NULL, na.rm = FALSE, native_draw_params=list(),
      axis_line_params = list(), axis_text_params = list(), axis_title_params = list(), axis_check_overlap = FALSE) {
    for (role in names(track_keys)) data[[role]] <- data[[track_keys[[role]]]]
    polys <- data[data$.component %in% "track_poly", , drop = FALSE]
    lines <- data[data$.component %in% "track_line", , drop = FALSE]
    points <- data[data$.component %in% "track_point", , drop = FALSE]
    labels <- data[data$.component %in% "track_text", , drop = FALSE]
    axis <- data[data$.component %in% "track_axis_line", , drop = FALSE]
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
    sequences <- ggchord_sequence_reference_data(plot)
    for (i in which(pending))
      plot$layers[[i]] <- ggchord_track_native_input(plot$layers[[i]], plot, sequences)
  }
  plot
}

ggchord_prepare_track_layers <- function(plot, sequences = ggchord_sequence_reference_data(plot)) {
  displayed <- plot$ggchord$sequence_layout$order
  subset <- !is.null(displayed) && !all(sequences$accver %in% displayed)
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
    for (name in setdiff(names(inherited), c("x", "y"))) if (!name %in% names(mapping)) mapping[[name]] <- inherited[[name]]
    for (role in c("fill", "colour")) {
      key <- ggchord_track_aes(role, layer$ggchord_params$track_id)
      public <- paste0("track_", role)
      for (name in c(role, public)) if (name %in% names(mapping)) {
        if (!key %in% names(mapping)) mapping[[key]] <- mapping[[name]]
        mapping[[name]] <- NULL
      }
    }
    for (name in intersect(names(layer$aes_params), names(mapping))) mapping[[name]] <- NULL
    layer <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    layer$ggchord_input_mapping <- ggchord_track_stage_mapping(mapping, layer$ggchord_params$track_id)
    layer$inherit.aes <- FALSE
    input <- ggchord_resolve_layer_input(layer, plot$ggchord$data$seq_data, data = source)
    ggchord_track_reference_index(input, sequences, caller)
    layer$ggchord_input_data <- input
    if (subset) layer$ggchord_params$track_sources <- source
    for (role in intersect(layer$ggchord_role_aes, names(mapping))) layer$ggchord_input_mapping[[role]] <- NULL
    layer
  })
  if (subset) {
    # Source mapping and statistics have already used the complete reference.
    # Only the displayed results feed instance expansion, facets and scales.
    for (i in seq_along(plot$layers)) {
      layer <- plot$layers[[i]]
      if (!identical(layer$ggchord_type, "track")) next
      input <- layer$ggchord_input_data
      keep <- input$accver %in% displayed
      layer$ggchord_params$visual_source <- ggchord_visual_rows(layer$ggchord_params$visual_source, which(keep))
      layer$ggchord_input_data <- input[keep, , drop = FALSE]
    }
  }
  numeric_layers <- Filter(function(l) identical(l$ggchord_type, "track") &&
    !identical(l$ggchord_params$track_kind, "tile"), plot$layers)
  ids <- vapply(numeric_layers, function(l) l$ggchord_params$track_id, character(1))
  for (id in unique(ids)) {
    members <- numeric_layers[ids == id]
    measures <- vapply(members, function(l) l$ggchord_params$measure %||% NA_character_, character(1))
    if (length(unique(measures[!is.na(measures)])) > 1L) {
      owners <- vapply(members, function(l) l$ggchord_layer_id, character(1))
      ggchord_stop("Track '", id, "' mixes height measures: ",
        paste(paste0(owners, "=", measures), collapse = ", "), ". Use separate bands for different quantities.")
    }
  }
  plot
}

#' Set a named track's numerical height scale
#'
#' Uses an independent continuous scale for each track ID. The same ID shares
#' its domain across sequences, windows and layers. The scale controls height,
#' not colour. A numerical axis is drawn by default; `guide = "none"` hides it.
#' @param track Track ID from [track_band()], or `"default"`.
#' @param name Axis title, default the track ID.
#' @param limits,breaks,labels,transform,oob Standard continuous-scale controls.
#' @param unit Optional data-unit label retained in track exports.
#' @param facet_domain With `"shared"`, all height layers on this band train
#'   one domain using the panels selected by each layer. With `"panel"`,
#'   each panel trains its own domain and requires a visible track axis.
#' @param guide `"axis"` or `"none"` for the numerical track axis, or [guide_track_axis()]
#'   with a local theme (axis.line/text/title), angle and
#'   check.overlap. Other guide_axis options fail explicitly because this
#'   sequence-local numerical reference is not a Cartesian panel axis. Local
#'   guide styling overrides global ggchord.track.axis.* theme elements.
#' @param ... Additional ggplot2 continuous-scale arguments.
#' @return A ggplot2 scale, scoped to this track.
#' @md
#' @export
#' @examples
#' scale_track_value_continuous("coverage", limits = c(0, 100), unit = "reads")
#' scale_track_value_continuous("depth", facet_domain = "panel", unit = "reads")
#' scale_track_value_continuous("gc", limits = c(0, 1),
#'   guide = guide_track_axis(theme = ggplot2::theme(
#'     axis.text = ggplot2::element_text(size = 9, colour = "#345F58"))))
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

#' Find common numerical limits for explicitly selected tracks
#'
#' Each argument supplies one plot before composition. Repeat a plot to select
#' two different bands. The returned limits are a snapshot of the trained
#' ranges. Use them on each participating scale before assembling a Fusion.
#' @param ... ggchord plots containing the tracks to compare.
#' @param tracks One track ID for all plots, or one ID per plot in order.
#' @return Two numerical limits in the original data units. All participating
#'   height layers must identify their measure. Units, transformation, explicit
#'   limits and missing/out-of-bounds policies must agree.
#' @export
#' @examples
#' ref <- data.frame(accver = "A", length = 20)
#' depth_plot <- function(values) {
#'   ggchord(ref) +
#'     geom_track(aes(x = pos, y = depth), measure = "read_depth",
#'       data = data.frame(accver = "A", pos = c(3, 9), depth = values)) +
#'     coord_collinear()
#' }
#' first <- depth_plot(c(2, 5))
#' second <- depth_plot(c(8, 12))
#' limits <- shared_track_limits(first, second)
#' first + scale_track_value_continuous("default", limits = limits)
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
