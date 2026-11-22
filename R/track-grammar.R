# Native statistical grammar runs on source coordinates, before any projection.
#' Draw a sequence-local data layer with ggplot2 grammar
#'
#' Statistics and data positions are calculated independently within each
#' reference and facet before projection into the physical track band. Supported
#' geoms are point, line, path, step, col, bar, tile, area and interval.
#' Supported statistics are identity,
#' count, bin, summary, smooth, density and ecdf, plus the ggchord interval
#' coverage statistic; unsupported combinations fail
#' in the native layer. Smooth/density/ecdf retain each estimation group's
#' complete eligible source-row set, not new independent observations. These
#' sets describe the group submitted to the native estimator, not per-point
#' influence or a custom model's additional internal filtering decisions.
#' Count and summary provenance lists the eligible input rows at each native
#' group and position, before a data-position adjustment.
#' Source x stays in sequence units. The track value transformation applies
#' to input y before the statistic and to new y values at the after-stat stage.
#' Exported heights use the original measurement units.
#' ECDF requires explicit pad=FALSE because infinite anchors are not source
#' positions. Density's native bounds argument is not forwarded: bounds here
#' controls source-coordinate clipping, not density boundary correction.
#' New estimation combinations use position="identity": smooth with
#' line/path/point, density with line/path/point/area, and ecdf with
#' line/path/point/step. Other combinations are explicitly rejected.
#' Count counts identical positions; bin requires explicit binwidth or breaks.
#' An explicit group mapping or fixed group takes priority over discrete
#' style mappings. For example, `group = 1` pools the source observations.
#' A fixed aesthetic overrides the layer or inherited mapping before statistics
#' and grouping.
#' Point shape categories and their factor order are shared with other point
#' layers through the plot's shape scale. Use `scale_shape_manual()` to choose
#' symbols by category, or set a fixed `shape` on the layer.
#' Missing discrete styles form an implicit group and follow the visual
#' scale's missing-value settings. An explicit group must be non-missing.
#' When exported track geometry combines factor shape categories and fixed
#' numeric point symbols, the combined column uses their character values.
#' Per-layer statistics retain the factor levels of the mapped categories.
#' An explicit sequence_unit in track data must match the reference (bp or aa).
#' When a sequence layer supplies the backbone table, its accessions, lengths
#' and units define the reference for tracks and coverage windows.
#' Bin source membership follows the native statistic, including outer edges
#' and its numerical boundary tolerance, before data positions move the bars.
#' Mapping x/y refers to source position/local numerical value, not screen x/y.
#' @param mapping,data Standard layer input. Supply accver and source positions
#'   through x/pos, or interval start/end. Map y or supply the value column.
#'   Start-stage expressions read the layer's source columns, including the
#'   original accession when `accver` is mapped. A zero-row input must still
#'   provide the columns used by these expressions.
#' @param geom,stat Native geometry and statistic names, or `stat = "coverage"`
#'   for interval coverage. The seven native
#'   statistics also accept their built-in ggplot2 Stat objects, for example
#'   `ggplot2::StatBin`. Bar defaults to count, other geoms to identity.
#'   A ggplot2 Stat object may also compute positions and values for a track.
#'   Its input includes `source_row`, an integer index into this layer's source
#'   table. It can report the contributing indices as a `source_rows` list
#'   column (one integer vector per result). Unreported membership is exported
#'   as NULL with `source_status = "unknown"`. An explicit empty integer vector
#'   denotes a known empty set. Reported members must belong to the same
#'   accession and facet as the result.
#' @param track Character ID declared with [track_band()]. Omit for the shared default band.
#' @param position Native data adjustment such as `"identity"`, `"stack"` or
#'   `ggplot2::position_dodge(width = 10)`. Source x uses sequence units.
#'   Feature placement belongs on feature layers.
#' @param id Layer identifier for visual scales and exported rows. See [auto_id()].
#' @param show.legend Standard ggplot2 layer legend selection. A named logical
#'   vector such as `c(fill = FALSE, colour = TRUE)` controls this layer's
#'   contribution to each guide, including guides on named track scales.
#' @param inherit.aes Whether to inherit the plot mapping.
#' @param ... Native geom/stat parameters and fixed styles. Coverage requires
#'   a positive whole `binwidth` in sequence units or explicit `windows`.
#'   `measure` names the quantity represented by height. Different declared
#'   measures need separate bands. Coverage and mean depth identify their
#'   measures automatically. `bounds` controls whether out-of-range positions
#'   fail or are clipped, and a nonzero `baseline` requires identity placement.
#'   `outline` controls the border drawn by area tracks.
#'   A supplied parameter accepted by both an external Stat and the selected
#'   Geom raises an error. For example, name a custom Stat's search radius
#'   `radius` to use it alongside a column Geom's `width`.
#' @section Interval coverage:
#' Supply one inclusive interval per row with accver, start and end. Endpoint
#' columns can be mapped with `aes(start = ..., end = ...)`. Supply the
#' observation table to the layer. Compound features can declare their actual
#' intervals in `segments`. Each accession or focused window
#' starts at position 1. With binwidth 4, bins are 1-4, 5-8 and so on.
#' The final bin uses its actual length. Cross-bin intervals contribute only
#' their intersection with each bin. Missing endpoints are errors. The default
#' bounds policy rejects out-of-range intervals; clip intersects them with the
#' sequence first. For grouped data, compute a complete grid per observed
#' group within each accession. An accession/sample without observations
#' has no statistical rows. Explicit windows may overlap; see [stat_coverage()].
#'
#' The default height `after_stat(coverage)` is the fraction covered by the
#' union of intervals. `after_stat(mean_depth)` divides the sum of intersection
#' lengths by bin length, counting overlaps repeatedly. Other computed variables
#' are covered (distinct covered units), denominator (actual window length),
#' depth_sum (sum of intersection lengths), n_interval (intersecting input rows),
#' window_start and window_end (inclusive source coordinates). Empty bins have
#' zero values and an empty source_rows list element. Bars span the bin,
#' points/lines/steps use its centre, and tiles fill its full band and
#' interval. Coverage layers use identity position.
#' @return A ggplot2 layer using a shared sequence-local band.
#' @export
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track(aes(x = pos, y = value), data = signal, geom = "line") +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track <- function(mapping = NULL, data = NULL, geom = "point", stat = NULL,
    position = NULL, track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geoms <- c("point", "line", "path", "step", "col", "bar", "area", "tile", "interval")
  if (!is.character(geom) || length(geom) != 1L || !geom %in% geoms)
    ggchord_stop("geom_track(): unsupported geometry")
  dots <- ggchord_colour_dots(list(...), "geom_track()")
  bounds <- ggchord_match_arg(dots$bounds %||% "error", c("error", "clip"))
  baseline <- dots$baseline %||% 0
  measure <- dots$measure
  closed <- dots$closed %||% FALSE
  outline <- dots$outline %||% "upper"
  if (!is.logical(closed) || length(closed) != 1L || is.na(closed))
    ggchord_stop("geom_track(): closed must be TRUE or FALSE")
  if ("closed" %in% names(dots) && !geom %in% c("line", "path"))
    ggchord_stop("geom_track(): only line and path support closing")
  if ("baseline" %in% names(dots) && !geom %in% c("col", "bar", "area"))
    ggchord_stop("geom_track(): baseline requires col, bar or area")
  if ("outline" %in% names(dots) && geom != "area")
    ggchord_stop("geom_track(): outline applies to area")
  if (geom == "interval" && !is.null(measure))
    ggchord_stop("geom_track_interval(): measure requires a numerical height")
  outline <- ggchord_match_arg(outline, c("upper", "both", "full"))
  dots[c("bounds", "baseline", "measure", "closed", "outline")] <- NULL
  if (geom == "interval" && ggchord_identity_stat(stat %||% "identity")) {
    if (!ggchord_is_data_identity(position %||% "identity") ||
        any(c("x", "y", "ymin", "ymax", "height", "width") %in% c(names(mapping), names(dots))))
      ggchord_stop("geom_track_interval(): use interval endpoints and an identity position")
    return(ggchord_track_layer("tile", mapping, data, track,
      "grey60", NA, 1, 0.2, 1, bounds, 0, FALSE, show.legend, inherit.aes,
      c(fill = FALSE, colour = FALSE, alpha = FALSE, linewidth = FALSE, linetype = FALSE),
      c(list(na.rm = na.rm), dots), id))
  }
  interval <- geom == "interval"
  if (interval) geom <- "tile"
  stat <- stat %||% if (geom == "bar") "count" else "identity"
  stats <- c(identity = "StatIdentity", count = "StatCount", bin = "StatBin",
    summary = "StatSummary", smooth = "StatSmooth", density = "StatDensity", ecdf = "StatEcdf")
  custom_stat <- NULL
  if (inherits(stat, "Stat")) {
    builtin <- names(stats)[vapply(stats, function(name)
      identical(stat, getExportedValue("ggplot2", name)), logical(1))]
    if (!length(builtin)) custom_stat <- stat
    stat <- if (length(builtin)) builtin else "custom"
  }
  if (!is.character(stat) || length(stat) != 1L ||
      (!stat %in% c(names(stats), "coverage") && is.null(custom_stat)))
    ggchord_stop("geom_track(): use a ggplot2 Stat object or a supported statistic name")
  position <- position %||% if (geom %in% c("bar", "col", "area")) "stack" else "identity"
  estimates <- list(smooth = c("line", "path", "point"),
    density = c("line", "path", "point", "area"), ecdf = c("line", "path", "point", "step"))
  if (ggchord_local_position(position))
    ggchord_stop("geom_track(): position adjusts source data; use feature placement on a feature layer")
  identity_position <- ggchord_is_data_identity(position)
  if (stat %in% names(estimates) && (!geom %in% estimates[[stat]] || !identity_position))
    ggchord_stop("geom_track(): ", stat, " requires position='identity' and geom in ", paste(estimates[[stat]], collapse = ", "))
  fixed <- character()
  if (length(dots)) {
    # Resolve the source layer's Geom, Stat and Position just as ggplot2 will
    # for its build. Coverage consumes its public windows before that build.
    receiver_stat <- if (stat == "coverage") StatGgchordCoverage else custom_stat %||% stat
    geom_aesthetics <- ggchord_check_geom_params(dots, geom, "geom_track()", position, receiver_stat,
      extra = if (stat == "coverage") c("binwidth", "windows"))
    names(dots) <- ggplot2::standardise_aes_names(names(dots))
    fixed <- intersect(names(dots), geom_aesthetics)
  }
  if (interval && stat != "coverage")
    ggchord_stop("geom_track_interval(): this Stat does not supply source interval edges")
  if (stat == "coverage") {
    if (!geom %in% c("bar", "col", "point", "line", "step", "tile", "path", "area") ||
        (interval && !identity_position))
      ggchord_stop("Coverage requires a compatible track geometry and position")
    width <- dots[["binwidth"]]
    windows <- dots[["windows"]]
    if (is.null(width) == is.null(windows))
      ggchord_stop("Coverage requires exactly one of binwidth or windows")
    if (!is.null(windows) && !is.data.frame(windows))
      ggchord_stop("Coverage windows must be a table of accver, start and end")
    if (!is.null(width) && (!is.numeric(width) || length(width) != 1L ||
        !is.finite(width) || width < 1 || width != floor(width)))
      ggchord_stop("Coverage requires binwidth in positive whole sequence units")
    conflict <- intersect(names(dots), c("width", "breaks", "bins", "boundary", "center", "closed", "pad", "weight", "sequence_length", "window_grid"))
    if (length(conflict)) ggchord_stop("Coverage window specification is incompatible with parameter: ", conflict[1L])
  }
  if (stat == "ecdf" && !identical(dots$pad, FALSE))
    ggchord_stop("geom_track(stat='ecdf'): supply pad=FALSE; infinite ECDF anchors are not source coordinates")
  if (stat == "bin" && is.null(dots$binwidth) && is.null(dots$breaks))
    ggchord_stop("geom_track(stat='bin'): supply binwidth or breaks explicitly")
  if (!is.null(dots$orientation) && !identical(dots$orientation,"x"))
    ggchord_stop("geom_track(): source position uses x; orientation='y' is unsupported")
  if (!is.null(dots$size.unit) && !identical(dots$size.unit,"mm"))
    ggchord_stop("geom_track(): text size currently uses mm")
  if (!is.null(dots$outline.type) && !identical(dots$outline.type,"full"))
    ggchord_stop("geom_track(): only full polygon outlines are supported")
  mapping <- ggchord_normalize_mapping(mapping) %||% aes()
  mapping <- ggchord_track_mapping_aliases(mapping)
  mapping[fixed] <- NULL
  dots$na.rm <- na.rm
  if (!is.numeric(baseline) || length(baseline) != 1L || !is.finite(baseline)) ggchord_stop("Track baseline must be finite")
  if (baseline != 0 && !identity_position)
    ggchord_stop("geom_track(): a nonzero baseline requires position='identity'")
  l <- ggchord_track_layer("line", NULL, data, track, "#86AFB9",
    if (geom %in% c("point", "line", "path", "step", "segment", "linerange", "errorbar", "text")) "#387D88" else NA,
    1, .4, 1, bounds, 0, FALSE, show.legend, FALSE,
    c(fill = FALSE, colour = FALSE, alpha = FALSE, linewidth = FALSE, linetype = FALSE), dots[intersect(c("key_glyph", "name"), names(dots))])
  l$layout <- dots[["layout"]]
  dots$layout <- NULL
  l$ggchord_params$id <- id
  l$ggchord_params$native_geom <- geom
  l$ggchord_params$interval <- interval
  l$ggchord_params$closed <- closed
  l$ggchord_params$outline <- outline
  l$ggchord_params$native_stat <- stat
  l$ggchord_params$custom_stat <- custom_stat
  l$ggchord_params$measure <- ggchord_measure_name(measure)
  l$ggchord_params$native_position <- position
  l$ggchord_params$native_mapping <- mapping
  l$ggchord_params$native_params <- dots
  l$ggchord_params$fixed_aesthetics <- fixed
  l$geom$default_aes$fill <- "grey60"
  l$geom$default_aes$colour <- if (geom %in% c("point", "line", "path", "step")) "black" else NA
  l$geom$default_aes$linewidth <- if (geom %in% c("line", "path", "step")) .5 else
    if (geom == "tile") .2 else .3
  l$geom$default_aes$size <- 1.5
  l$ggchord_params$native_inherit <- inherit.aes
  l$ggchord_params$baseline <- baseline
  l
}

#' Draw a track point layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_point(aes(x = pos, y = value), data = signal, size = 2) +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track_point <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "point", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track line layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_line(aes(x = pos, y = value), data = signal) +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track_line <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "line", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track path layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_path(aes(x = pos, y = value), data = signal) +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track_path <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "path", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track step layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_step(aes(x = pos, y = value), data = signal, direction = "hv") +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track_step <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "step", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track col layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_col(aes(x = pos, y = value), data = signal, width = 40) +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track_col <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "stack", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "col", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track bar layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_bar(aes(x = pos), data = signal, width = 40) +
#'   scale_track_value(name = "Count") +
#'   coord_collinear()
geom_track_bar <- function(mapping = NULL, data = NULL, stat = "count",
    position = "stack", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "bar", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track area layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_area(aes(x = pos, y = value), data = signal, position = "identity") +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track_area <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "stack", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "area", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track tile layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
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
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_tile(aes(x = pos, y = value, fill = value), data = signal, width = 40, height = 0.5) +
#'   scale_track_value(name = "Signal") +
#'   coord_collinear()
geom_track_tile <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "tile", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

#' Draw a track interval layer
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Track observation table, or a function returning it. Numerical measurements use source positions and values.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param track ID of the physical band. The default band is available when first used.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Numerical tracks use `x` in sequence units and `y` in measurement units. They default to the `pos` and `value` columns. An interval fills the band across `start` and `end` without training a numerical height. Track bands set geometry in millimetres.
#' @examples
#' library(ggplot2)
#'
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
#' states <- genes[c("accver", "start", "end", "type")]
#'
#' gene_colours <- c(CDS = "#187B77", tRNA = "#CC7D3A")
#'
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_track_interval(aes(fill = type), data = states) +
#'   scale_fill_manual(values = gene_colours) +
#'   coord_collinear()
geom_track_interval <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", track = "default", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "interval", stat = stat, position = position,
    track = track, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}

ggchord_track_mapping_aliases <- function(mapping) {
  aliases <- c(pos = "x", track_value = "y", track_fill = "fill", track_colour = "colour")
  for (key in intersect(names(aliases), names(mapping))) {
    target <- aliases[[key]]
    if (target %in% names(mapping)) ggchord_stop("geom_track(): use only one of ", key, " and ", target)
    names(mapping)[names(mapping) == key] <- target
  }
  mapping
}

ggchord_track_reference_index <- function(input, sequences, caller) {
  ggchord_require_columns(input, "accver", caller)
  ids <- as.character(input$accver)
  index <- match(ids, sequences$accver)
  bad <- is.na(index) | is.na(ids) | !nzchar(ids)
  if (any(bad)) ggchord_stop(caller, ": unknown reference accession in accver ",
    paste(head(unique(ids[bad]), 5L), collapse = ", "), ggchord_input_rows(bad),
    ". For a focused reference, use its window IDs.")
  index
}

# Evaluate user start stages against the source rows once, including on
# empty input. Pass the values to ggplot2 in columns so they travel with
# their rows through faceting. Stat and scale expressions keep their masks.
ggchord_track_start_aesthetics <- function(mapping, roles, source, data, caller) {
  native_mapping <- mapping
  integer_weight <- FALSE
  for (key in intersect(roles, names(mapping))) {
    start <- ggchord_start_mapping(mapping[[key]])
    if (is.null(start)) next
    value <- ggchord_eval_input_mapping(key, start, source, caller)
    if (!is.null(value) && !NROW(value) %in% c(1L, nrow(source)))
      ggchord_stop(caller, ": mapped ", key, " must be scalar or supply one value per input row")
    # Direct positional mappings to the standard coordinate columns keep
    # their numeric-label conversion. Other expressions see the source type.
    if ((is.character(value) || is.factor(value)) &&
        key %in% c("x", "xmin", "xmax", "xend") &&
        as_label(start) %in% c("pos", "start", "end"))
      value <- as.numeric(as.character(value))
    if (key == "weight") integer_weight <- is.integer(value)
    column <- paste0(".ggchord_start_", key)
    replacement <- NULL
    if (!is.null(value)) {
      data[[column]] <- if (!nrow(data)) head(value, 0L) else value
      replacement <- sym(column)
    }
    expr <- get_expr(mapping[[key]])
    if (is_call(expr, "stage")) {
      expr <- match.call(stage, expr)
      expr["start"] <- list(replacement)
      native_mapping[[key]] <- quo_set_expr(mapping[[key]], expr)
    } else native_mapping[[key]] <- new_quosure(replacement, emptyenv())
  }
  list(data = data, mapping = native_mapping, integer_weight = integer_weight)
}

ggchord_track_native_input <- function(layer, plot, sequences) {
  coverage <- identical(layer$ggchord_params$native_stat, "coverage")
  fallback <- ggchord_layer_fallback(layer, plot$ggchord$data)
  raw <- ggchord_layer_source(layer, fallback)
  tryCatch(ggchord_compute_track_input(layer, plot, raw, sequences, coverage),
    error = function(e) ggchord_abort(ggchord_layer_input_context(layer), parent = e, call = NULL))
}

# Only an unmodified computed column identifies the units of a height.
ggchord_coverage_measure <- function(height) {
  if (is.null(height)) return("coverage")
  measures <- c("coverage", "mean_depth")
  if (is_symbol(height)) {
    legacy <- paste0("..", measures, "..")
    key <- as_string(height)
    if (key %in% legacy) return(measures[match(key, legacy)])
    return(NULL)
  }
  column <- if (is_call(height, "stage")) match.call(stage, height)$after_stat else
    if (is_call(height, c("after_stat", "stat"), n = 1L)) height[[2L]]
  while (is_call(column, "(", n = 1L)) column <- column[[2L]]
  if (is_call(column, c("$", "[["), n = 2L) && identical(column[[2L]], quote(.data))) {
    if (is_call(column, "[[") && !is_string(column[[3L]])) return(NULL)
    key <- column[[3L]]
  } else {
    if (!is_symbol(column)) return(NULL)
    key <- column
  }
  if (is_symbol(key)) key <- as_string(key)
  if (is_string(key) && key %in% measures) key
}

ggchord_compute_track_input <- function(layer, plot, raw, sequences, coverage) {
  p <- layer$ggchord_params
  facet_keys <- ggchord_facet_vars(plot)
  if (!isTRUE(plot$coordinates$ggchord_coord))
    ggchord_stop("geom_track(): use a sequence-path coordinate such as coord_collinear(axis = 'y')")
  if (coverage && is.null(raw)) ggchord_stop("Coverage needs interval data or a plot feature table")
  ggchord_require_columns(raw, character(), "geom_track()")
  value_scale <- plot$scales$get_scales(ggchord_track_aes("value", p$track_id))
  value_transform <- if (!is.null(value_scale) && !isTRUE(p$interval))
    value_scale$get_transformation()
  transformed_height <- !is.null(value_transform) &&
    (value_transform$name != "identity" || isTRUE(value_scale$ggchord_custom_transform))
  if (transformed_height) {
    stacked <- identical(p$native_position, "stack") || inherits(p$native_position, "PositionStack")
    if (stacked && value_transform$name != "reverse")
      ggchord_stop("Track '", p$track_id, "' cannot stack values on a nonlinear scale")
    if (p$native_geom %in% c("bar", "col", "area") &&
        !is.finite(suppressWarnings(value_transform$transform(p$baseline))))
      ggchord_stop("Track '", p$track_id, "' has a baseline incompatible with its scale transform")
  }
  original_source <- raw
  facet_rows <- ggchord_facet_stat_input(raw, plot, layer)
  original_index <- if (is.null(facet_rows)) seq_len(nrow(raw)) else facet_rows$.ggchord_input_index
  raw <- raw[original_index, , drop = FALSE]
  if (!is.null(facet_rows)) raw$.ggchord_stat_panel <- as.integer(facet_rows$PANEL)
  source_data <- raw
  layer$ggchord_params$measure <- p$measure <- ggchord_measure_name(p$measure %||%
    attr(raw, "measure", exact = TRUE))
  layer$ggchord_params$value_unit <- ggchord_measure_name(attr(raw, "value_unit", exact = TRUE) %||%
    attr(raw, "unit", exact = TRUE))
  capture <- p$stage_capture
  caller <- if (coverage) "Coverage intervals" else "Track input"
  ggchord_check_interval_unit(raw, ggchord_sequence_unit(sequences), caller)
  mapping <- p$native_mapping
  fixed <- p$fixed_aesthetics
  if (is.null(fixed)) {
    # Layers saved before fixed styles were recorded still have their original
    # public mapping and named parameters.
    mapping <- ggchord_track_mapping_aliases(mapping)
    if (length(p$native_params)) {
      geom <- ggplot2::layer(geom = p$native_geom, stat = "identity", position = "identity")$geom
      names(p$native_params) <- ggplot2::standardise_aes_names(names(p$native_params))
      fixed <- intersect(names(p$native_params), geom$aesthetics())
      mapping[intersect(names(mapping), fixed)] <- NULL
    }
  }
  if (isTRUE(p$native_inherit)) {
    inherited <- ggchord_track_mapping_aliases(plot$mapping)
    if (isTRUE(p$interval) && any(c("y", "ymin", "ymax") %in% names(inherited))) {
      warning("geom_track_interval(): ignored inherited y aesthetics", call. = FALSE)
      inherited[intersect(c("y", "ymin", "ymax"), names(inherited))] <- NULL
    }
    for (key in setdiff(names(inherited), c(names(mapping), fixed)))
      mapping[[key]] <- inherited[[key]]
  }
  for (key in intersect(c("accver", "group"), names(mapping))) {
    if (ggchord_has_stage(mapping[[key]]))
      ggchord_stop("geom_track(): ", key, " must be available before the source statistic")
  }
  native_roles <- setdiff(names(mapping), c("accver", "group", if (coverage) c("start", "end")))
  reference_data <- raw
  if ("accver" %in% names(mapping)) {
    ids <- ggchord_eval_input_mapping("accver", mapping$accver, raw, caller)
    if (!length(ids) %in% c(1L, nrow(raw)))
      ggchord_stop("geom_track(): mapped accver must be scalar or supply one accession per input row")
    reference_data$accver <- rep_len(ids, nrow(raw))
    if ("accver" %in% facet_keys && !identical(as.character(reference_data$accver), as.character(raw$accver)))
      ggchord_stop("Facet accver conflicts with the mapped reference. Copy the source facet to a sample column.")
    mapping$accver <- NULL
  }
  reference_index <- ggchord_track_reference_index(reference_data, sequences, caller)
  if (coverage) {
    height <- get_expr(mapping$y)
    known <- ggchord_coverage_measure(height)
    if (isTRUE(p$interval) && (!is.null(height) || "height" %in% c(names(mapping), names(p$native_params))))
      ggchord_stop("Coverage tiles fill the whole band; map fill to after_stat(coverage) or after_stat(mean_depth)")
    if (!is.null(known) && !is.null(p$measure) && !identical(known, p$measure))
      ggchord_stop("The mapped height represents ", known, "; the declared measure is ", p$measure)
    layer$ggchord_params$measure <- p$measure %||% known
    forbidden <- intersect(names(mapping), c("x", "xmin", "xmax", "ymin", "ymax", "width", "weight"))
    if (length(forbidden) || ("y" %in% names(mapping) &&
        (!is.null(ggchord_start_mapping(mapping$y)) || ggchord_has_stage(mapping$y, "after_scale"))))
      ggchord_stop("Coverage maps start and end; map y with after_stat() and set window size with binwidth")
    raw <- ggchord_coverage_input(raw, mapping, sequences, p$bounds, reference_index)
    raw$.track_reference <- reference_index
    mapping$start <- mapping$end <- NULL
    mapping$xmin <- quo(.coverage_start)
    mapping$xmax <- quo(.coverage_end)
  }
  for(field in if (!coverage) intersect(c("pos","start","end"),names(raw))) {
    if(is.factor(raw[[field]]) || is.character(raw[[field]])) {
      value <- suppressWarnings(as.numeric(as.character(raw[[field]])))
      if(any(is.na(value) & !is.na(raw[[field]])))ggchord_stop("geom_track(): ",field," must contain numeric coordinate labels")
      raw[[field]]<-value
    }
  }
  interval_x <- FALSE
  if (!coverage && !"x" %in% names(mapping)) {
    if ("pos" %in% names(raw)) mapping$x <- quo(pos) else
      ggchord_stop("geom_track(): map x or supply the pos column")
  }
  if (!coverage && p$native_stat != "count" && !"y" %in% names(mapping) &&
      "value" %in% names(raw)) mapping$y <- quo(value)
  if (p$native_stat == "count" && "y" %in% names(mapping) &&
      !is.null(ggchord_start_mapping(mapping$y))) ggchord_stop(
    "geom_track_bar(): count does not accept supplied heights; use geom_track_col() or explicit stat='identity'")
  layer$ggchord_params$value_label <- if (!is.null(mapping$y)) as_label(mapping$y) else
    if (p$native_stat == "count") "count"
  late <- names(mapping)[vapply(mapping, ggchord_has_stage, logical(1), stages = "after_scale")]
  late_mapping <- mapping[late]
  for (key in late) {
    expr <- get_expr(mapping[[key]])
    if (is_call(expr, "stage")) {
      expr <- match.call(stage, expr)
      expr$after_scale <- NULL
      if (!is.null(expr$after_stat) || !is.null(expr$start)) {
        mapping[[key]] <- quo_set_expr(mapping[[key]], expr)
        next
      }
    }
    mapping[[key]] <- NULL
  }
  params <- p$native_params
  fixed_group <- params[["group"]]
  implicit_group <- !"group" %in% names(mapping) && is.null(fixed_group)
  if (!implicit_group) {
    source_group <- fixed_group %||% ggchord_eval_input_mapping("group", mapping$group, source_data, caller)
    if (!length(source_group) %in% c(1L, nrow(raw)))
      ggchord_stop("geom_track(): group must be scalar or supply one value per input row")
    if (anyNA(source_group)) ggchord_stop("geom_track(): group values must be non-missing")
    raw$.track_group <- rep_len(source_group, nrow(raw))
    mapping$group <- quo(.track_group)
    params$group <- NULL
  }
  if (!coverage && !is.null(raw[[".track_width"]]) && is.null(params$width) &&
      is.null(mapping$width) && p$native_geom %in% c("bar", "col", "tile")) {
    mapping$width <- quo(.track_width)
  }
  # Partition source statistics by accession and facet values in ONE native
  # layer. ggplot2 evaluates after_stat() after all panels have been combined.
  domain_keys <- c("accver", if (!is.null(facet_rows)) ".ggchord_stat_panel")
  reference_data$.ggchord_stat_panel <- raw$.ggchord_stat_panel
  domain_id <- do.call(interaction, c(unname(lapply(reference_data[domain_keys], function(x)
    match(x, unique(x)))), list(drop = TRUE, lex.order = TRUE)))
  raw$.track_domain <- as.integer(domain_id)
  input_domain <- raw$.track_domain
  raw$.track_input_index <- seq_len(nrow(raw))
  domains <- reference_data[domain_keys]
  domains$.track_domain <- raw$.track_domain
  domains <- domains[!duplicated(domains$.track_domain), , drop = FALSE]
  stat <- p$custom_stat %||% p$native_stat
  if (!is.null(p$custom_stat) && any(c("source_row", "source_rows") %in% names(mapping)))
    ggchord_stop("Custom track statistics reserve source_row and source_rows for provenance")
  if (coverage) {
    stat <- StatGgchordCoverage
    if (isTRUE(p$interval))
      stat <- ggproto(NULL, stat, default_aes = aes(y = after_stat(0.5)))
    params$window_grid <- ggchord_coverage_grid(sequences, params$binwidth, params$windows,
      seq_len(nrow(sequences)))
    params$binwidth <- params$windows <- NULL
    mapping$.ggchord_reference <- quo(.track_reference)
  }
  source <- source_data
  start_aes <- ggchord_track_start_aesthetics(mapping, native_roles, source, raw, caller)
  if (!is.null(capture) && "mapped" %in% capture$requested) {
    mapped <- source
    mapped$accver <- reference_data$accver
    for (key in native_roles) {
      value <- start_aes$data[[paste0(".ggchord_start_", key)]]
      if (!is.null(value)) mapped[[key]] <- value
    }
    if (coverage) {
      mapped$start <- raw$.coverage_start
      mapped$end <- raw$.coverage_end
    }
    mapped$source_row <- original_index
    capture$mapped <- mapped
  }
  if (coverage && nrow(raw) && !is.null(raw$segments)) {
    expanded <- lapply(seq_len(nrow(raw)), function(i) {
      parts <- raw$segments[[i]]
      if (is.null(parts) || !nrow(parts)) return(start_aes$data[i, , drop = FALSE])
      ggchord_require_columns(parts, c("start", "end"), "Coverage segments")
      if ("segment_type" %in% names(parts))
        parts <- parts[is.na(parts$segment_type) | parts$segment_type != "gap", , drop = FALSE]
      start <- pmax(parts$start, raw$.coverage_start[i]); end <- pmin(parts$end, raw$.coverage_end[i])
      keep <- start <= end
      result <- start_aes$data[rep(i, sum(keep)), , drop = FALSE]
      result$.coverage_start <- start[keep]; result$.coverage_end <- end[keep]
      result
    })
    start_aes$data <- ggchord_rbind_fill(expanded)
  }
  if (coverage && implicit_group) {
    style_columns <- paste0(".ggchord_start_", c("fill", "colour", "linetype", "shape"))
    style_columns <- intersect(style_columns, names(start_aes$data))
    grouped <- any(vapply(start_aes$data[style_columns], function(x) !is.numeric(x), logical(1)))
    if (!grouped) {
      reference_members <- ggchord_facet_stat_input(sequences, plot, layer)
      if (is.null(reference_members)) {
        missing <- data.frame(accver = sequences$accver)
      } else {
        missing <- data.frame(accver = sequences$accver[reference_members$.ggchord_input_index],
          .ggchord_stat_panel = as.integer(reference_members$PANEL))
      }
      missing <- missing[!duplicated(missing), , drop = FALSE]
      domain_key <- function(x) paste(x$accver, x$.ggchord_stat_panel %||% 1L, sep = "\r")
      missing <- missing[!domain_key(missing) %in% domain_key(domains), , drop = FALSE]
      if (nrow(missing)) {
        missing$.track_domain <- nrow(domains) + seq_len(nrow(missing))
        domains <- rbind(domains, missing)
        # A zero-width sentinel before the first residue trains the missing
        # reference domain, without contributing an interval to any window.
        zeros <- start_aes$data[rep(NA_integer_, nrow(missing)), , drop = FALSE]
        zeros$.track_domain <- missing$.track_domain
        zeros$.track_reference <- match(missing$accver, sequences$accver)
        zeros$.coverage_start <- zeros$.coverage_end <- rep(0, nrow(missing))
        zeros$.track_input_index <- rep(0L, nrow(missing))
        start_aes$data <- rbind(start_aes$data, zeros)
      }
    }
  }
  if (nrow(start_aes$data)) {
    # Identity scales retain source aesthetic values; final role scales train
    # only after projection. No device geometry enters these statistics.
    native <- do.call(getExportedValue("ggplot2", paste0("geom_", p$native_geom)),
      c(list(mapping = start_aes$mapping, data = start_aes$data, stat = stat,
        position = p$native_position, inherit.aes = FALSE), params))
    if (p$native_stat != "identity") {
      # Preserve Stat coordinates before a data Position moves the display
      # centre/edges. Source membership must not follow a dodge displacement.
      source_geom <- native$geom
      native$geom <- ggproto(NULL,source_geom)
      native$geom$default_aes$.ggchord_stat_x <- NA_real_
      native$mapping$.ggchord_stat_x <- aes(value = after_stat(x))$value
    }
    lineage_key <- if (is.null(p$custom_stat)) ".ggchord_source_index" else "source_row"
    native$mapping[[lineage_key]] <- quo(.track_input_index)
    source_stat <- native$stat
    native$stat <- ggproto(NULL, source_stat,
      required_aes = unique(c(source_stat$required_aes, lineage_key)),
      dropped_aes = unique(c(source_stat$dropped_aes, lineage_key)))
    if (p$native_stat == "bin") native$stat <- ggchord_track_bin_stat(native$stat)
    if (p$native_stat %in% c("smooth", "density", "ecdf"))
      native$stat <- ggchord_track_group_stat(native$stat)
    if (p$native_stat %in% c("count", "summary"))
      native$stat <- ggchord_track_x_stat(native$stat)
    failures <- capture %||% new.env(parent = emptyenv())
    native <- ggchord_track_source_layer(native, failures, implicit_group, interval_x, p$bounds)
    # A normal identity scale converts factors to character in map().
    # Keep their level order for the final plot's shared point-shape scale.
    shape_scale <- ggplot2::scale_discrete_identity(aesthetics = "shape")
    shape_scale$map <- identity
    q <- ggplot() + native + scale_fill_identity() +
      ggplot2::scale_colour_identity() + ggplot2::scale_alpha_identity() +
      ggplot2::scale_size_identity() + ggplot2::scale_linewidth_identity() + shape_scale +
      ggplot2::scale_continuous_identity(aesthetics=".ggchord_stat_x")
    if (!is.null(value_transform)) {
      limits <- value_scale$ggchord_domain_policy$limits
      # Data-dependent limits use the final shared or facet-specific domain.
      q <- q + scale_y_continuous(transform = value_transform,
        limits = if (is.numeric(limits)) limits, oob = value_scale$oob)
    }
    if (nrow(domains) > 1L)
      q <- q + ggplot2::facet_wrap(ggplot2::vars(.track_domain), scales = "free")
    built <- ggplot_build(q)
    if (is.null(capture) && length(failures$statistic_errors))
      ggchord_stop("Statistic failed: ", paste(failures$statistic_errors, collapse = "\n"))
    computed <- built$data[[1L]]
    if (transformed_height) for (field in intersect(c("y", "ymin", "ymax"), names(computed)))
      computed[[field]] <- value_transform$inverse(computed[[field]])
    panel <- built$layout$layout
    domain <- if (nrow(domains) > 1L) panel$.track_domain[match(computed$PANEL, panel$PANEL)] else
      rep(domains$.track_domain[1L], nrow(computed))
    at <- match(domain, domains$.track_domain)
    computed$accver <- domains$accver[at]
    computed$PANEL <- NULL
    computed$pos <- computed$x %||% ((computed$xmin + computed$xmax) / 2)
    computed$track_value <- computed$y %||% computed$ymax
    computed$start <- computed$xmin %||% computed$pos
    computed$end <- computed$xmax %||% computed$pos
    computed$track_ymin <- computed$ymin %||% computed$track_value
    computed$track_ymax <- computed$ymax %||% computed$track_value
    if (p$baseline != 0) {
      computed$track_ymin <- pmin(p$baseline, computed$track_value)
      computed$track_ymax <- pmax(p$baseline, computed$track_value)
    }
    # Positions move display coordinates; memberships index the full layer input.
    members <- switch(p$native_stat,
      identity = as.list(computed$.ggchord_source_index),
      bin = computed$.ggchord_bin_rows,
      custom = ggchord_track_custom_members(computed$source_rows, domain, input_domain),
      computed$.ggchord_group_rows)
    members <- lapply(members, function(rows) {
      if (is.null(rows)) return(NULL)
      original_index[as.integer(rows)]
    })
    computed$source_rows <- I(members)
    if (coverage) computed$source_members <- I(members)
    if (!is.null(p$custom_stat)) {
      computed$source_row <- NULL
      computed$source_status <- ifelse(vapply(members, is.null, logical(1)), "unknown", "reported")
    }
    if (coverage && nrow(computed)) {
      computed$start <- computed$window_start - if (isTRUE(p$interval)) 0 else 1
      computed$end <- computed$window_end
      computed$width <- computed$denominator
    }
    if (p$native_stat %in% c("count", "summary", "coverage")) computed$.ggchord_source_index <- NULL
    computed$.ggchord_bin_rows <- NULL
    computed$.ggchord_group_rows <- NULL
    output <- computed
  } else output <- NULL
  if (is.null(output) || !nrow(output)) {
    if (is.null(output)) output <- domains[FALSE, "accver", drop = FALSE]
    # ggplot2 skips a Stat with no input. Supply its documented variables
    # without invoking an estimator or guessing extension-specific columns.
    numeric_fields <- switch(p$native_stat,
      count = c("x", "y", "count", "prop", "width"),
      bin = c("x", "y", "xmin", "xmax", "width", "count", "density", "ncount", "ndensity"),
      density = c("x", "y", "density", "scaled", "ndensity", "count", "wdensity"),
      ecdf = c("x", "y", "ecdf"),
      smooth = c("x", "y", if (!identical(params$se, FALSE)) c("ymin", "ymax", "se")),
      summary = c("x", "y", if (is.null(params$fun.data)) c("ymin", "ymax")),
      identity = intersect(names(mapping), c("x", "y", "xmin", "xmax", "ymin", "ymax", "xend", "yend")),
      coverage = c("x", "y", "xmin", "xmax"))
    if (p$native_stat == "count" && !"count" %in% names(output) && start_aes$integer_weight)
      output$count <- integer()
    for (key in setdiff(numeric_fields, names(output))) output[[key]] <- numeric()
    if (p$native_stat == "density" && !"n" %in% names(output)) output$n <- integer()
    output$pos <- output$start <- output$end <- numeric()
    output$track_value <- output$track_ymin <- output$track_ymax <- numeric()
    output$group <- integer()
    if (coverage) {
      grid <- params$window_grid[[1L]][FALSE, , drop = FALSE]
      for (key in names(grid)) output[[key]] <- grid[[key]]
      output$width <- output$covered <- output$depth_sum <- numeric()
      output$coverage <- output$mean_depth <- numeric()
      output$n_interval <- integer()
    }
    output$source_rows <- I(list())
    if (!is.null(p$custom_stat)) output$source_status <- character()
    mapping <- late_mapping <- aes()
  }
  if (!is.null(facet_rows)) {
    panel <- if (nrow(output)) domains$.ggchord_stat_panel[at] else integer()
    output$.ggchord_stat_panel <- panel
    panels <- plot$ggchord$source_facet_setup$panels
    for (key in setdiff(facet_keys, "accver")) output[[key]] <- panels[[key]][panel]
  }
  if (!is.null(capture) && "statistic" %in% capture$requested) {
    statistic <- output
    # Coverage positions here are closed source windows. The preceding
    # boundary retained by the renderer belongs only to the displayed bar.
    if (coverage) statistic$start <- statistic$window_start
    capture$statistic <- statistic
  }
  visual <- aes()
  for (key in intersect(c("fill", "colour", "alpha", "linewidth", "linetype", "size", "shape", "label", "angle", "hjust", "vjust", "family", "fontface", "lineheight", "stroke"), names(mapping))) {
    column <- paste0(".track_", key)
    output[[column]] <- output[[key]]
    target <- key
    visual[[target]] <- new_quosure(sym(column), emptyenv())
  }
  for (key in names(late_mapping)) {
    target <- key
    q <- ggchord_track_stage_mapping(late_mapping[key], p$track_id)[[key]]
    expr <- get_expr(q)
    if (is_call(expr, "stage")) {
      expr$after_stat <- NULL
      expr$start <- if (key %in% names(mapping)) sym(paste0(".track_", key)) else NULL
      q <- quo_set_expr(q, expr)
    }
    visual[[target]] <- q
  }
  # Fixed styles have already been interpreted by the native geom.
  for (key in intersect(c("fill", "colour", "alpha", "linewidth", "linetype", "size", "shape", "angle", "hjust", "vjust", "family", "fontface", "lineheight", "stroke"), names(params))) {
    target <- key
    layer$aes_params[[target]] <- params[[key]]
  }
  layer$ggchord_input_mapping <- visual
  mapping_label <- function(q) {
    expr <- get_expr(q)
    if (is_call(expr, "stage")) {
      stages <- match.call(stage, expr)
      expr <- stages$after_stat %||% stages$start
    } else if (is_call(expr, c("after_stat", "stat"))) expr <- expr[[2L]]
    as_label(expr)
  }
  layer$ggchord_params$track_mapping_labels <- setNames(
    vapply(mapping, mapping_label, character(1)),
    names(mapping))
  layer$ggchord_input_data <- output
  layer$ggchord_params$track_kind <- if (coverage && isTRUE(p$interval)) "tile" else
    if (p$native_geom %in% c("bar", "col")) "bar" else "line"
  layer$ggchord_params$source_edges <- TRUE
  layer$ggchord_params$track_sources <- original_source
  layer$ggchord_params$native_prepared <- TRUE
  layer$geom_params$native_draw_params <- params[intersect(names(params),c("arrow","lineend","linejoin","linemitre","parse","check_overlap","rule"))]
  layer$geom_params$track_kind <- p$native_geom
  layer
}

# Inspect the native mapping result before scales and statistics. Catch Stat
# errors at their source, including failures ggplot2 turns into panel warnings.
ggchord_track_source_layer <- function(layer, capture, implicit_group, interval_x, bounds) {
  record <- function(error) {
    capture$statistic_errors <- union(capture$statistic_errors, conditionMessage(error))
  }
  panel_error <- NULL
  stat <- layer$stat
  layer$stat <- ggproto(NULL, stat,
    parameters = function(extra = FALSE) stat$parameters(extra),
    compute_panel = function(self, data, scales, ...) {
      panel_error <<- NULL
      withCallingHandlers(ggproto_parent(stat, self)$compute_panel(data, scales, ...),
        error = function(e) panel_error <<- e)
    })
  # ggplot2 attaches the escaping panel error as its warning's parent. An
  # extension may instead recover in compute_layer() and return valid data.
  ggproto(NULL, layer,
    compute_aesthetics = function(self, data, plot) {
      mapped <- ggproto_parent(layer, self)$compute_aesthetics(data, plot)
      x <- mapped[["x"]]
      if (!is.null(x)) {
        if (!is.numeric(x)) ggchord_stop("geom_track(): x must contain numeric source positions; convert factor labels explicitly")
        if (any(is.infinite(x))) ggchord_stop("geom_track(): source positions must be finite")
        if (!interval_x && identical(bounds, "error") && any(x < 1, na.rm = TRUE))
          ggchord_stop("geom_track(): source positions are one-based; x/pos must be at least 1")
      }
      if (implicit_group) {
        styles <- mapped[intersect(c("fill", "colour", "linetype", "shape"), names(mapped))]
        discrete <- Filter(function(value) !is.numeric(value), styles)
        # Missing style values form an implicit group. Leave the aesthetic
        # itself untouched so the final scale controls its missing category.
        codes <- lapply(discrete, function(value) as.integer(addNA(value, ifany = TRUE)))
        mapped$group <- if (length(codes)) as.integer(do.call(interaction,
          c(unname(codes), list(drop = TRUE, lex.order = TRUE)))) else rep(1L, nrow(mapped))
      }
      mapped
    },
    compute_statistic = function(self, data, layout) {
      tryCatch(withCallingHandlers(
        ggproto_parent(layer, self)$compute_statistic(data, layout),
        warning = function(w) {
          if (!is.null(panel_error) && identical(w$parent, panel_error)) record(panel_error)
        }),
        error = function(e) {
          record(e)
          data.frame()
        })
    })
}

ggchord_track_custom_members <- function(members, domain, input_domain) {
  if (is.null(members)) return(vector("list", length(domain)))
  if (!is.list(members) || length(members) != length(domain))
    ggchord_stop("Custom track Stat source_rows must be a list with one entry per output row")
  lapply(seq_along(members), function(i) {
    rows <- members[[i]]
    if (is.null(rows)) return(NULL)
    if (!is.numeric(rows) || any(!is.finite(rows) | rows != floor(rows) |
        rows < 1 | rows > length(input_domain)) || anyDuplicated(rows) ||
        any(input_domain[rows] != domain[i]))
      ggchord_stop("Custom track Stat reported invalid source_rows at output row ", i,
        "; use the source_row indices supplied for this reference and facet")
    as.integer(rows)
  })
}

# Count and summary aggregate the eligible rows at a particular x and group.
# Capture native input after Stat's missing-value policy, before dodge or nudge.
ggchord_track_x_stat <- function(stat) {
  ggproto(NULL, stat,
    dropped_aes = unique(c(stat$dropped_aes, ".ggchord_source_index")),
    parameters = function(self, extra = FALSE) stat$parameters(extra = extra),
    compute_panel = function(data, scales, ...) {
      inputs <- data
      inputs$.ggchord_source_index <- NULL
      out <- stat$compute_panel(inputs, scales, ...)
      positions <- unique(data$x)
      key <- function(group, x) paste(group, match(x, positions), sep = "/")
      # Native weighted counts ignore missing weights inside rowsum().
      if (inherits(stat, "StatCount") && !is.null(data$weight))
        data <- data[!is.na(data$weight), , drop = FALSE]
      members <- split(data$.ggchord_source_index, key(data$group, data$x))
      out$.ggchord_group_rows <- I(unname(members[key(out$group, out$x)]))
      out
    })
}

ggchord_track_group_stat <- function(stat) {
  # Preserve native parameter discovery while wrapping the actual eligible
  # compute_group input, after native missing-value removal and before Position.
  # Model/density estimates depend on this group, not just the output x value.
  ggproto(NULL, stat,
    dropped_aes = unique(c(stat$dropped_aes, ".ggchord_source_index")),
    parameters = function(self, extra = FALSE) stat$parameters(extra = extra),
    compute_group = function(data, scales, ...) {
      result <- stat$compute_group(data, scales, ...)
      result$.ggchord_group_rows <- I(rep(list(data$.ggchord_source_index), nrow(result)))
      result
    })
}

ggchord_track_bin_stat <- function(stat) {
  # Keep the dependency on native bin internals in this adapter. Their fuzzy
  # edges/outer inclusion must be shared with StatBin, not approximated here.
  # Verified against ggplot2 >= 4.0.0's R/bin.R and R/stat-bin.R:
  # https://github.com/tidyverse/ggplot2/tree/v4.0.0/R
  compute_bins <- utils::getFromNamespace("compute_bins", "ggplot2")
  bin_cut <- utils::getFromNamespace("bin_cut", "ggplot2")
  ggproto(NULL, stat,
    dropped_aes = unique(c(stat$dropped_aes, ".ggchord_source_index")),
    compute_group = function(data, scales, binwidth = NULL, bins = NULL,
        center = NULL, boundary = NULL, closed = c("right", "left"),
        pad = FALSE, breaks = NULL, flipped_aes = FALSE, drop = "none") {
      resolved <- compute_bins(data$x, scales$x, breaks = breaks,
        binwidth = binwidth, bins = bins, center = center, boundary = boundary,
        closed = closed)
      membership <- as.integer(bin_cut(data$x, resolved))
      # Resolve user callbacks only once. Track's native x scale is identity.
      result <- stat$compute_group(data, scales, breaks = resolved$breaks,
        closed = closed, pad = pad, flipped_aes = flipped_aes, drop = drop)
      centers <- (head(resolved$breaks, -1L) + tail(resolved$breaks, -1L)) / 2
      at <- match(result$x, centers)
      result$.ggchord_bin_rows <- I(lapply(at, function(i) {
        if (is.na(i)) integer() else data$.ggchord_source_index[which(membership == i)]
      }))
      result
    })
}

ggchord_track_native_geometry <- function(d, p, layout, sc, at_value, draw_rows) {
  pieces <- list(); group <- 0L
  map_value <- function(x) sc$map(suppressWarnings(sc$transform(x)))
  add <- function(pos, value, rows, component) {
    if (!length(pos)) return()
    mapped <- map_value(value)
    if (any(!is.finite(pos)) || any(!is.finite(mapped))) return()
    if (any(mapped < 0 | mapped > 1)) ggchord_stop("Track values extend outside their physical band")
    acc <- d$accver[rows[1L]]
    len <- layout$sequence_reference$lens[[acc]]
    if (any(pos < 0 | pos > len)) {
      if (p$bounds == "error") ggchord_stop("Track source geometry outside sequence boundaries; set bounds='clip' or adjust width")
      pos <- pmax(0, pmin(pos, len))
    }
    xy <- at_value(pos, acc, mapped)
    group <<- group + 1L
    pieces[[length(pieces) + 1L]] <<- data.frame(x = xy[, 1], y = xy[, 2],
      group = group, source_row = rep_len(rows, length(pos)), accver = acc,
      track_id = p$track_id,
      label = if (p$native_geom == "text") as.character(d$label[rep_len(rows, length(pos))]) else NA_character_,
      .component = component)
  }
  path <- function(x, y, rows) {
    if (length(x) < 2L) return(list(x = x, y = y, rows = rows))
    xx <- yy <- numeric(); rr <- integer()
    for (i in seq_len(length(x) - 1L)) {
      n <- if (isTRUE(layout$collinear)) 2L else 16L
      xx <- c(xx, seq(x[i], x[i+1L], length.out = n))
      yy <- c(yy, seq(y[i], y[i+1L], length.out = n))
      rr <- c(rr, rep(rows[i], n-1L), rows[i+1L])
    }
    list(x = xx, y = yy, rows = rr)
  }
  geom <- p$native_geom
  if (geom %in% c("line", "path", "step", "area", "ribbon")) {
    series <- split(seq_len(nrow(d)), interaction(d$accver, d$group, drop = TRUE))
    for (rows in series) {
      if (!identical(geom, "path")) rows <- rows[order(d$pos[rows])]
      valid <- is.finite(d$pos[rows]) & is.finite(d$track_ymin[rows]) & is.finite(d$track_ymax[rows])
      runs <- split(rows[valid], cumsum(!valid)[valid])
      for (run in runs) {
        if (length(run) < 2L) next
        x <- d$pos[run]
        upper <- if (geom %in% c("area", "ribbon")) d$track_ymax[run] else d$track_value[run]
        lower <- d$track_ymin[run]
        if (geom == "step") {
          direction <- p$native_params$direction %||% "hv"
          if (direction == "mid") {
            mid <- (head(x,-1)+tail(x,-1))/2
            x <- c(x[1], rep(mid, each=2), tail(x,1))
            upper <- rep(upper, each=2); run <- rep(run, each=2)
          } else {
            indices <- if (direction == "vh") c(1, rep(seq.int(2,length(x)), each=2)) else c(rep(seq_len(length(x)-1),each=2),length(x))
            x <- if (direction == "vh") c(rep(head(x,-1),each=2),tail(x,1)) else c(x[1],rep(tail(x,-1),each=2))
            upper <- upper[indices]; run <- run[indices]
          }
        }
        top <- path(x, upper, run)
        if (geom %in% c("area", "ribbon")) {
          bottom <- path(rev(d$pos[run]), rev(lower), rev(run))
          partial <- geom == "area" && !identical(p$outline, "full")
          add(c(top$x,bottom$x), c(top$y,bottom$y), c(top$rows,bottom$rows),
            if (partial) "track_area_fill" else "track_poly")
          if (partial) add(top$x, top$y, top$rows, "track_line")
          if (partial && identical(p$outline, "both"))
            add(bottom$x, bottom$y, bottom$rows, "track_line")
        } else {
          if (isTRUE(p$closed) && length(run) == length(rows) && all(valid)) {
            len <- layout$sequence_reference$lens[[d$accver[run[1L]]]]
            seam <- seq(tail(x, 1L), x[1L] + len, length.out = 24L)
            top$x <- c(top$x, seam %% len)
            top$y <- c(top$y, seq(tail(upper, 1L), upper[1L], length.out = 24L))
            top$rows <- c(top$rows, rep(tail(run, 1L), 24L))
          }
          add(top$x, top$y, top$rows, "track_line")
        }
      }
    }
  } else for (i in draw_rows) {
    if (geom %in% c("point", "text")) add(d$pos[i], d$track_value[i], i,
      if (geom == "point") "track_point" else "track_text") else if (geom %in% c("bar", "col", "rect", "tile")) {
      n <- if (isTRUE(layout$collinear)) 2L else 32L
      pos <- seq(d$start[i], d$end[i], length.out = n)
      add(c(pos, rev(pos)), c(rep(d$track_ymin[i],n),rep(d$track_ymax[i],n)), i, "track_poly")
    } else if (geom == "segment") {
      v <- path(c(d$pos[i], d$xend[i]), c(d$track_value[i], d$yend[i]), c(i,i))
      add(v$x,v$y,v$rows,"track_line")
    } else {
      add(rep(d$pos[i],2), c(d$track_ymin[i],d$track_ymax[i]), i, "track_line")
      if (geom == "errorbar") for (y in c(d$track_ymin[i],d$track_ymax[i])) {
        v <- path(c(d$start[i],d$end[i]),c(y,y),c(i,i)); add(v$x,v$y,v$rows,"track_line")
      }
    }
  }
  if (length(pieces)) ggchord_rbind_fill(pieces) else ggchord_empty_track_geometry()
}
