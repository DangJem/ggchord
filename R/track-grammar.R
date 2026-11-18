# Native statistical grammar runs on source coordinates, before any projection.
#' Draw a sequence-local data layer with ggplot2 grammar
#'
#' Statistics and data positions are calculated independently within each
#' accession/window before projection into the physical track band. Supported
#' geoms are point, line, path, step, col, bar, tile, rect, area, ribbon,
#' segment, linerange, errorbar and text. Supported statistics are identity,
#' count, bin, summary, smooth, density and ecdf, plus the ggchord interval
#' coverage statistic; unsupported combinations fail
#' in the native layer. Smooth/density/ecdf retain each estimation group's
#' complete eligible source-row set, not new independent observations. These
#' sets describe the group submitted to the native estimator, not per-point
#' influence or a custom model's additional internal filtering decisions.
#' Count and summary provenance lists the eligible input rows at each native
#' group and position, before a data-position adjustment.
#' Statistics operate on raw source x/y before Track numerical-scale transforms.
#' ECDF requires explicit pad=FALSE because infinite anchors are not source
#' positions. Density's native bounds argument is not forwarded: bounds here
#' controls source-coordinate clipping, not density boundary correction.
#' New estimation combinations use position="identity": smooth with
#' line/path/point/ribbon, density with line/path/point/area, and ecdf with
#' line/path/point/step. Other combinations are explicitly rejected.
#' Count counts identical positions; bin requires explicit binwidth or breaks.
#' An explicit group mapping or fixed group takes priority over discrete
#' style mappings. For example, `group = 1` pools the source observations.
#' A fixed aesthetic overrides the inherited mapping before statistics and
#' grouping. Mapping and fixing the same aesthetic in one layer is an error.
#' An explicit sequence_unit in track data must match the reference (bp or aa).
#' When a sequence layer supplies the backbone table, its accessions, lengths
#' and units define the reference for tracks and coverage windows.
#' Bin source membership follows the native statistic, including outer edges
#' and its numerical boundary tolerance, before data positions move the bars.
#' Mapping x/y refers to source position/local numerical value, not screen x/y.
#' Existing track_value/track_fill/track_colour roles remain accepted.
#' @param mapping,data Standard layer input. Supply accver and source positions
#'   through x/pos, or interval start/end. Map y/track_value for identity values.
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
#' @param id Optional unique layer ID for exported statistics and source rows.
#' @param measure Optional name of the quantity represented by height, e.g.
#'   `"read_depth"`. Layers with different declared measures need separate
#'   bands. Unmodified coverage and mean_depth outputs identify themselves,
#'   including direct references such as `after_stat(.data[["mean_depth"]])`.
#' @param bounds Coordinate error/clip policy.
#' @param baseline Numeric identity-column baseline. Nonzero baselines require
#'   position='identity' so stacking does not change their meaning.
#' @param fill,colour,alpha,linewidth,linetype Fixed styles for the count shortcut.
#' @param show.legend,inherit.aes Standard layer controls.
#' @param ... Native geom/stat parameters and fixed styles. Coverage requires
#'   a positive whole `binwidth` in sequence units or explicit `windows`.
#'   A supplied parameter accepted by both an external Stat and the selected
#'   Geom raises an error. For example, name a custom Stat's search radius
#'   `radius` to use it alongside a column Geom's `width`.
#' @section Interval coverage:
#' Supply one inclusive interval per row with accver, start and end. Endpoint
#' columns can be mapped with `aes(start = ..., end = ...)`. Without layer data
#' the statistic uses the plot's feature table. Split compound or origin-crossing
#' intervals into separate records first. Each accession or focused window
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
#' @examples
#' # Synthetic observations illustrate counting, not biological measurements.
#' events <- data.frame(accver = "A", pos = c(100, 100, 400))
#' ggchord(data.frame(accver = "A", length = 1000)) +
#'   geom_seq() +
#'   geom_track_bar(data = events, stat = "count") +
#'   coord_collinear()
#' intervals <- data.frame(accver = "A", start = c(3, 4), end = c(6, 5))
#' ggchord(data.frame(accver = "A", length = 10)) +
#'   geom_track_bar(aes(y = after_stat(mean_depth)), data = intervals, id = "depth",
#'     stat = "coverage", binwidth = 4) +
#'   coord_collinear()
#' # Fit existing GC windows in source space; estimates are not new data.
#' path <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
#' signal <- read_bedgraph(path)
#' signal$pos <- (signal$start + signal$end) / 2
#' ggchord(data.frame(accver = "pUC19", length = 2686)) +
#'   geom_seq() +
#'   geom_track(aes(x = pos, y = value), data = signal, geom = "line",
#'     stat = "smooth", method = "lm", formula = y ~ x, se = FALSE) +
#'   coord_collinear()
geom_track <- function(mapping = NULL, data = NULL, geom = "point", stat = NULL,
    track = NULL, position = NULL, bounds = c("error", "clip"), baseline = 0,
    show.legend = NA, inherit.aes = TRUE, ..., id = NULL, measure = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geoms <- c("point", "line", "path", "step", "col", "bar", "tile", "rect", "area",
    "ribbon", "segment", "linerange", "errorbar", "text")
  if (!is.character(geom) || length(geom) != 1L || !geom %in% geoms)
    ggchord_stop("geom_track(): unsupported geom; see the documented support matrix")
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
  position <- position %||% if (stat != "coverage" && geom %in% c("bar", "col", "area")) "stack" else "identity"
  estimates <- list(smooth = c("line", "path", "point", "ribbon"),
    density = c("line", "path", "point", "area"), ecdf = c("line", "path", "point", "step"))
  if (ggchord_local_position(position))
    ggchord_stop("geom_track(): position adjusts source data; use feature placement on a feature layer")
  identity_position <- ggchord_is_data_identity(position)
  if (stat %in% names(estimates) && (!geom %in% estimates[[stat]] || !identity_position))
    ggchord_stop("geom_track(): ", stat, " requires position='identity' and geom in ", paste(estimates[[stat]], collapse = ", "))
  dots <- ggchord_colour_dots(list(...), "geom_track()")
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
  if (stat == "coverage") {
    if (!geom %in% c("bar", "col", "point", "line", "step", "tile") || !identity_position)
      ggchord_stop("Coverage requires position='identity' and geom bar, col, point, line, step or tile")
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
  mapping <- ggchord_normalize_mapping(mapping) %||% ggplot2::aes()
  mapping <- ggchord_track_mapping_aliases(mapping)
  ggchord_check_track_styles(names(mapping), fixed, id %||% paste0("geom_track(", geom, ")"))
  if (!is.numeric(baseline) || length(baseline) != 1L || !is.finite(baseline)) ggchord_stop("Track baseline must be finite")
  if (baseline != 0 && !(geom == "col" && stat == "identity" && identity_position))
    ggchord_stop("geom_track(): a nonzero baseline needs col/identity and position='identity'")
  l <- ggchord_track_layer("line", NULL, data, track, "#86AFB9",
    if (geom %in% c("point", "line", "path", "step", "segment", "linerange", "errorbar", "text")) "#387D88" else NA,
    1, .4, 1, bounds, 0, FALSE, show.legend, FALSE,
    c(fill = FALSE, colour = FALSE, alpha = FALSE, linewidth = FALSE, linetype = FALSE), dots[intersect(c("key_glyph", "name"), names(dots))])
  l$layout <- dots[["layout"]]
  dots$layout <- NULL
  l$ggchord_params$id <- id
  l$ggchord_params$native_geom <- geom
  l$ggchord_params$native_stat <- stat
  l$ggchord_params$custom_stat <- custom_stat
  l$ggchord_params$measure <- ggchord_measure_name(measure)
  l$ggchord_params$native_position <- position
  l$ggchord_params$native_mapping <- mapping
  l$ggchord_params$native_params <- dots
  l$ggchord_params$fixed_aesthetics <- fixed
  l$ggchord_params$native_inherit <- inherit.aes
  l$ggchord_params$baseline <- baseline
  l
}

#' @rdname geom_track
#' @export
geom_track_col <- function(mapping = NULL, data = NULL, position = NULL,
    track = NULL, ..., stat = "identity", id = NULL, show.legend = NA, inherit.aes = TRUE) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping, data, geom = "col", stat = stat, track = track, id = id,
    position = position, ..., show.legend = show.legend, inherit.aes = inherit.aes)
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

ggchord_check_track_styles <- function(mapped, fixed, id) {
  both <- intersect(mapped, fixed)
  if (length(both)) ggchord_stop("Layer '", id, "' both maps and fixes ", paste(both, collapse = ", "))
}

ggchord_track_reference_index <- function(input, sequences, caller) {
  ggchord_require_columns(input, "accver", caller)
  ids <- as.character(input$accver)
  index <- match(ids, sequences$accver)
  bad <- is.na(index) | is.na(ids) | !nzchar(ids)
  if (any(bad)) ggchord_stop(caller, ": unknown reference accession in accver ",
    paste(utils::head(unique(ids[bad]), 5L), collapse = ", "), ggchord_input_rows(bad),
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
        rlang::as_label(start) %in% c("pos", "start", "end"))
      value <- as.numeric(as.character(value))
    if (key == "weight") integer_weight <- is.integer(value)
    column <- paste0(".ggchord_start_", key)
    replacement <- NULL
    if (!is.null(value)) {
      data[[column]] <- if (!nrow(data)) utils::head(value, 0L) else value
      replacement <- rlang::sym(column)
    }
    expr <- rlang::get_expr(mapping[[key]])
    if (rlang::is_call(expr, "stage")) {
      expr <- match.call(ggplot2::stage, expr)
      expr["start"] <- list(replacement)
      native_mapping[[key]] <- rlang::quo_set_expr(mapping[[key]], expr)
    } else native_mapping[[key]] <- rlang::new_quosure(replacement, emptyenv())
  }
  list(data = data, mapping = native_mapping, integer_weight = integer_weight)
}

ggchord_track_native_input <- function(layer, plot, sequences) {
  coverage <- identical(layer$ggchord_params$native_stat, "coverage")
  fallback <- ggchord_layer_fallback(layer, plot$ggchord$data)
  raw <- ggchord_layer_source(layer, fallback)
  tryCatch(ggchord_compute_track_input(layer, plot, raw, sequences, coverage),
    error = function(e) rlang::abort(ggchord_layer_input_context(layer), parent = e, call = NULL))
}

# Only an unmodified computed column identifies the units of a height.
ggchord_coverage_measure <- function(height) {
  if (is.null(height)) return("coverage")
  measures <- c("coverage", "mean_depth")
  if (rlang::is_symbol(height)) {
    legacy <- paste0("..", measures, "..")
    key <- rlang::as_string(height)
    if (key %in% legacy) return(measures[match(key, legacy)])
    return(NULL)
  }
  column <- if (rlang::is_call(height, "stage")) match.call(ggplot2::stage, height)$after_stat else
    if (rlang::is_call(height, c("after_stat", "stat"), n = 1L)) height[[2L]]
  while (rlang::is_call(column, "(", n = 1L)) column <- column[[2L]]
  if (rlang::is_call(column, c("$", "[["), n = 2L) && identical(column[[2L]], quote(.data))) {
    if (rlang::is_call(column, "[[") && !rlang::is_string(column[[3L]])) return(NULL)
    key <- column[[3L]]
  } else {
    if (!rlang::is_symbol(column)) return(NULL)
    key <- column
  }
  if (rlang::is_symbol(key)) key <- rlang::as_string(key)
  if (rlang::is_string(key) && key %in% measures) key
}

ggchord_compute_track_input <- function(layer, plot, raw, sequences, coverage) {
  p <- layer$ggchord_params
  facet_keys <- ggchord_facet_vars(plot)
  if(!isTRUE(plot$coordinates$ggchord_coord))ggchord_stop("geom_track(): use a ggchord coordinate; for vertical backbones use coord_collinear(seq_axis='vertical') instead of coord_flip()")
  if (coverage && is.null(raw)) ggchord_stop("Coverage needs interval data or a plot feature table")
  ggchord_require_columns(raw, character(), "geom_track()")
  ggchord_require_columns(raw, facet_keys, "Faceted track observations")
  source_data <- raw
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
      ggchord_check_track_styles(names(mapping), fixed, layer$ggchord_layer_id)
    }
  }
  if (isTRUE(p$native_inherit)) {
    inherited <- ggchord_track_mapping_aliases(plot$mapping)
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
    height <- rlang::get_expr(mapping$y)
    known <- ggchord_coverage_measure(height)
    if (p$native_geom == "tile" && (!is.null(height) || "height" %in% c(names(mapping), names(p$native_params))))
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
    mapping$xmin <- rlang::quo(.coverage_start)
    mapping$xmax <- rlang::quo(.coverage_end)
  }
  for(field in if (!coverage) intersect(c("pos","start","end"),names(raw))) {
    if(is.factor(raw[[field]]) || is.character(raw[[field]])) {
      value <- suppressWarnings(as.numeric(as.character(raw[[field]])))
      if(any(is.na(value) & !is.na(raw[[field]])))ggchord_stop("geom_track(): ",field," must contain numeric coordinate labels")
      raw[[field]]<-value
    }
  }
  interval_x <- FALSE
  if (!coverage && !"x" %in% names(mapping) && p$native_geom != "rect") {
    if ("pos" %in% names(raw)) mapping$x <- rlang::quo(pos) else if (all(c("start", "end") %in% names(raw))) {
      interval_x <- TRUE
      raw$.track_x <- (raw$start + raw$end - 1) / 2
      raw$.track_width <- raw$end - raw$start + 1
      mapping$x <- rlang::quo(.track_x)
    } else ggchord_stop("geom_track(): map x or supply pos/start/end")
  }
  if (!coverage && !"y" %in% names(mapping) && "track_value" %in% names(raw)) mapping$y <- rlang::quo(track_value)
  if (p$native_stat == "count" && "y" %in% names(mapping) &&
      !is.null(ggchord_start_mapping(mapping$y))) ggchord_stop(
    "geom_track_bar(): count does not accept supplied heights; use geom_track_col() or explicit stat='identity'")
  if (p$native_geom == "tile" && !coverage && !"y" %in% names(mapping)) mapping$y <- rlang::quo(0.5)
  late <- names(mapping)[vapply(mapping, ggchord_has_stage, logical(1), stages = "after_scale")]
  late_mapping <- mapping[late]
  for (key in late) {
    expr <- rlang::get_expr(mapping[[key]])
    if (rlang::is_call(expr, "stage")) {
      expr <- match.call(ggplot2::stage, expr)
      expr$after_scale <- NULL
      if (!is.null(expr$after_stat) || !is.null(expr$start)) {
        mapping[[key]] <- rlang::quo_set_expr(mapping[[key]], expr)
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
    mapping$group <- rlang::quo(.track_group)
    params$group <- NULL
  }
  if (!coverage && !is.null(raw[[".track_width"]]) && is.null(params$width) &&
      is.null(mapping$width) && p$native_geom %in% c("bar", "col", "tile")) {
    mapping$width <- rlang::quo(.track_width)
  }
  # Partition source statistics by accession and facet values in ONE native
  # layer. ggplot2 evaluates after_stat() after all panels have been combined.
  domain_keys <- unique(c("accver", facet_keys))
  domain_id <- do.call(interaction, c(unname(lapply(reference_data[domain_keys], function(x)
    match(x, unique(x)))), list(drop = TRUE, lex.order = TRUE)))
  raw$.track_domain <- as.integer(domain_id)
  input_domain <- raw$.track_domain
  raw$.track_input_index <- seq_len(nrow(raw))
  if (length(facet_keys) && !is.null(layer$layout)) {
    ggchord_check_panel_layout(layer)
    panel_map <- ggchord_facet_panels(plot, ggchord_facet_frames(plot))
    rows <- raw[facet_keys]
    attr(rows, "layout") <- layer$layout
    rows$.input <- raw$.track_input_index
    kept <- plot$facet$map_data(rows, panel_map$panels, panel_map$params)$.input
    raw <- raw[kept, , drop = FALSE]
  }
  domains <- reference_data[raw$.track_input_index, domain_keys, drop = FALSE]
  domains$.track_domain <- raw$.track_domain
  domains <- domains[!duplicated(domains$.track_domain), , drop = FALSE]
  stat <- p$custom_stat %||% p$native_stat
  if (!is.null(p$custom_stat) && any(c("source_row", "source_rows") %in% names(mapping)))
    ggchord_stop("Custom track statistics reserve source_row and source_rows for provenance")
  if (coverage) {
    stat <- StatGgchordCoverage
    params$window_grid <- ggchord_coverage_grid(sequences, params$binwidth, params$windows,
      unique(raw$.track_reference))
    params$binwidth <- params$windows <- NULL
    mapping$.ggchord_reference <- rlang::quo(.track_reference)
  }
  source <- if (length(facet_keys) && !is.null(layer$layout))
    source_data[raw$.track_input_index, , drop = FALSE] else source_data
  start_aes <- ggchord_track_start_aesthetics(mapping, native_roles, source, raw, caller)
  if (nrow(raw)) {
    # Identity scales retain source aesthetic values; final role scales train
    # only after projection. No device geometry enters these statistics.
    native <- do.call(getExportedValue("ggplot2", paste0("geom_", p$native_geom)),
      c(list(mapping = start_aes$mapping, data = start_aes$data, stat = stat,
        position = p$native_position, inherit.aes = FALSE), params))
    if (p$native_stat != "identity") {
      # Preserve Stat coordinates before a data Position moves the display
      # centre/edges. Source membership must not follow a dodge displacement.
      source_geom <- native$geom
      native$geom <- ggplot2::ggproto(NULL,source_geom)
      native$geom$default_aes$.ggchord_stat_x <- NA_real_
      native$mapping$.ggchord_stat_x <- ggplot2::aes(value=after_stat(x))$value
    }
    lineage_key <- if (is.null(p$custom_stat)) ".ggchord_source_index" else "source_row"
    native$mapping[[lineage_key]] <- rlang::quo(.track_input_index)
    source_stat <- native$stat
    native$stat <- ggplot2::ggproto(NULL, source_stat,
      required_aes = unique(c(source_stat$required_aes, lineage_key)),
      dropped_aes = unique(c(source_stat$dropped_aes, lineage_key)))
    if (p$native_stat == "bin") native$stat <- ggchord_track_bin_stat(native$stat)
    if (p$native_stat %in% c("smooth", "density", "ecdf"))
      native$stat <- ggchord_track_group_stat(native$stat)
    if (p$native_stat %in% c("count", "summary"))
      native$stat <- ggchord_track_x_stat(native$stat)
    failures <- capture %||% new.env(parent = emptyenv())
    native <- ggchord_track_source_layer(native, failures, implicit_group, interval_x, p$bounds)
    q <- ggplot2::ggplot() + native + ggplot2::scale_fill_identity() +
      ggplot2::scale_colour_identity() + ggplot2::scale_alpha_identity() +
      ggplot2::scale_size_identity() + ggplot2::scale_linewidth_identity() +
      ggplot2::scale_continuous_identity(aesthetics=".ggchord_stat_x")
    if (nrow(domains) > 1L)
      q <- q + ggplot2::facet_wrap(ggplot2::vars(.track_domain), scales = "free")
    built <- ggplot2::ggplot_build(q)
    if (is.null(capture) && length(failures$statistic_errors))
      ggchord_stop("Statistic failed: ", paste(failures$statistic_errors, collapse = "\n"))
    computed <- built$data[[1L]]
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
    if (is.null(p$custom_stat)) members <- lapply(members, as.integer)
    computed$source_rows <- I(members)
    if (!is.null(p$custom_stat)) {
      computed$source_row <- NULL
      computed$source_status <- ifelse(vapply(members, is.null, logical(1)), "unknown", "reported")
    }
    if (coverage && nrow(computed)) {
      computed$start <- computed$window_start - 1
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
    mapping <- late_mapping <- ggplot2::aes()
  }
  facet_fields <- setdiff(facet_keys, "accver")
  conflicts <- intersect(facet_fields, c(names(output), ".track_domain", ".track_input_index"))
  if (length(conflicts)) ggchord_stop("Facet column conflicts with track output: ",
    paste(conflicts, collapse = ", "), ". Copy it to a sample column before faceting.")
  for (key in facet_fields)
    output[[key]] <- if (nrow(output)) domains[[key]][at] else domains[[key]][FALSE]
  visual <- ggplot2::aes()
  for (key in intersect(c("fill", "colour", "alpha", "linewidth", "linetype", "size", "shape", "label", "angle", "hjust", "vjust", "family", "fontface", "lineheight", "stroke"), names(mapping))) {
    column <- paste0(".track_", key)
    output[[column]] <- output[[key]]
    target <- if (key %in% c("fill", "colour")) ggchord_track_aes(key, p$track_id) else key
    visual[[target]] <- rlang::new_quosure(rlang::sym(column), emptyenv())
  }
  for (key in names(late_mapping)) {
    target <- if (key %in% c("fill", "colour")) ggchord_track_aes(key, p$track_id) else key
    q <- ggchord_track_stage_mapping(late_mapping[key], p$track_id)[[key]]
    expr <- rlang::get_expr(q)
    if (rlang::is_call(expr, "stage")) {
      expr$after_stat <- NULL
      expr$start <- if (key %in% names(mapping)) rlang::sym(paste0(".track_", key)) else NULL
      q <- rlang::quo_set_expr(q, expr)
    }
    visual[[target]] <- q
  }
  # Fixed styles have already been interpreted by the native geom.
  for (key in intersect(c("fill", "colour", "alpha", "linewidth", "linetype", "size", "shape", "angle", "hjust", "vjust", "family", "fontface", "lineheight", "stroke"), names(params))) {
    target <- if (key %in% c("fill", "colour")) ggchord_track_aes(key, p$track_id) else key
    layer$aes_params[[target]] <- params[[key]]
  }
  layer$ggchord_input_mapping <- visual
  layer$ggchord_params$track_mapping_labels <- stats::setNames(
    vapply(mapping, rlang::as_label, character(1)),
    vapply(names(mapping), function(key) if (key %in% c("fill", "colour"))
      ggchord_track_aes(key, p$track_id) else key, character(1)))
  layer$ggchord_input_data <- output
  layer$ggchord_params$track_kind <- if (coverage && p$native_geom == "tile") "tile" else
    if (p$native_geom %in% c("bar", "col")) "bar" else "line"
  layer$ggchord_params$source_edges <- TRUE
  layer$ggchord_params$track_sources <- source_data
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
  layer$stat <- ggplot2::ggproto(NULL, stat,
    parameters = function(extra = FALSE) stat$parameters(extra),
    compute_panel = function(self, data, scales, ...) {
      panel_error <<- NULL
      withCallingHandlers(ggplot2::ggproto_parent(stat, self)$compute_panel(data, scales, ...),
        error = function(e) panel_error <<- e)
    })
  # ggplot2 attaches the escaping panel error as its warning's parent. An
  # extension may instead recover in compute_layer() and return valid data.
  ggplot2::ggproto(NULL, layer,
    compute_aesthetics = function(self, data, plot) {
      mapped <- ggplot2::ggproto_parent(layer, self)$compute_aesthetics(data, plot)
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
        if (any(vapply(discrete, anyNA, logical(1))))
          ggchord_stop("geom_track(): group values must be non-missing")
        # Factor codes preserve level order without joining user labels with a separator.
        codes <- lapply(discrete, function(value) as.integer(as.factor(value)))
        mapped$group <- if (length(codes)) as.integer(do.call(interaction,
          c(unname(codes), list(drop = TRUE, lex.order = TRUE)))) else rep(1L, nrow(mapped))
      }
      mapped
    },
    compute_statistic = function(self, data, layout) {
      tryCatch(withCallingHandlers(
        ggplot2::ggproto_parent(layer, self)$compute_statistic(data, layout),
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
  ggplot2::ggproto(NULL, stat,
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
  ggplot2::ggproto(NULL, stat,
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
  ggplot2::ggproto(NULL, stat,
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
      centers <- (utils::head(resolved$breaks, -1L) + utils::tail(resolved$breaks, -1L)) / 2
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
          add(c(top$x,bottom$x), c(top$y,bottom$y), c(top$rows,bottom$rows), "track_poly")
        } else add(top$x, top$y, top$rows, "track_line")
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
