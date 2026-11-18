# Native statistical grammar runs on source coordinates, before any projection.
#' Draw a sequence-local data layer with ggplot2 grammar
#'
#' Statistics and data positions are calculated independently within each
#' accession/window before projection into the physical track band. Supported
#' geoms are point, line, path, step, col, bar, tile, rect, area, ribbon,
#' segment, linerange, errorbar and text. Supported statistics are identity,
#' count, bin, summary, smooth, density and ecdf; unsupported combinations fail
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
#' An explicit group mapping takes priority over discrete style mappings.
#' Bin source membership follows the native statistic, including outer edges
#' and its numerical boundary tolerance, before data positions move the bars.
#' Mapping x/y refers to source position/local numerical value, not screen x/y.
#' Existing track_value/track_fill/track_colour roles remain accepted.
#' @param mapping,data Standard layer input. Supply accver and source positions
#'   through x/pos, or interval start/end. Map y/track_value for identity values.
#' @param geom,stat Native geometry and statistic names. The seven supported
#'   statistics also accept their built-in ggplot2 Stat objects, for example
#'   `ggplot2::StatBin`. Bar defaults to count, other geoms to identity.
#'   Custom Stat subclasses require a source-membership rule and are unsupported.
#' @param track A physical [position_track()] declaration.
#' @param position Native data adjustment such as `"identity"`, `"stack"` or
#'   `ggplot2::position_dodge(width = 10)`. Source x uses sequence units.
#'   Feature placement belongs on feature layers. The legacy
#'   `position = position_track(...)` spelling still declares the physical band.
#' @param bounds Coordinate error/clip policy.
#' @param baseline Numeric identity-column baseline. Nonzero baselines require
#'   position='identity' so stacking does not change their meaning.
#' @param fill,colour,alpha,linewidth,linetype Fixed styles for the count shortcut.
#' @param show.legend,inherit.aes Standard layer controls.
#' @param ... Native geom/stat parameters and fixed styles.
#' @return A ggplot2 layer using a shared sequence-local band.
#' @export
#' @examples
#' # Synthetic observations illustrate counting, not biological measurements.
#' events <- data.frame(accver = "A", pos = c(100, 100, 400))
#' ggchord(data.frame(accver = "A", length = 1000)) +
#'   geom_seq() +
#'   geom_track_bar(data = events, track = position_track("count"), stat = "count") +
#'   coord_collinear()
#' # Fit existing GC windows in source space; estimates are not new data.
#' path <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
#' signal <- read_bedgraph(path)
#' signal$pos <- (signal$start + signal$end) / 2
#' ggchord(data.frame(accver = "pUC19", length = 2686)) +
#'   geom_seq() +
#'   geom_track(aes(x = pos, y = value), data = signal, geom = "line",
#'     stat = "smooth", method = "lm", formula = y ~ x, se = FALSE,
#'     track = position_track("GC")) +
#'   coord_collinear()
geom_track <- function(mapping = NULL, data = NULL, geom = "point", stat = NULL,
    track = NULL, position = NULL, bounds = c("error", "clip"), baseline = 0,
    show.legend = NA, inherit.aes = FALSE, ...) {
  geoms <- c("point", "line", "path", "step", "col", "bar", "tile", "rect", "area",
    "ribbon", "segment", "linerange", "errorbar", "text")
  if (!is.character(geom) || length(geom) != 1L || !geom %in% geoms)
    ggchord_stop("geom_track(): unsupported geom; see the documented support matrix")
  stat <- stat %||% if (geom == "bar") "count" else "identity"
  stats <- c(identity = "StatIdentity", count = "StatCount", bin = "StatBin",
    summary = "StatSummary", smooth = "StatSmooth", density = "StatDensity", ecdf = "StatEcdf")
  if (inherits(stat, "Stat")) {
    stat <- names(stats)[vapply(stats, function(name)
      identical(stat, getExportedValue("ggplot2", name)), logical(1))]
  }
  if (!is.character(stat) || length(stat) != 1L || !stat %in% names(stats))
    ggchord_stop("geom_track(): supported stats are identity, count, bin, summary, smooth, density, ecdf")
  if (inherits(position, "Position") && isTRUE(position$ggchord_track_position)) {
    if (!is.null(track)) ggchord_stop("geom_track(): declare the physical track only once")
    track <- position; position <- NULL
  }
  if (is.null(track)) ggchord_stop("geom_track(): supply track = position_track(...)")
  position <- position %||% if (geom %in% c("bar", "col", "area")) "stack" else "identity"
  estimates <- list(smooth = c("line", "path", "point", "ribbon"),
    density = c("line", "path", "point", "area"), ecdf = c("line", "path", "point", "step"))
  if (ggchord_local_position(position))
    ggchord_stop("geom_track(): position adjusts source data; use feature placement on a feature layer")
  identity_position <- ggchord_is_data_identity(position)
  if (stat %in% names(estimates) && (!geom %in% estimates[[stat]] || !identity_position))
    ggchord_stop("geom_track(): ", stat, " requires position='identity' and geom in ", paste(estimates[[stat]], collapse = ", "))
  dots <- ggchord_colour_dots(list(...), "geom_track()")
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
  if (!is.numeric(baseline) || length(baseline) != 1L || !is.finite(baseline)) ggchord_stop("Track baseline must be finite")
  if (baseline != 0 && !(geom == "col" && stat == "identity" && identity_position))
    ggchord_stop("geom_track(): a nonzero baseline needs col/identity and position='identity'")
  l <- ggchord_track_layer("line", NULL, data, track, "#86AFB9",
    if (geom %in% c("point", "line", "path", "step", "segment", "linerange", "errorbar", "text")) "#387D88" else NA,
    1, .4, 1, bounds, 0, FALSE, show.legend, FALSE,
    c(fill = FALSE, colour = FALSE, alpha = FALSE, linewidth = FALSE, linetype = FALSE), list())
  l$ggchord_params$native_geom <- geom
  l$ggchord_params$native_stat <- stat
  l$ggchord_params$native_position <- position
  l$ggchord_params$native_mapping <- mapping
  l$ggchord_params$native_params <- dots
  l$ggchord_params$native_inherit <- inherit.aes
  l$ggchord_params$baseline <- baseline
  l
}

#' @rdname geom_track
#' @export
geom_track_col <- function(mapping = NULL, data = NULL, position = NULL,
    track = NULL, ..., show.legend = NA, inherit.aes = FALSE) {
  geom_track(mapping, data, geom = "col", stat = "identity", track = track,
    position = position, ..., show.legend = show.legend, inherit.aes = inherit.aes)
}

ggchord_track_native_input <- function(layer, plot) {
  p <- layer$ggchord_params
  if(!inherits(plot$facet,"FacetNull"))ggchord_stop("geom_track(): native tracks currently require one ggplot panel; use separate ggchord plots for faceted statistics")
  if(!isTRUE(plot$coordinates$ggchord_coord))ggchord_stop("geom_track(): use a ggchord coordinate; for vertical backbones use coord_collinear(seq_axis='vertical') instead of coord_flip()")
  raw <- layer$ggchord_input_data %||% plot$ggchord$data$seq_data
  if (is.function(raw)) raw <- raw(plot$ggchord$data$seq_data)
  ggchord_require_columns(raw, character(), "geom_track()")
  capture <- p$stage_capture
  if (!is.null(capture) && "source" %in% capture$requested && is.null(capture$source))
    capture$source <- raw
  mapping <- p$native_mapping
  if (isTRUE(p$native_inherit)) for (key in setdiff(names(plot$mapping), names(mapping)))
    mapping[[key]] <- plot$mapping[[key]]
  for (key in intersect(c("accver", "group"), names(mapping))) {
    if (ggchord_has_stage(mapping[[key]]))
      ggchord_stop("geom_track(): ", key, " must be available before the source statistic")
  }
  if ("accver" %in% names(mapping)) {
    ids <- rlang::eval_tidy(mapping$accver, raw)
    if (!length(ids) %in% c(1L, nrow(raw)))
      ggchord_stop("geom_track(): mapped accver must be scalar or supply one accession per input row")
    raw$accver <- rep_len(ids, nrow(raw))
    mapping$accver <- NULL
  }
  ggchord_require_columns(raw, "accver", "geom_track()")
  aliases <- c(pos = "x", track_value = "y", track_fill = "fill", track_colour = "colour")
  for (key in names(aliases)) if (key %in% names(mapping)) {
    target <- aliases[[key]]
    if (target %in% names(mapping)) ggchord_stop("geom_track(): use only one of ", key, " and ", target)
    mapping[[target]] <- mapping[[key]]; mapping[[key]] <- NULL
  }
  for(field in intersect(c("pos","start","end"),names(raw))) {
    if(is.factor(raw[[field]]) || is.character(raw[[field]])) {
      value <- suppressWarnings(as.numeric(as.character(raw[[field]])))
      if(any(is.na(value) & !is.na(raw[[field]])))ggchord_stop("geom_track(): ",field," must contain numeric coordinate labels")
      raw[[field]]<-value
    }
  }
  interval_x <- FALSE
  if (!"x" %in% names(mapping) && p$native_geom != "rect") {
    if ("pos" %in% names(raw)) mapping$x <- rlang::quo(pos) else if (all(c("start", "end") %in% names(raw))) {
      interval_x <- TRUE
      raw$.track_x <- (raw$start + raw$end - 1) / 2
      raw$.track_width <- raw$end - raw$start + 1
      mapping$x <- rlang::quo(.track_x)
    } else ggchord_stop("geom_track(): map x or supply pos/start/end")
  }
  if("x" %in% names(mapping)) {
    source_x <- tryCatch(rlang::eval_tidy(mapping$x,raw),error=function(e)NULL)
    if(!is.null(source_x)) {
      if(!is.numeric(source_x))ggchord_stop("geom_track(): x must contain numeric source positions; convert factor labels explicitly")
      if(any(is.infinite(source_x)))ggchord_stop("geom_track(): source positions must be finite")
      if(!interval_x && identical(p$bounds,"error") && any(source_x<1,na.rm=TRUE))
        ggchord_stop("geom_track(): source positions are one-based; x/pos must be at least 1")
    }
  }
  if (!"y" %in% names(mapping) && "track_value" %in% names(raw)) mapping$y <- rlang::quo(track_value)
  if (p$native_stat == "count" && "y" %in% names(mapping) &&
      !is.null(ggchord_start_mapping(mapping$y))) ggchord_stop(
    "geom_track_bar(): count does not accept supplied heights; use geom_track_col() or explicit stat='identity'")
  if (p$native_geom == "tile" && !"y" %in% names(mapping)) mapping$y <- rlang::quo(0.5)
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
  source_group <- if ("group" %in% names(mapping)) rlang::eval_tidy(mapping$group, raw) else rep(1L, nrow(raw))
  if (!length(source_group) %in% c(1L, nrow(raw)))
    ggchord_stop("geom_track(): group must be scalar or supply one value per input row")
  group_fields <- list(rep_len(source_group, nrow(raw)))
  for (key in if (!"group" %in% names(mapping))
      intersect(c("fill", "colour", "linetype", "shape"), names(mapping)) else character()) {
    start <- ggchord_start_mapping(mapping[[key]])
    if (is.null(start)) next
    value <- rlang::eval_tidy(start, raw)
    if (!is.numeric(value)) group_fields[[length(group_fields) + 1L]] <- rep_len(value, nrow(raw))
  }
  raw$.track_group <- as.integer(do.call(interaction, c(group_fields, list(drop = TRUE, lex.order = TRUE))))
  if (anyNA(raw$.track_group)) ggchord_stop("geom_track(): group values must be non-missing")
  mapping$group <- rlang::quo(.track_group)
  params <- p$native_params
  if (!is.null(raw[[".track_width"]]) && is.null(params$width) && p$native_geom %in% c("bar", "col", "tile")) {
    mapping$width <- rlang::quo(.track_width)
  }
  parts <- list()
  for (seq in unique(as.character(raw$accver))) {
    rows <- which(as.character(raw$accver) == seq)
    d <- raw[rows, , drop = FALSE]
    # Identity scales retain source aesthetic values; final role scales train
    # only after projection. No device geometry enters these statistics.
    native <- do.call(getExportedValue("ggplot2", paste0("geom_", p$native_geom)),
      c(list(mapping = mapping, data = d, stat = p$native_stat,
        position = p$native_position, inherit.aes = FALSE), params))
    if (p$native_stat != "identity") {
      # Preserve Stat coordinates before a data Position moves the display
      # centre/edges. Source membership must not follow a dodge displacement.
      source_geom <- native$geom
      native$geom <- ggplot2::ggproto(NULL,source_geom)
      native$geom$default_aes$.ggchord_stat_x <- NA_real_
      native$mapping$.ggchord_stat_x <- ggplot2::aes(value=after_stat(x))$value
    }
    native$data$.track_input_index <- seq_len(nrow(d))
    native$mapping$.ggchord_source_index <- rlang::quo(.track_input_index)
    source_stat <- native$stat
    required <- unique(c(source_stat$required_aes, ".ggchord_source_index"))
    native$stat <- ggplot2::ggproto(NULL, source_stat, required_aes = required)
    if (p$native_stat == "bin") native$stat <- ggchord_track_bin_stat(native$stat)
    if (p$native_stat %in% c("smooth", "density", "ecdf"))
      native$stat <- ggchord_track_group_stat(native$stat)
    if (p$native_stat %in% c("count", "summary"))
      native$stat <- ggchord_track_x_stat(native$stat)
    q <- ggplot2::ggplot() + native + ggplot2::scale_fill_identity() +
      ggplot2::scale_colour_identity() + ggplot2::scale_alpha_identity() +
      ggplot2::scale_size_identity() + ggplot2::scale_linewidth_identity() +
      ggplot2::scale_continuous_identity(aesthetics=".ggchord_stat_x")
    computed <- ggplot2::ggplot_build(q)$data[[1L]]
    if (!nrow(computed)) next
    computed$accver <- seq
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
    # Stat adapters return local input indices before any data-position move.
    members <- switch(p$native_stat,
      identity = as.list(computed$.ggchord_source_index),
      bin = computed$.ggchord_bin_rows,
      computed$.ggchord_group_rows)
    computed$source_rows <- I(lapply(members, function(local) rows[local]))
    if (p$native_stat %in% c("count", "summary")) computed$.ggchord_source_index <- NULL
    computed$.ggchord_bin_rows <- NULL
    computed$.ggchord_group_rows <- NULL
    parts[[length(parts) + 1L]] <- computed
  }
  output <- ggchord_rbind_fill(parts)
  if (is.null(output)) output <- data.frame(accver = character(), pos = numeric(), start = numeric(),
    end = numeric(), track_value = numeric(), track_ymin = numeric(), track_ymax = numeric(), group = integer())
  visual <- ggplot2::aes()
  for (key in intersect(c("fill", "colour", "alpha", "linewidth", "linetype", "size", "shape", "label", "angle", "hjust", "vjust", "family", "fontface", "lineheight", "stroke"), names(mapping))) {
    column <- paste0(".track_", key)
    output[[column]] <- output[[key]]
    target <- if (key %in% c("fill", "colour")) ggchord_track_aes(key, p$track_position$track) else key
    visual[[target]] <- rlang::new_quosure(rlang::sym(column), emptyenv())
  }
  for (key in names(late_mapping)) {
    target <- if (key %in% c("fill", "colour")) ggchord_track_aes(key, p$track_position$track) else key
    q <- ggchord_track_stage_mapping(late_mapping[key], p$track_position$track)[[key]]
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
    target <- if (key %in% c("fill", "colour")) ggchord_track_aes(key, p$track_position$track) else key
    layer$aes_params[[target]] <- params[[key]]
  }
  layer$ggchord_input_mapping <- visual
  layer$ggchord_params$track_mapping_labels <- stats::setNames(
    vapply(mapping, rlang::as_label, character(1)),
    vapply(names(mapping), function(key) if (key %in% c("fill", "colour"))
      ggchord_track_aes(key, p$track_position$track) else key, character(1)))
  layer$ggchord_input_data <- output
  layer$ggchord_params$track_kind <- if (p$native_geom %in% c("bar", "col")) "bar" else "line"
  layer$ggchord_params$source_edges <- TRUE
  layer$ggchord_params$track_sources <- raw
  layer$ggchord_params$native_prepared <- TRUE
  layer$geom_params$native_draw_params <- params[intersect(names(params),c("arrow","lineend","linejoin","linemitre","parse","check_overlap","rule"))]
  layer$geom_params$track_kind <- p$native_geom
  layer
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

ggchord_track_native_geometry <- function(d, p, layout, sc, at_value) {
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
      track_id = p$track_position$track,
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
  } else for (i in seq_len(nrow(d))) {
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
