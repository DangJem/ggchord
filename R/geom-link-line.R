ggchord_link_route <- function(route, dots, caller) {
  route <- ggchord_match_arg(route, c("curve", "straight", "shared"))
  curvature <- dots$curvature %||% 0.5
  bulge <- dots$bulge %||% 0.2
  corridor <- ggchord_physical_mm(dots$corridor_width %||% 2, "corridor_width", TRUE)
  if (!is.numeric(curvature) || length(curvature) != 1L || !is.finite(curvature) ||
      curvature < 0 || curvature > 1 || !is.numeric(bulge) || length(bulge) != 1L ||
      !is.finite(bulge) || bulge < 0 || length(corridor) != 1L || corridor <= 0)
    ggchord_stop(caller, ": invalid route curvature, bulge, or corridor width")
  if (route == "straight" && any(c("curvature", "bulge") %in% names(dots)))
    ggchord_stop(caller, ": straight links do not use curvature or bulge")
  if (route != "shared" && "corridor_width" %in% names(dots))
    ggchord_stop(caller, ": corridor_width requires route = 'shared'")
  args <- dots[setdiff(names(dots), c("curvature", "bulge", "corridor_width"))]
  structure(list(route = route, curvature = curvature, bulge = bulge,
    explicit_bulge = "bulge" %in% names(dots), corridor_width = corridor, dots = args),
    class = "ggchord_link_route")
}

#' Draw point-to-point sequence links
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param route Choose curved or straight routes, or a shared corridor where supported.
#'   On a curved collinear route, an explicit positive `bulge` clears the high
#'   coordinate end of intervening rows by that fraction of the longest visible
#'   backbone. Radial curves with positive `bulge` clear intermediate spoke tips.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Source references are `qaccver` and `saccver`. Point links consume `qpos` and `spos`. Interval links and alignment segments consume closed interval endpoints on both references. Supply reliable `qstrand` and `sstrand` or explicit relative direction when the geometry requires it.
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' points <- data.frame(qaccver = "A", qpos = 200, saccver = "B", spos = 220)
#'
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_link_line(data = points, route = "curve", colour = "#187B77") +
#'   coord_chord()
geom_link_line <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", route = "curve", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  caller <- "geom_link_line()"
  controls <- ggchord_link_route(route, ggchord_colour_dots(list(...), caller), caller)
  dots <- controls$dots
  controls$dots <- NULL
  position <- ggchord_as_feature_position(position, caller)
  if (ggchord_local_position(position) && !identical(position$placement, "offset"))
    ggchord_stop(caller, ": use identity or position_offset()")
  ggchord_check_geom_params(dots, GeomChordLinkLine, caller, stat = stat)
  parameters <- ggchord_split_source_stat(stat, dots, GeomChordLinkLine)
  layer <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), group = integer()),
    mapping = aes(x = x, y = y, group = group), geom = GeomChordLinkLine,
    stat = "identity", position = "identity", inherit.aes = inherit.aes,
    show.legend = show.legend, check.aes = FALSE,
    params = c(list(na.rm = na.rm), parameters$geom)
  )
  layer$ggchord_type <- "link"
  layer$ggchord_params <- list(id = id, type = "link", route = controls,
    link_type = if (controls$route == "straight") "straight" else "curve",
    link_gap = 0, link_avoid = "none", feature_position = position,
    physical_layout = ggchord_local_position(position) || controls$route == "shared",
    link_branch = if (controls$route == "shared") "query" else "none",
    link_branch_fraction = 0.15, source_stat = parameters$stat)
  ggchord_capture_layer_input(layer, data, mapping, c("qaccver", "saccver", "qpos", "spos"))
}

GeomChordLinkLine <- ggproto("GeomChordLinkLine", ggplot2::GeomPath,
  rename_size = FALSE,
  default_aes = aes(colour = "grey40", alpha = 1, linewidth = .5, linetype = 1),
  draw_panel = function(data, panel_params, coord, arrow = NULL, lineend = "round",
                        linejoin = "round", linemitre = 10, na.rm = FALSE) {
    for (nm in c("colour", "alpha", "linewidth", "linetype")) data[[nm]] <- data[[nm]] %||% data[[paste0("link_", nm)]]
    if (".arrow_first" %in% names(data)) {
      pieces <- lapply(split(data, data$group), function(d) {
        ar <- arrow
        if (!is.null(ar)) {
          first <- d$.arrow_first[1] && ar$ends %in% c(1,3)
          last <- d$.arrow_last[1] && ar$ends %in% c(2,3)
          if (!first && !last) ar <- NULL else ar$ends <- if (first && last) 3L else if (first) 1L else 2L
        }
        ggplot2::GeomPath$draw_panel(d, panel_params, coord, arrow=ar,
          lineend=lineend,linejoin=linejoin,linemitre=linemitre,na.rm=na.rm)
      })
      return(do.call(grobTree,pieces))
    }
    ggplot2::GeomPath$draw_panel(data, panel_params, coord, arrow = arrow,
      lineend = lineend, linejoin = linejoin, linemitre = linemitre, na.rm = na.rm)
  },
  draw_key = function(data, params, size) {
    for (nm in c("colour", "alpha", "linewidth", "linetype")) data[[nm]] <- data[[nm]] %||% data[[paste0("link_", nm)]]
    ggplot2::draw_key_path(data, params, size)
  })

ggchord_check_branch_fraction <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x > 0.5)
    ggchord_stop("link_branch_fraction must be in (0, 0.5]")
}

ggchord_rotate_points <- function(x, degrees) {
  a <- degrees * pi / 180
  cbind(x[, 1] * cos(a) - x[, 2] * sin(a), x[, 1] * sin(a) + x[, 2] * cos(a))
}

ggchord_link_controls <- function(control, row, fallback) {
  cp <- if (is.list(control)) control[[min(row, length(control))]] else control
  if (is.list(cp)) cp <- unlist(cp, use.names = FALSE)
  if (is.null(cp)) return(list(fallback, fallback))
  if (!is.numeric(cp) || !length(cp) %in% c(2L, 4L) || any(!is.finite(cp)))
    ggchord_stop("link_ctrl_point must contain two or four finite numbers per link")
  list(cp[1:2], if (length(cp) == 2L) cp else cp[3:4])
}

# A vector or a one-element list applies to every row. A shorter list extends
# its last entry, as it does when rendering the unsliced layer.
ggchord_link_control_rows <- function(control, rows) {
  if (is.list(control) && !inherits(control, "ggchord_link_route") && length(control) > 1L) control[pmin(rows, length(control))]
  else control
}

ggchord_connection_offsets <- function(params, sequences, axis = NULL) {
  offset <- ggchord_resolve_position_offset(params$feature_position$offset %||% 0,
    sequences, "offset", "connection position")
  units <- ggchord_units_per_mm(axis)
  vapply(offset, function(value) unname(value[1L]) * units, numeric(1))
}

ggchord_interval_front <- function(start, end) {
  if (start <= end) c(start - 1, end) else c(start, end - 1)
}

ggchord_link_geometry <- function(data, params, layout) {
  if (is.null(data)) ggchord_stop("geom_link_line(): supply data with qaccver/saccver/qpos/spos")
  ggchord_require_columns(data, c("qaccver", "saccver", "qpos", "spos"), "geom_link_line()")
  ref <- layout$sequence_reference
  rows <- which(data$qaccver %in% names(ref$lens) & data$saccver %in% names(ref$lens))
  for (pair in list(c("qaccver", "qpos"), c("saccver", "spos"))) {
    ids <- as.character(data[[pair[1]]][rows]); pos <- data[[pair[2]]][rows]
    if (!is.numeric(pos) || any(!is.finite(pos)) || any(pos < 1 | pos > ref$lens[ids]))
      ggchord_stop("geom_link_line(): positions must be finite and within sequence length")
  }
  if (isTRUE(layout$collinear)) return(ggchord_collinear_link_geometry(data, params, layout, rows))
  gaps <- ggchord_connection_offsets(params, layout$seqs)
  if (any(!is.finite(unlist(gaps)))) ggchord_stop("link_gap must contain finite numbers")
  endpoints <- function(ids, positions) {
    xy <- matrix(NA_real_, nrow = length(ids), ncol = 2L)
    for (id in unique(ids)) {
      at <- which(ids == id)
      gap <- unname(gaps[[id]])
      if (is.null(params$link_gap) && params$link_avoid != "none") {
        gap <- vapply(positions[at], ggchord_gap_profile, numeric(1),
          id = id, obstacles = layout$obstacles, avoid = "uniform", close_gap = gap)
      }
      xy[at, ] <- ggchord_track_point(positions[at], id, gap, 0, layout)
    }
    xy
  }
  query <- endpoints(as.character(data$qaccver[rows]), data$qpos[rows])
  subject <- endpoints(as.character(data$saccver[rows]), data$spos[rows])
  out <- lapply(seq_along(rows), function(j) {
    i <- rows[j]
    q <- query[j, ]; s <- subject[j, ]
    curvature <- params$route$curvature %||% .5
    cp <- (q + s) * (1 - curvature) / 2
    xy <- if (params$link_type == "straight") rbind(q, s) else bezier_pts(q, s, cp, cp, n = 60)
    data.frame(x = xy[,1], y = xy[,2], group = i, source_row = i, .component = "branch")
  })
  if (!length(out)) return(data.frame(x = numeric(), y = numeric(), group = integer(), source_row = integer()))
  do.call(rbind, out)
}
