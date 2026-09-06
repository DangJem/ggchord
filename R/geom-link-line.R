#' Draw point-to-point sequence links
#'
#' Requires qaccver, saccver, qpos and spos. Positions use the same genomic
#' coordinates as sequence tracks. Curves run from query to subject; an arrow
#' supplied by [arrow()] therefore points to the subject by default.
#' @param mapping,data Standard ggplot2 mapping and layer data, including data functions.
#' @param link_type Either "curve" (default) or "straight".
#' @param link_gap Explicit endpoint spacing. With the default
#'   \code{link_avoid = "none"}, \code{NULL} uses a fixed spacing of 0.035.
#'   Set \code{link_avoid} to \code{"smooth"} or \code{"uniform"} to opt in
#'   to entity-aware endpoint spacing.
#' @param link_ctrl_point Bezier control points: two numbers, four numbers, or a list per input row.
#' @param link_avoid One of "none" (default), "smooth", or "uniform".
#'   Explicit gap takes precedence.
#' @param link_branch Share equal query or subject endpoints: "none", "query", "subject".
#' @param link_branch_fraction Shared length relative to the shortest member, in (0, 0.5].
#' @param arrow A grid arrow specification, or NULL.
#' @param colour,alpha,linewidth,linetype Fixed styles; mappings use link_* aesthetics.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional path parameters.
#' @return A ggplot2 layer.
#' @export
#' @examples
#' s <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
#' d <- data.frame(qaccver = "A", saccver = "B", qpos = 200, spos = 600)
#' ggchord(s) + geom_seq() + geom_link_line(data = d,
#'   arrow = arrow(length = unit(2, "mm")))
geom_link_line <- function(mapping = NULL, data = NULL,
    link_type = c("curve", "straight"), link_gap = NULL, link_ctrl_point = NULL,
    link_avoid = c("none", "smooth", "uniform"),
    link_branch = c("none", "query", "subject"), link_branch_fraction = 0.15,
    arrow = NULL, colour = "#59636D", alpha = 0.6, linewidth = 0.4,
    linetype = 1, position = "identity", show.legend = NA, inherit.aes = FALSE, ...) {
  link_type <- match.arg(link_type)
  link_avoid <- match.arg(link_avoid)
  link_branch <- match.arg(link_branch)
  ggchord_check_branch_fraction(link_branch_fraction)
  if (link_type == "straight" && (!is.null(link_ctrl_point) || link_branch != "none"))
    ggchord_stop("Straight links cannot use link_ctrl_point or link_branch; use link_type='curve'")
  if (!is.null(arrow) && !inherits(arrow, "arrow")) ggchord_stop("arrow must be a grid::arrow() object")
  mapping <- ggchord_normalize_mapping(mapping)
  alias <- ggchord_colour_alias(colour, list(...), "geom_link_line()", sys.call())
  styles <- list(link_colour = alias$colour, link_alpha = alpha,
                 link_linewidth = linewidth, link_linetype = linetype)
  styles[names(styles) %in% names(mapping)] <- NULL
  lyr <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), group = integer()),
    mapping = ggplot2::aes(x = x, y = y, group = group), geom = GeomChordLinkLine,
    stat = "identity", position = position, inherit.aes = inherit.aes,
    show.legend = show.legend, check.aes = FALSE,
    params = c(styles, list(arrow = arrow), alias$dots))
  lyr$ggchord_type <- "link"
  lyr$ggchord_params <- list(type = "link", link_type = link_type, link_gap = link_gap,
    link_ctrl_point = link_ctrl_point, link_avoid = link_avoid,
    link_branch = link_branch, link_branch_fraction = link_branch_fraction)
  ggchord_capture_layer_input(lyr, data, mapping, c("qaccver", "saccver", "qpos", "spos"))
}

GeomChordLinkLine <- ggplot2::ggproto("GeomChordLinkLine", ggplot2::GeomPath,
  rename_size = FALSE,
  default_aes = ggplot2::aes(link_colour = "#59636D", link_alpha = 0.6,
    link_linewidth = 0.4, link_linetype = 1),
  draw_panel = function(data, panel_params, coord, arrow = NULL, lineend = "round",
                        linejoin = "round", linemitre = 10, na.rm = FALSE) {
    for (nm in c("colour", "alpha", "linewidth", "linetype")) data[[nm]] <- data[[paste0("link_", nm)]]
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
      return(do.call(grid::grobTree,pieces))
    }
    ggplot2::GeomPath$draw_panel(data, panel_params, coord, arrow = arrow,
      lineend = lineend, linejoin = linejoin, linemitre = linemitre, na.rm = na.rm)
  },
  draw_key = function(data, params, size) {
    for (nm in c("colour", "alpha", "linewidth", "linetype")) data[[nm]] <- data[[paste0("link_", nm)]]
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

ggchord_link_geometry <- function(data, params, layout) {
  if (is.null(data)) ggchord_stop("geom_link_line(): supply data with qaccver/saccver/qpos/spos")
  ggchord_require_columns(data, c("qaccver", "saccver", "qpos", "spos"), "geom_link_line()")
  ref <- layout$sequence_reference
  for (pair in list(c("qaccver", "qpos"), c("saccver", "spos"))) {
    ids <- as.character(data[[pair[1]]]); pos <- data[[pair[2]]]
    if (anyNA(ids) || any(!ids %in% names(ref$lens))) ggchord_stop("geom_link_line(): unknown sequence accession")
    if (!is.numeric(pos) || any(!is.finite(pos)) || any(pos < 1 | pos > ref$lens[ids]))
      ggchord_stop("geom_link_line(): positions must be finite and within sequence length")
  }
  gaps <- process_sequence_param(params$link_gap %||% 0.035, layout$seqs, "link_gap", 0.035)
  if (any(!is.finite(unlist(gaps)))) ggchord_stop("link_gap must contain finite numbers")
  endpoint <- function(id, pos) {
    f <- (pos - 1) / ref$lens[id]
    if (ref$orientation[id] != 1) f <- 1 - f
    a <- ref$starts[id] + f * (ref$ends[id] - ref$starts[id])
    gap <- unname(gaps[[id]])
    if (is.null(params$link_gap) && params$link_avoid != "none") {
      gap <- ggchord_gap_profile(pos, id, layout$obstacles, "uniform", close_gap = gap)
    }
    map_to_curve(a, ref$radius[id] + gap, ref$refs[[id]])
  }
  out <- lapply(seq_len(nrow(data)), function(i) {
    q <- endpoint(as.character(data$qaccver[i]), data$qpos[i])
    s <- endpoint(as.character(data$saccver[i]), data$spos[i])
    cp <- ggchord_link_controls(params$link_ctrl_point, i, (q + s) / 4)
    xy <- if (params$link_type == "straight") rbind(q, s) else bezier_pts(q, s, cp[[1]], cp[[2]], n = 60)
    xy <- ggchord_rotate_points(xy, layout$rotation)
    data.frame(x = xy[,1], y = xy[,2], group = i, source_row = i, .component = "branch")
  })
  if (!length(out)) return(data.frame(x = numeric(), y = numeric(), group = integer(), source_row = integer()))
  do.call(rbind, out)
}
