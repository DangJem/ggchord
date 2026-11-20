# coord-chord.R - chord diagram coordinate system

#' Chord-diagram coordinate ggproto
#' @noRd
CoordGgchord <- ggplot2::ggproto(
  "CoordGgchord", ggplot2::CoordCartesian,
  ggchord_coord = TRUE,
  render_axis_h = function(panel_params, theme) list(top = grid::nullGrob(), bottom = grid::nullGrob()),
  render_axis_v = function(panel_params, theme) list(left = grid::nullGrob(), right = grid::nullGrob())
)

CoordChord <- ggplot2::ggproto("CoordChord", CoordGgchord)

new_ggchord_coord <- function(proto, rotation, ratio, xlim, ylim, expand,
                              clip, fit, user_xlim = xlim,
                              user_ylim = ylim) {
  ggplot2::ggproto(
    NULL, proto,
    limits = list(x = xlim, y = ylim), reverse = "none",
    expand = expand, default = FALSE, clip = clip, ratio = ratio,
    rotation = rotation, fit = fit,
    user_xlim = user_xlim, user_ylim = user_ylim
  )
}

#' Construct a CoordChord instance
#' @noRd
new_coord_chord <- function(rotation, ratio, xlim, ylim, expand, clip, fit,
                            user_xlim = xlim, user_ylim = ylim) {
  new_ggchord_coord(
    CoordChord, rotation, ratio, xlim, ylim, expand, clip, fit,
    user_xlim, user_ylim
  )
}

#' Project reference sequences onto arcs
#'
#' @inheritParams coord_dotplot
#' @param start_angle Initial angle in degrees counterclockwise from the right.
#' @param direction Direction of the positive reference frame.
#' @return A coordinate system for a ggchord plot.
#' @export
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
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' ggchord(data = seqs) +
#'   geom_seq() +
#'   geom_gene(data = genes, aes(fill = type)) +
#'   geom_link_ribbon(data = links) +
#'   coord_chord(start_angle = 90)
coord_chord <- function(start_angle = 0, direction = "counterclockwise",
    xlim = NULL, ylim = NULL, expand = FALSE, clip = "on", fit = "geometry") {
  ggchord_check_argument_names(sys.call(), sys.function())
  direction <- ggchord_coord_direction(start_angle, direction)
  settings <- ggchord_coord_settings(1, xlim, ylim, expand, clip, fit)
  ggplot2::ggproto(NULL, CoordChord, !!!settings,
    rotation = start_angle %% 360, reference_direction = direction)
}

ggchord_coord_direction <- function(start_angle, direction) {
  if (!is.numeric(start_angle) || length(start_angle) != 1L || !is.finite(start_angle))
    ggchord_stop("start_angle must be one finite angle in degrees")
  ggchord_match_choice(direction, c("clockwise", "counterclockwise"), "direction")
}

ggchord_match_choice <- function(value, choices, field) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !value %in% choices)
    ggchord_stop(field, " must be one of ", paste(choices, collapse = ", "))
  value
}

ggchord_coord_settings <- function(ratio, xlim, ylim, expand, clip, fit) {
  if (!is.null(ratio) && (!is.numeric(ratio) || length(ratio) != 1L ||
      !is.finite(ratio) || ratio <= 0))
    ggchord_stop("ratio must be NULL or a positive finite number")
  for (field in c("xlim", "ylim")) {
    limit <- get(field)
    if (!is.null(limit) && (!is.numeric(limit) || length(limit) != 2L ||
        any(!is.finite(limit)) || limit[1L] >= limit[2L]))
      ggchord_stop(field, " must be NULL or two increasing finite limits")
  }
  if (!is.logical(expand) || length(expand) != 1L || is.na(expand))
    ggchord_stop("expand must be TRUE or FALSE")
  clip <- ggchord_match_choice(clip, c("on", "off"), "clip")
  fit <- ggchord_match_choice(fit, c("geometry", "labels", "manual"), "fit")
  if (fit == "manual" && (is.null(xlim) || is.null(ylim)))
    ggchord_stop("fit = 'manual' requires both xlim and ylim")
  expand <- if (expand) c(top = is.null(ylim), right = is.null(xlim),
    bottom = is.null(ylim), left = is.null(xlim)) else FALSE
  list(ratio = ratio, limits = list(x = xlim, y = ylim),
    expand = expand, clip = clip, fit = fit, default = FALSE, reverse = "none",
    user_xlim = xlim, user_ylim = ylim)
}
