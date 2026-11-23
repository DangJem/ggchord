#' Attach a track value guide to a reference path
#'
#' @param position Display start or end of the reference path.
#' @param sequences Source accessions, or a list of source and instance IDs.
#' @param angle Tick label angle or the default waiver.
#' @param check.overlap Hide overlapping tick labels.
#' @param theme Optional local axis theme.
#' @return A guide for [`scale_y_continuous()`][ggplot2::scale_y_continuous].
#' @details For ordinary track value axes, use [ggplot2::guide_axis()] with
#'   `theme`, `angle`, `check.overlap`, `minor.ticks`, `cap`, and an optional
#'   text `title` (`NULL` removes it). The track uses its own reference path,
#'   so native panel-side `position`, non-default `n.dodge` and `order` are
#'   unsupported. Use this wrapper to choose the path endpoint or references.
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
#' ggchord(data = seqs) + layout_track("signal") +
#'   geom_track_line(data = signal, aes(x = pos, y = value), track = "signal") +
#'   scale_y_continuous(guide = guide_track_axis(
#'     position = "end", sequences = "A", angle = 0)) +
#'   coord_collinear()
guide_track_axis <- function(position = "start", sequences = NULL,
    angle = waiver(), check.overlap = FALSE, theme = NULL) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  position <- ggchord_match_choice(position, c("start", "end"), "position")
  ggchord_check_axis_selection(sequences)
  ggchord_check_axis_theme(theme, check.overlap)
  if (!is_waiver(angle) &&
      (!is.numeric(angle) || length(angle) != 1L || !is.finite(angle)))
    ggchord_stop("angle must be waiver() or finite degrees")
  guide <- guide_axis(angle = angle, check.overlap = check.overlap, theme = theme)
  guide$ggchord_position <- position
  guide$ggchord_sequences <- sequences
  guide
}

# Track value scales own their own axis Guide declaration. Theme affects only
# decorations; the named scale and physical position retain their old domains.
ggchord_track_axis_guide <- function(guide) {
  if (is_waiver(guide)) guide <- "axis"
  if (is.character(guide)) return(list(show = ggchord_match_arg(guide, c("axis", "none")) == "axis", guide = NULL))
  if (!inherits(guide, "GuideAxis"))
    ggchord_stop("Track value guide must be 'axis', 'none', guide_axis(), or guide_track_axis()")
  params <- guide$params
  defaults <- guide_axis()$params
  supported <- c("theme", "angle", "check.overlap", "title", "minor.ticks", "cap")
  if (!is_waiver(params$position) && !is.null(params$position))
    ggchord_stop("Track guide_axis(): 'position' is a panel side; use ",
      "guide_track_axis(position = 'start') or guide_track_axis(position = 'end')")
  for (name in setdiff(names(params), c(supported, "position"))) {
    if (!isTRUE(all.equal(params[[name]], defaults[[name]])))
      ggchord_stop("Track guide_axis(): unsupported option '", name,
        "'; use theme, angle, check.overlap, title, minor.ticks or cap")
  }
  title <- params$title
  if (!is_waiver(title) && !is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title)))
    ggchord_stop("Track guide_axis(): title must be one character string or NULL")
  local <- params$theme
  allowed <- c("axis.line", "axis.ticks", "axis.ticks.length",
    "axis.minor.ticks.length", "axis.text", "axis.title")
  if (!is.null(local) && any(!names(local) %in% allowed)) ggchord_stop(
    "Track guide_axis(): local theme supports ", paste(allowed, collapse = ", "))
  list(show = TRUE, guide = guide)
}

ggchord_track_axis_style <- function(plot, id) {
  sc <- ggchord_track_y_scale(plot, id)
  guide <- sc$ggchord_track_axis_guide
  local <- guide$params$theme
  themed <- plot
  if (!is.null(local)) {
    decorations <- as.list(local)[setdiff(names(local), "axis.minor.ticks.length")]
    overrides <- setNames(decorations, paste0("ggchord.track.", names(decorations)))
    themed$theme <- plot$theme + do.call(theme, overrides)
  }
  element <- function(role) ggchord_theme_element(themed, paste0("ggchord.track.axis.", role))
  line <- element("line"); text <- element("text"); title <- element("title"); tick <- element("ticks")
  angle <- guide$params$angle
  if (!is.null(angle) && !is_waiver(angle)) {
    if (!is.numeric(angle) || length(angle) != 1L || !is.finite(angle)) ggchord_stop("Track guide_axis(): angle must be finite")
    if (!inherits(text, "element_blank")) text@angle <- angle
  }
  as_params <- function(el) {
    if (inherits(el, "element_blank")) return(list(hidden = TRUE))
    if (inherits(el, "element_line")) list(colour = el@colour, linewidth = el@linewidth, linetype = el@linetype,
      lineend = el@lineend, linejoin = el@linejoin) else list(colour = el@colour,
        size = el@size / ggplot2::.pt, family = el@family, fontface = el@face,
        angle = el@angle %||% 0, hjust = el@hjust %||% .5, vjust = el@vjust %||% .5,
        lineheight = el@lineheight %||% 1.2)
  }
  offset <- function(name, fallback) {
    key <- paste0("ggchord.track.axis.", name, ".offset")
    value <- plot$theme[[key]] %||% unit(fallback, "mm")
    ggchord_check_axis_unit(value, key, nonnegative = TRUE)
    ggchord_unit_inches(value) * 25.4
  }
  tick_length <- themed$theme[["ggchord.track.axis.ticks.length"]] %||% unit(1.2, "mm")
  ggchord_check_axis_unit(tick_length, "Track axis major tick length")
  minor <- NULL
  if (isTRUE(guide$params$minor.ticks)) {
    minor <- local[["axis.minor.ticks.length"]] %||% plot$theme$axis.minor.ticks.length %||% rel(.75)
    minor <- if (inherits(minor, "rel")) as.numeric(minor) * tick_length else minor
    ggchord_check_axis_unit(minor, "Track axis minor tick length")
  }
  list(axis_line_params = as_params(line), axis_tick_params = as_params(tick),
    tick_length_mm = ggchord_unit_inches(tick_length) * 25.4,
    minor_tick_length_mm = if (!is.null(minor)) ggchord_unit_inches(minor) * 25.4,
    axis_text_params = as_params(text),
    axis_title_params = as_params(title), axis_check_overlap = isTRUE(guide$params$check.overlap),
    text_offset_mm = offset("text",2.5), title_offset_mm = offset("title",5.5),
    text_offset_locked = !is.null(plot$ggchord$theme_explicit[["ggchord.track.axis.text.offset"]]),
    title_offset_locked = !is.null(plot$ggchord$theme_explicit[["ggchord.track.axis.title.offset"]]))
}

# Physical half-extent along the actual screen-space offset direction.
ggchord_track_axis_extent <- function(labels, style, tangent) {
  if (!length(labels) || isTRUE(style$hidden)) return(0)
  d <- data.frame(text_x = 0, text_y = 0, text = as.character(labels))
  for (name in c("size", "family", "fontface", "lineheight", "hjust", "vjust"))
    if (length(style[[name]])) d[[name]] <- style[[name]]
  d$text_angle <- style$angle %||% 0
  boxes <- ggchord_text_boxes(d, units_per_inch = 25.4)
  max(abs(c(tangent[1] * boxes$xmin + tangent[2] * boxes$ymin,
    tangent[1] * boxes$xmin + tangent[2] * boxes$ymax,
    tangent[1] * boxes$xmax + tangent[2] * boxes$ymin,
    tangent[1] * boxes$xmax + tangent[2] * boxes$ymax)))
}
