# Track value scales own their own axis Guide declaration. Theme affects only
# decorations; the named scale and physical position retain their old domains.
ggchord_track_axis_guide <- function(guide) {
  if (is.character(guide)) return(list(show = match.arg(guide, c("axis", "none")) == "axis", guide = NULL))
  if (!inherits(guide, "GuideAxis"))
    ggchord_stop("Track value guide must be 'axis', 'none', or guide_axis(theme = ...)")
  params <- guide$params; defaults <- ggplot2::guide_axis()$params
  supported <- c("theme", "angle", "check.overlap")
  for (name in setdiff(names(params), supported))
    if (!identical(params[[name]], defaults[[name]])) ggchord_stop(
      "Track guide_axis(): unsupported option '", name, "'; use theme, angle or check.overlap")
  local <- params$theme
  allowed <- c("axis.line", "axis.text", "axis.title")
  if (!is.null(local) && any(!names(local) %in% allowed)) ggchord_stop(
    "Track guide_axis(): local theme supports axis.line, axis.text and axis.title")
  list(show = TRUE, guide = guide)
}

ggchord_track_axis_style <- function(plot, id) {
  sc <- plot$scales$get_scales(ggchord_track_aes("value", id))
  guide <- sc$ggchord_track_axis_guide
  local <- guide$params$theme
  themed <- plot
  if (!is.null(local)) {
    overrides <- stats::setNames(as.list(local), paste0("ggchord.track.", names(local)))
    themed$theme <- plot$theme + do.call(ggplot2::theme, overrides)
  }
  element <- function(role) ggchord_theme_element(themed, paste0("ggchord.track.axis.", role))
  line <- element("line"); text <- element("text"); title <- element("title")
  angle <- guide$params$angle
  if (!is.null(angle) && !ggplot2::is_waiver(angle)) {
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
    value <- plot$theme[[paste0("ggchord.track.axis.", name, ".offset")]] %||% grid::unit(fallback,"mm")
    if (!grid::is.unit(value) || length(value) != 1L || !is.finite(as.numeric(value)) || as.numeric(value) < 0)
      ggchord_stop("Track axis offsets require nonnegative scalar absolute units")
    ggchord_unit_inches(value) * 25.4
  }
  list(axis_line_params = as_params(line), axis_text_params = as_params(text),
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
