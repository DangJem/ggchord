#' Specify feature geometry
#'
#' @param name `"prototype"` or a preset supported by [geom_feature()] or
#'   [geom_site()].
#' @param head_length,tail_length Requested end length along the sequence in
#'   millimetres or absolute grid units. The full-span triangle uses the
#'   interval length.
#' @param head_width,tail_width End width across the sequence in millimetres
#'   or absolute grid units. A notch and the head of `"arrow_flush"` use
#'   the body width. For a circle or diamond, `head_width` sets the
#'   nominal symbol diameter or cross-sequence width.
#' @param short_feature End treatment when the visible interval is too short.
#'   Choose `"shrink"` (default), `"block"`, or `"error"`.
#' @param body_width Body width in millimetres or absolute grid units. For a
#'   block at a site or cut boundary it sets the side of the square. For a
#'   lollipop it sets stem thickness. For tick and cut it sets the full
#'   normal extent, with each cut stroke using half.
#' @param body_fill,body_colour,head_fill,head_colour,tail_fill,tail_colour
#'   Fill and outline of each part. An explicit layer value or mapped
#'   aesthetic takes priority for the corresponding channel. `NA` suppresses
#'   the specified channel.
#' @param body_shape For the prototype, `"block"` (default), `"line"` for a
#'   thin body, `"none"` with a full-span triangle or a single head, or
#'   `"tick"` and `"cut"` for normal strokes.
#' @param head_shape,tail_shape For the prototype, `"none"` (flat),
#'   `"triangle"`, `"cap"`, or `"notch"`. A head can also be
#'   `"circle"` or `"diamond"` for anchored marks.
#' @param cap_length,cap_width Common length and width for ends that use
#'   `"cap"`. A dimension supplied for a particular end takes priority.
#' @return A feature shape specification. Supply it to the layer's `shape`
#'   argument or as a value in [ggplot2::scale_shape_manual()].
#' @details The prototype starts as a block with two flat ends. A preset
#'   selects a predefined outline and accepts dimensions and styles for its
#'   active parts. Use the prototype to change the outline. The visible
#'   interval sets the total length. The body fills the space between both
#'   ends. The strand determines where the head and tail go on directional
#'   shapes. Legend keys show a representative sample of the same profile.
#'
#' `"arrow"` adds a shouldered triangular head. `"arrow_flush"` uses a
#' head as wide as its body. `"arrow_triangle"` is a full-span triangle.
#' `"arrow_arc"` pairs a narrow body with a triangular head. `"chevron"`
#' adds a triangular head and notched tail. `"capped_line"` has a narrow body
#' with a cap at each end. Settings for absent parts raise an error.
#'
#' The same block specification spans an interval, draws a square centred on
#' a site, or draws a square on the corresponding side of each known cut
#' boundary. With no known cut boundary, the block is omitted. When no
#' spec width is supplied, the layer's thickness or size supplies the width.
#' Block legend keys show a representative location, including two squares
#' for a staggered cut.
#'
#' Use a circular head with a line body for `"lollipop"`, or with no
#' body for `"point"`. A diamond head with no body gives `"diamond"`.
#' Use a tick or cut body with flat ends for `"tick"` or `"cut"`.
#' The same lollipop specification starts from the sequence backbone at a
#' visible interval midpoint, from the recorded `pos` at a site, or
#' from each recorded chain boundary at a cut. On sites and cuts the
#' stem extends 1.5 times the layer size. Interval stems reach the feature
#' track. The layer controls placement in each location.
#' @examples
#' library(ggplot2)
#' seqs <- data.frame(accver = "A", length = 1000, sequence_unit = "bp")
#' intervals <- data.frame(
#'   accver = "A", start = c(100, 580), end = c(430, 900),
#'   strand = "+", kind = c("short", "standard")
#' )
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_feature(aes(shape = kind), data = intervals,
#'     position = position_offset(4), fill = "#77B5AE") +
#'   scale_shape_manual(name = "Arrow geometry", values = c(
#'     short = shape_feature("arrow", head_length = 1),
#'     standard = "arrow"
#'   )) +
#'   coord_collinear()
#'
#' block_reference <- data.frame(
#'   accver = "A", length = 1000L, topology = "circular", sequence_unit = "bp"
#' )
#' block_intervals <- data.frame(accver = "A", start = 120L, end = 310L, kind = "Interval")
#' block_sites <- data.frame(accver = "A", pos = 490L, kind = "Site")
#' block_cuts <- data.frame(
#'   accver = "A", pos = 700L, cut_top = 700L, cut_bottom = 810L, kind = "Cut"
#' )
#' prototype <- shape_feature("prototype",
#'   body_width = 4, body_fill = "#7ABCB7", body_colour = "#235354"
#' )
#'
#' block_plot <- ggchord(block_reference) +
#'   geom_seq(geometry = "line", colour = "#89979A", show.legend = FALSE) +
#'   geom_feature(aes(shape = kind), data = block_intervals) +
#'   geom_feature(aes(shape = kind), data = block_sites) +
#'   geom_feature(aes(shape = kind), data = block_cuts) +
#'   scale_shape_manual(
#'     name = "Location type",
#'     values = c(Interval = prototype, Site = prototype, Cut = prototype),
#'     breaks = c("Interval", "Site", "Cut"),
#'     guide = guide_legend(nrow = 1)
#'   ) +
#'   scale_seq_continuous(guide = "none") +
#'   theme_ggchord_minimal() +
#'   theme(
#'     legend.position = "bottom",
#'     legend.key.width = grid::unit(10, "mm"),
#'     legend.key.height = grid::unit(8, "mm"),
#'     plot.margin = margin(12, 14, 12, 14)
#'   )
#'
#' block_plot + coord_collinear() +
#'   labs(title = "One block specification in collinear coordinates")
#' block_plot + coord_circular(fit = "labels") +
#'   labs(title = "One block specification in circular coordinates")
#'
#' double_cap <- shape_feature("prototype",
#'   body_shape = "line", head_shape = "cap", tail_shape = "cap",
#'   body_width = 0.4, cap_length = 1, cap_width = 2, tail_width = 3,
#'   head_fill = "#D98245", tail_fill = "#5885B2"
#' )
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_feature(data = intervals, shape = double_cap,
#'     position = position_offset(4)) +
#'   coord_collinear()
#'
#' lollipop <- shape_feature("prototype",
#'   body_shape = "line", head_shape = "circle",
#'   body_width = 0.4, body_colour = "#274F80",
#'   head_width = 3, head_fill = "#EF9F55"
#' )
#' ggchord(block_reference) +
#'   geom_seq() +
#'   geom_feature(data = block_sites, shape = lollipop) +
#'   geom_feature(data = block_cuts, shape = lollipop) +
#'   coord_circular()
#' @export
shape_feature <- function(name, head_length = NULL, head_width = NULL,
    short_feature = NULL, body_width = NULL, body_fill = NULL, body_colour = NULL,
    body_shape = NULL, head_shape = NULL, tail_shape = NULL,
    tail_length = NULL, tail_width = NULL, cap_length = NULL, cap_width = NULL,
    head_fill = NULL, head_colour = NULL, tail_fill = NULL, tail_colour = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  name <- ggchord_match_choice(name, unique(c("prototype", ggchord_feature_shapes,
    ggchord_site_shapes, "cut")), "shape_feature(): name")
  options <- list(body_shape = body_shape, head_shape = head_shape, tail_shape = tail_shape,
    short_feature = short_feature, body_width = body_width, head_length = head_length,
    head_width = head_width, tail_length = tail_length, tail_width = tail_width,
    cap_length = cap_length, cap_width = cap_width, body_fill = body_fill,
    body_colour = body_colour, head_fill = head_fill, head_colour = head_colour,
    tail_fill = tail_fill, tail_colour = tail_colour)
  options <- Filter(Negate(is.null), options)
  if (!length(options)) return(if (name == "prototype") "block" else name)
  sizes <- intersect(names(options), c("body_width", "head_length", "head_width",
    "tail_length", "tail_width", "cap_length", "cap_width"))
  for (field in sizes)
    options[[field]] <- ggchord_feature_shape_mm(options[[field]], field)
  paints <- grep("_(fill|colour)$", names(options), value = TRUE)
  options[paints] <- lapply(options[paints], ggchord_encode_body_colour)
  ggchord_resolve_feature_profile(name, options)
  paste("ggchord", "parts", name, paste(names(options), options, sep = "=", collapse = ";"),
    sep = ":")
}

ggchord_feature_parts_fields <- function() {
  c("body_shape", "head_shape", "tail_shape", "short_feature", "body_width",
    "head_length", "head_width", "tail_length", "tail_width", "cap_length",
    "cap_width", "body_fill", "body_colour", "head_fill", "head_colour",
    "tail_fill", "tail_colour")
}

ggchord_midpoint_shapes <- c("diamond", "lollipop", "point", "tick")

ggchord_resolve_feature_profile <- function(name, options = list()) {
  profile <- switch(name,
    prototype = list(body_shape = "block", head_shape = "none", tail_shape = "none"),
    block = list(body_shape = "block", head_shape = "none", tail_shape = "none"),
    arrow = list(body_shape = "block", head_shape = "triangle", tail_shape = "none"),
    arrow_flush = list(body_shape = "block", head_shape = "triangle", tail_shape = "none"),
    arrow_triangle = list(body_shape = "none", head_shape = "triangle", tail_shape = "none"),
    arrow_arc = list(body_shape = "line", head_shape = "triangle", tail_shape = "none"),
    chevron = list(body_shape = "block", head_shape = "triangle", tail_shape = "notch"),
    capped_line = list(body_shape = "line", head_shape = "cap", tail_shape = "cap"),
    diamond = list(body_shape = "none", head_shape = "diamond", tail_shape = "none"),
    point = list(body_shape = "none", head_shape = "circle", tail_shape = "none"),
    lollipop = list(body_shape = "line", head_shape = "circle", tail_shape = "none"),
    tick = list(body_shape = "tick", head_shape = "none", tail_shape = "none"),
    cut = list(body_shape = "cut", head_shape = "none", tail_shape = "none"),
    NULL)
  if (is.null(profile)) {
    if (length(options)) ggchord_stop("shape_feature(): these settings require an interval profile")
    return(list(name = name))
  }
  topology <- intersect(c("body_shape", "head_shape", "tail_shape"), names(options))
  if (length(topology) && name != "prototype")
    ggchord_stop("shape_feature(): use 'prototype' to configure the body, head or tail shape")
  profile[names(options)] <- options
  profile$body_shape <- ggchord_match_choice(profile$body_shape, c("block", "line", "none", "tick", "cut"),
    "shape_feature(): body_shape")
  ends <- c("none", "triangle", "cap", "notch")
  profile$head_shape <- ggchord_match_choice(profile$head_shape, c(ends, "circle", "diamond"),
    "shape_feature(): head_shape")
  profile$tail_shape <- ggchord_match_choice(profile$tail_shape, ends, "shape_feature(): tail_shape")
  profile$short_feature <- ggchord_match_choice(profile$short_feature %||% "shrink",
    c("shrink", "block", "error"), "shape_feature(): short_feature")
  if (name == "prototype") {
    name <- switch(paste(profile$body_shape, profile$head_shape, profile$tail_shape),
      "block none none" = "block", "line none none" = "prototype",
      "block triangle none" = "arrow", "line triangle none" = "arrow_arc",
      "none triangle none" = "arrow_triangle", "block triangle notch" = "chevron",
      "line cap cap" = "capped_line", "none diamond none" = "diamond",
      "none circle none" = "point", "line circle none" = "lollipop",
      "tick none none" = "tick", "cut none none" = "cut", "prototype")
  }
  anchored <- name %in% c(ggchord_midpoint_shapes, "cut")
  if (!anchored && (profile$body_shape %in% c("tick", "cut") ||
      profile$head_shape %in% c("circle", "diamond") ||
      (profile$body_shape == "none" &&
        !(profile$head_shape == "triangle" && profile$tail_shape == "none"))))
    ggchord_stop("shape_feature(): use an admissible body, head and tail combination")
  allowed <- c(topology)
  if (profile$body_shape != "none")
    allowed <- c(allowed, "body_width", "body_fill", "body_colour")
  if (profile$head_shape != "none")
    allowed <- c(allowed, "head_length", "head_width", "head_fill", "head_colour")
  if (profile$tail_shape != "none")
    allowed <- c(allowed, "tail_length", "tail_width", "tail_fill", "tail_colour")
  if (profile$head_shape == "cap" || profile$tail_shape == "cap")
    allowed <- c(allowed, "cap_length", "cap_width")
  if (profile$body_shape != "none" &&
      (profile$head_shape != "none" || profile$tail_shape != "none"))
    allowed <- c(allowed, "short_feature")
  if (name == "arrow_flush") allowed <- setdiff(allowed, "head_width")
  if (profile$body_shape == "none") allowed <- setdiff(allowed, "head_length")
  if (profile$head_shape == "notch") allowed <- setdiff(allowed, "head_width")
  if (profile$tail_shape == "notch") allowed <- setdiff(allowed, "tail_width")
  if (anchored) {
    allowed <- c(topology, switch(name,
      diamond = , point = c("head_width", "head_fill", "head_colour"),
      lollipop = c("body_width", "body_colour", "head_width", "head_fill", "head_colour"),
      tick = , cut = c("body_width", "body_colour")))
  }
  invalid <- setdiff(names(options), allowed)
  if (length(invalid))
    ggchord_stop("shape_feature(): ", invalid[1L], " does not apply to this profile")
  for (end in c("head", "tail")) {
    if (profile[[paste0(end, "_shape")]] != "cap") next
    for (dimension in c("length", "width")) {
      field <- paste0(end, "_", dimension)
      profile[[field]] <- profile[[field]] %||% profile[[paste0("cap_", dimension)]]
    }
  }
  for (field in grep("_(fill|colour)$", names(profile), value = TRUE)) {
    if (identical(profile[[field]], "NA")) profile[[field]] <- NA_character_
  }
  # Use the existing whole-outline draw for unmodified presets. A configured
  # end, component size or part colour uses the common component geometry.
  active <- setdiff(names(options), topology)
  changes_parts <- if (name %in% c("arrow", "arrow_flush"))
    setdiff(active, c("head_length", "head_width", "short_feature")) else active
  profile$components <- !anchored && (identical(name, "prototype") ||
    (name != "block" && length(changes_parts) > 0L))
  profile$name <- name
  profile
}

ggchord_encode_body_colour <- function(value) {
  if (!is.character(value) || length(value) != 1L)
    ggchord_stop("shape_feature(): colours must be single colour strings or NA")
  if (is.na(value)) return("NA")
  rgb <- grDevices::col2rgb(value, alpha = TRUE)
  grDevices::rgb(rgb[1L], rgb[2L], rgb[3L], rgb[4L], maxColorValue = 255)
}

ggchord_feature_shape_mm <- function(value, name) {
  if (is.null(value)) return(NULL)
  mm <- ggchord_physical_mm(value, name, TRUE)
  if (!is.numeric(mm) || length(mm) != 1L || mm <= 0)
    ggchord_stop("shape_feature(): ", name, " must be a positive physical size")
  mm
}

ggchord_parse_feature_shape <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value))
    ggchord_stop("shape must be a preset name or a shape_feature() specification")
  if (!startsWith(value, "ggchord:")) {
    if (identical(value, "primer_arc"))
      ggchord_stop("Use shape = 'arrow_arc' in place of 'primer_arc'")
    if (value %in% c("compact_arrow", "promoter_arrow", "primer_arrow"))
      ggchord_stop("Use shape_feature('arrow', head_length = 1) for a short-headed arrow")
    return(ggchord_resolve_feature_profile(value))
  }
  fields <- strsplit(value, ":", fixed = TRUE)[[1L]]
  if (length(fields) == 4L && identical(fields[2L], "parts")) {
    pairs <- strsplit(strsplit(fields[4L], ";", fixed = TRUE)[[1L]], "=", fixed = TRUE)
    if (any(lengths(pairs) != 2L)) ggchord_stop("Invalid shape_feature() specification")
    keys <- vapply(pairs, `[[`, character(1), 1L)
    if (anyDuplicated(keys) || any(!keys %in% ggchord_feature_parts_fields()))
      ggchord_stop("Invalid shape_feature() specification")
    options <- setNames(lapply(pairs, `[[`, 2L), keys)
    sizes <- grep("_(length|width)$", keys, value = TRUE)
    for (field in sizes) options[[field]] <- ggchord_feature_shape_mm(
      suppressWarnings(as.numeric(options[[field]])), field)
    paints <- grep("_(fill|colour)$", keys, value = TRUE)
    for (field in paints)
      if (options[[field]] != "NA") options[[field]] <- ggchord_encode_body_colour(options[[field]])
    result <- ggchord_resolve_feature_profile(fields[3L], options)
  } else if (length(fields) == 6L && fields[6L] %in% c("shape", "block")) {
    # Specifications captured by an earlier development snapshot remain
    # usable in bundled data and saved plots.
    legacy <- if (fields[6L] == "block") c("body_width", "body_fill", "body_colour") else
      c("head_length", "head_width", "short_feature")
    options <- as.list(setNames(fields[3:5], legacy))
    options <- Filter(nzchar, options)
    sizes <- grep("_(length|width)$", names(options), value = TRUE)
    for (field in sizes) options[[field]] <- ggchord_feature_shape_mm(
      suppressWarnings(as.numeric(options[[field]])), field)
    result <- ggchord_resolve_feature_profile(fields[2L], options)
  } else {
    ggchord_stop("Invalid shape_feature() specification")
  }
  result$value <- value
  result
}

ggchord_shape_body_styles <- function(data, specifications, explicit) {
  ggchord_shape_part_styles(data, specifications, "body", explicit)
}

ggchord_shape_part_styles <- function(data, specifications, part, explicit) {
  if (!length(specifications)) return(data)
  values <- unique(specifications[which(startsWith(specifications, "ggchord:"))])
  roles <- rep_len(part, length(specifications))
  for (value in values) {
    spec <- ggchord_parse_feature_shape(value)
    for (channel in setdiff(c("fill", "colour"), explicit)) {
      for (role in unique(roles[specifications %in% value])) {
        colour <- spec[[paste(role, channel, sep = "_")]]
        rows <- which(specifications == value & roles == role)
        if (!is.null(colour)) data[[channel]][rows] <- colour
      }
    }
  }
  data
}

ggchord_reject_shape_dots <- function(dots) {
  if ("head_style" %in% names(dots))
    ggchord_stop("head_style was replaced by shape = 'arrow', 'arrow_flush', or 'arrow_triangle'")
  moved <- intersect(c("head_length", "head_width", "short_feature"), names(dots))
  if (length(moved))
    ggchord_stop(moved[1L], " now belongs in shape = shape_feature(...)")
}
