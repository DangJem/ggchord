# Role-specific scales for ggchord. Keeping roles separate lets one plot use
# independent sequence, ribbon, gene, feature and region scales without the
# usual fill/colour collisions.

#' @noRd
scale_seq_colour_manual <- function(..., values) {
  scale_colour_manual(..., values = values, aesthetics = "seq_colour")
}

#' @noRd
scale_seq_color_manual <- scale_seq_colour_manual

#' @noRd
scale_seq_ring_manual <- function(..., values, name = "Ring",
                                  limits = NULL, guide = "none") {
  if (!is.numeric(values) || length(values) == 0L || anyNA(values) ||
      any(!is.finite(values)) || any(values <= 0)) {
    ggchord_stop(
      "scale_seq_ring_manual(): values must be finite positive radii"
    )
  }
  if (is.null(limits) && !is.null(names(values))) limits <- names(values)
  discrete_scale(
    aesthetics = "seq_ring", palette = scales::manual_pal(values),
    ..., name = name, limits = limits, na.value = NA_real_, guide = guide
  )
}

#' @noRd
scale_ribbon_fill_stepsn <- function(
    ..., colours, values = NULL, colors,
    guide = ggplot2::guide_coloursteps(available_aes = "ribbon_fill")) {
  has_colours <- !missing(colours)
  has_colors <- !missing(colors)
  if (has_colours && has_colors)
    ggchord_stop("scale_ribbon_fill_stepsn(): use only one of colours and colors")
  if (!has_colours && !has_colors)
    ggchord_stop("scale_ribbon_fill_stepsn(): colours (or colors) is required")
  if (!has_colours) colours <- colors
  ggplot2::scale_fill_stepsn(
    ..., colours = colours, values = values, aesthetics = "ribbon_fill",
    guide = guide
  )
}

#' @noRd
scale_ribbon_fill_gradientn <- function(
    ..., colours, values = NULL, colors,
    guide = ggplot2::guide_colourbar(available_aes = "ribbon_fill")) {
  has_colours <- !missing(colours)
  has_colors <- !missing(colors)
  if (has_colours && has_colors)
    ggchord_stop("scale_ribbon_fill_gradientn(): use only one of colours and colors")
  if (!has_colours && !has_colors)
    ggchord_stop("scale_ribbon_fill_gradientn(): colours (or colors) is required")
  if (!has_colours) colours <- colors
  scale_fill_gradientn(
    ..., colours = colours, values = values, aesthetics = "ribbon_fill",
    guide = guide
  )
}

#' @noRd
scale_ribbon_fill_manual <- function(..., values) {
  scale_fill_manual(..., values = values, aesthetics = "ribbon_fill")
}

#' @noRd
scale_ribbon_fill_identity <- function(...) {
  scale_fill_identity(..., aesthetics = "ribbon_fill")
}

#' @noRd
scale_ribbon_alpha_continuous <- function(..., range = NULL) {
  # Older 4.0.x finds the default palette by aesthetic name, so an extension
  # role needs the same explicit fallback range as the native alpha scale.
  if (is.null(range) && !"fallback.palette" %in% names(formals(continuous_scale))) {
    range <- c(0.1, 1)
  }
  ggplot2::scale_alpha_continuous(
    ..., range = range, aesthetics = "ribbon_alpha"
  )
}

#' @noRd
scale_ribbon_alpha_manual <- function(..., values) {
  ggplot2::scale_alpha_manual(..., values = values, aesthetics = "ribbon_alpha")
}

#' @noRd
scale_ribbon_colour_manual <- function(..., values) {
  scale_colour_manual(
    ..., values = values, aesthetics = "ribbon_colour"
  )
}

#' @noRd
scale_ribbon_color_manual <- scale_ribbon_colour_manual

#' @noRd
scale_ribbon_linetype_manual <- function(..., values) {
  ggplot2::scale_linetype_manual(
    ..., values = values, aesthetics = "ribbon_linetype"
  )
}

#' @noRd
scale_gene_fill_manual <- function(..., values) {
  scale_fill_manual(..., values = values, aesthetics = "gene_fill")
}

#' @noRd
scale_feature_fill_manual <- function(..., values, limits = NULL) {
  if (is.null(limits) && !is.null(names(values))) limits <- names(values)
  out <- scale_fill_manual(..., values = values, limits = limits,
    aesthetics = "feature_fill")
  attr(out, "ggchord_scale_priority") <- 2L
  out
}

#' @noRd
scale_feature_shape_manual <- function(..., values, name = "Feature",
                                       limits = NULL, guide = "none") {
  allowed <- c(
    "arrow", "compact_arrow", "promoter_arrow", "primer_arrow", "marker",
    "block", "capped_line", "primer_arc", "chevron", "lollipop"
  )
  if (!is.character(values) || length(values) == 0L || anyNA(values) ||
      any(!values %in% allowed)) {
    ggchord_stop(
      "scale_feature_shape_manual(): values contain an unknown geometry"
    )
  }
  if (is.null(limits) && !is.null(names(values))) limits <- names(values)
  out <- discrete_scale(
    aesthetics = "feature_shape",
    palette = scales::manual_pal(values),
    ..., name = name, limits = limits, na.value = "arrow", guide = guide
  )
  attr(out, "ggchord_scale_priority") <- 2L
  out
}

ggchord_plasmid_feature_colours <- function() {
  c(
    CDS="#993366", gene="#993366", resistance_gene="#CCFFCC",
    selection_marker="#CCFFCC", reporter="#05FD14",
    peptide="#CC99B2", tag="#CC99B2",
    promoter="#FFFFFF", rep_origin="#FFFF00", replication_origin="#FFFF00",
    ori="#FFFF00", origin="#FFFF00", primer_bind="#A020F0",
    primer="#A020F0", terminator="#FFFFFF", enhancer="#FFFFFF",
    protein_bind="#31849B", binding_site="#31849B", operator="#31849B",
    RBS="#A6ACB3", polyA_signal="#A6ACB3", poly_a_signal="#A6ACB3",
    regulatory="#A6ACB3", LTR="#FFE4C4", repeat_region="#FFE4C4",
    misc_RNA="#00CCFF", MCS="#99CCFF", misc_feature="#A6ACB3"
  )
}

#' Lighten a feature colour for an external callout
#' @noRd
ggchord_feature_callout_fill <- function(fill, luminance = .56) {
  vapply(as.character(fill), function(value) {
    rgb <- tryCatch(
      col2rgb(value, alpha = TRUE)[, 1L] / 255,
      error = function(e) c(184, 189, 195, 255) / 255
    )
    brightness <- function(colour) {
      linear <- ifelse(colour <= .04045, colour / 12.92,
        ((colour + .055) / 1.055)^2.4)
      sum(linear * c(.2126, .7152, .0722))
    }
    # Preserve already light colours. Dark fills are mixed with just enough
    # white for readable dark callout text; one large fixed blend washed out
    # operators and reporter colours into nearly identical grey boxes.
    low <- 0
    high <- 1
    if (brightness(rgb[1:3]) < luminance) {
      for (i in seq_len(18L)) {
        amount <- (low + high) / 2
        mixed <- rgb[1:3] + (1 - rgb[1:3]) * amount
        if (brightness(mixed) < luminance) low <- amount else high <- amount
      }
    } else high <- 0
    mixed <- rgb[1:3] + (1 - rgb[1:3]) * high
    grDevices::rgb(mixed[1L], mixed[2L], mixed[3L], alpha = rgb[4L])
  }, character(1L), USE.NAMES = FALSE)
}

ggchord_contrast_colour <- function(fill) {
  vapply(as.character(fill), function(value) {
    rgb <- tryCatch(col2rgb(value)[, 1L] / 255,
      error = function(e) c(1, 1, 1))
    linear <- ifelse(rgb <= .04045, rgb / 12.92,
      ((rgb + .055) / 1.055)^2.4)
    if (sum(linear * c(.2126, .7152, .0722)) < .36) "#FFFFFF" else "#202020"
  }, character(1), USE.NAMES = FALSE)
}

#' @noRd
scale_feature_fill_plasmid <- function(..., limits = NULL) {
  values <- ggchord_plasmid_feature_colours()
  if (is.null(limits)) limits <- names(values)
  out <- scale_fill_manual(
    ..., values=values, limits=limits, aesthetics="feature_fill", na.value="#B8BDC3"
  )
  attr(out,"ggchord_scale_priority") <- 1L
  out
}

#' @noRd
scale_feature_shape_plasmid <- function(..., limits = NULL, guide = "none") {
  values <- c(
    CDS="arrow", gene="arrow", resistance_gene="arrow",
    selection_marker="arrow", reporter="arrow", peptide="compact_arrow",
    tag="compact_arrow",
    rep_origin="arrow", replication_origin="arrow", ori="arrow",
    promoter="promoter_arrow", primer_bind="primer_arrow",
    primer="primer_arrow", protein_bind="block", binding_site="block",
    operator="block", terminator="block", enhancer="block",
    regulatory="compact_arrow", RBS="block", polyA_signal="block",
    poly_a_signal="block", LTR="block", repeat_region="block",
    misc_RNA="block", intron="capped_line", MCS="block", misc_feature="block"
  )
  if (is.null(limits)) limits <- names(values)
  out <- discrete_scale(
    aesthetics="feature_shape",palette=scales::manual_pal(values),
    ...,limits=limits,na.value="block",guide=guide
  )
  attr(out,"ggchord_scale_priority") <- 1L
  out
}

#' @noRd
scale_region_fill_manual <- function(..., values) {
  scale_fill_manual(..., values = values, aesthetics = "region_fill")
}

#' @noRd
scale_seq_position_continuous <- function(
    name = ggplot2::waiver(), breaks = ggplot2::waiver(),
    minor_breaks = ggplot2::waiver(), labels = ggplot2::waiver(),
    limits = NULL, expand = ggplot2::waiver(), oob = scales::censor,
    transform = "identity", n.breaks = 6,
    coordinates = c("relative", "source"), guide = guide_seq_axis(), ...) {
  ggchord_check_argument_names(sys.call(), sys.function(), continuous_scale)
  coordinates <- ggchord_match_arg(coordinates)
  if ("sequences" %in% names(list(...)))
    ggchord_stop("Use guide = guide_seq_axis(sequences = ...) to select displayed axes")
  if (identical(guide, "none")) guide <- guide_seq_axis(character())
  if (!inherits(guide, "ggchord_seq_axis")) ggchord_stop("Use guide_seq_axis() or 'none'")
  sc <- continuous_scale(
    aesthetics = "seq_position",
    palette = function(x) x,
    name = name, breaks = breaks, minor_breaks = minor_breaks,
    labels = labels, limits = limits, expand = expand, oob = oob,
    transform = transform, n.breaks = n.breaks, guide = "none", ...
  )
  sc$ggchord_position_coordinates <- coordinates
  sc$ggchord_axis_sequences <- guide$sequences
  sc$ggchord_check_overlap <- guide$check.overlap
  sc
}

#' Configure a sequence position guide
#'
#' @param sequences Source accessions, or a list with source and instance IDs.
#' @param side Positive or negative local side of the reference path.
#' @param angle Automatic label orientation, or a numeric angle in degrees.
#' @param check.overlap Whether to adjust crowded sequence-end labels and hide
#'   overlapping interior labels. Endpoints keep their ticks and values. On
#'   projected axes, competing endpoint labels can move along their own
#'   sequence paths. Set to `FALSE` to keep the original positions.
#' @param theme Optional local axis theme.
#' @return A guide for [scale_seq_position()].
#' @export
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' ggchord(data = seqs) + geom_seq() +
#'   scale_seq_position(breaks = scales::breaks_width(200),
#'     guide = guide_seq_axis(sequences = "A", side = "negative",
#'       theme = theme(axis.text = element_text(size = 8)))) +
#'   coord_collinear()
guide_seq_axis <- function(sequences = NULL, side = "positive", angle = "auto",
    check.overlap = TRUE, theme = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_check_axis_selection(sequences)
  side <- ggchord_match_choice(side, c("positive", "negative"), "side")
  if (!identical(angle, "auto") && (!is.numeric(angle) || length(angle) != 1L || !is.finite(angle)))
    ggchord_stop("angle must be 'auto' or finite degrees")
  ggchord_check_axis_theme(theme, check.overlap)
  structure(list(sequences = sequences, side = side, angle = angle,
    check.overlap = check.overlap, theme = theme), class = "ggchord_seq_axis")
}

ggchord_check_axis_selection <- function(selection) {
  if (is.null(selection) || identical(selection, character())) return(invisible(NULL))
  if (is.list(selection)) {
    if (is.null(names(selection)) || anyDuplicated(names(selection)) ||
        any(!names(selection) %in% c("source", "instance")))
      ggchord_stop("Axis selection accepts list(source = ..., instance = ...)")
    for (ids in selection) ggchord_id_vector(ids, "Guide sequences")
  } else ggchord_id_vector(selection, "Guide sequences")
}

ggchord_axis_selection <- function(selection, sources, instances = NULL) {
  if (is.null(selection)) return(NULL)
  source_ids <- instances$accver %||% sources
  keys <- instances$key %||% sources
  source <- if (is.list(selection)) selection$source else selection
  local <- if (is.list(selection)) selection$instance
  if (any(!source %in% source_ids) || any(!local %in% (instances$instance_id %||% keys)))
    ggchord_stop("Axis guide selects an unknown source or instance")
  keys[source_ids %in% source | (instances$instance_id %||% keys) %in% local]
}

ggchord_check_axis_theme <- function(theme, overlap) {
  if (!is.logical(overlap) || length(overlap) != 1L || is.na(overlap))
    ggchord_stop("check.overlap must be TRUE or FALSE")
  if (!is.null(theme) && (!inherits(theme, "theme") ||
      any(!names(theme) %in% c("axis.line", "axis.ticks", "axis.ticks.length", "axis.text", "axis.title"))))
    ggchord_stop("The local guide theme can set axis.line, axis.ticks, axis.ticks.length, axis.text and axis.title")
}

#' Define the linear scale for sequence positions
#'
#' @param name Optional axis title.
#' @param breaks,minor_breaks Position vectors, functions or named lists.
#' @param labels Labels or a labelling function, optionally by source or instance.
#' @param coordinates Use local or original source coordinates for axis labels.
#' @param guide A sequence axis guide or `"none"`.
#'   When omitted for a circular reference, ticks and labels face inward.
#' @return A domain scale for use with sequence-path coordinates.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' ggchord(data = seqs) + geom_seq() +
#'   scale_seq_position(
#'     breaks = scales::breaks_width(200),
#'     labels = scales::label_number(suffix = " bp")
#'   ) + coord_collinear()
scale_seq_position <- function(name = waiver(), breaks = waiver(),
    minor_breaks = waiver(), labels = waiver(), coordinates = "local",
    guide = guide_seq_axis()) {
  ggchord_check_argument_names(sys.call(), sys.function())
  guide_explicit <- !missing(guide)
  coordinates <- ggchord_match_choice(coordinates, c("local", "source"), "coordinates")
  if (!is.null(name) && !is_waiver(name) &&
      (!(is.character(name) || is.expression(name)) || length(name) != 1L ||
       (is.character(name) && is.na(name))))
    ggchord_stop("scale_seq_position(): name must be one title or NULL")
  if (identical(guide, "none")) guide <- guide_seq_axis(character())
  if (!inherits(guide, "ggchord_seq_axis")) ggchord_stop("Use guide_seq_axis() or 'none'")
  out <- continuous_scale("seq_position", palette = identity,
    name = name, breaks = if (is.list(breaks)) waiver() else breaks,
    minor_breaks = if (is.list(minor_breaks)) waiver() else minor_breaks,
    labels = if (is.list(labels)) waiver() else labels, guide = "none")
  out$ggchord_axis_spec <- list(breaks = breaks, minor_breaks = minor_breaks, labels = labels)
  out$ggchord_position_coordinates <- coordinates
  out$ggchord_axis_guide <- guide
  out$ggchord_axis_guide_explicit <- guide_explicit
  out$ggchord_axis_sequences <- guide$sequences
  out$ggchord_check_overlap <- guide$check.overlap
  out
}

#' Map measurements to height within a physical track
#'
#' @param track Track ID.
#' @param name Axis title.
#' @param limits,breaks,labels,oob Continuous scale settings.
#' @param transform A continuous transformation name or object. An explicit
#'   transformation object must be monotone and invertible over the values
#'   used by the track scale.
#' @param unit Optional unit of measurement for the shared track.
#' @param domain Train a shared scale or train separately by panel.
#' @param guide A track axis guide or `"none"`.
#' @return A numerical height scale for the chosen track.
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
#' ggchord(data = seqs) +
#'   track_band("signal", thickness = 6) +
#'   geom_track_line(data = signal, aes(x = pos, y = value), track = "signal") +
#'   scale_track_value("signal", limits = c(0, 6), breaks = c(0, 3, 6),
#'                     name = "Signal", unit = "a.u.") +
#'   coord_collinear()
scale_track_value <- function(track = "default", name = waiver(), limits = NULL,
    breaks = waiver(), labels = waiver(), transform = "identity", oob = scales::censor,
    unit = NULL, domain = "shared", guide = guide_track_axis()) {
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_track_id(track)
  domain <- ggchord_match_choice(domain, c("shared", "panel"), "domain")
  if (!is.null(unit) && (!is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(unit)))
    ggchord_stop("unit must be one nonempty measurement unit")
  if (!is.null(limits) && !is.function(limits) && (!is.numeric(limits) ||
      length(limits) != 2L || all(is.na(limits)) || any(is.infinite(limits) | is.nan(limits)) ||
      isTRUE(limits[1L] > limits[2L])))
    ggchord_stop("limits must be an increasing range with at most one automatic endpoint")
  if (!is.function(oob)) ggchord_stop("oob must be a function")
  axis <- ggchord_track_axis_guide(guide)
  custom <- inherits(transform, "transform")
  if (custom) transform <- ggchord_check_track_transform(transform, track)
  sc <- continuous_scale(ggchord_track_aes("value", track), palette = identity,
    name = name, limits = limits, breaks = breaks, labels = labels,
    transform = transform, oob = oob, guide = "none")
  sc$ggchord_facet_domain <- domain
  sc$ggchord_domain_policy <- list(transform = sc$get_transformation()$name, oob = oob, limits = limits)
  sc$ggchord_track_axis <- axis$show
  sc$ggchord_track_axis_guide <- axis$guide
  sc$ggchord_track_unit <- unit
  sc$ggchord_track_id <- track
  sc$ggchord_custom_transform <- custom
  sc
}

ggchord_check_track_transform <- function(transformation, track) {
  original <- transformation
  transformation$transform <- function(x) {
    result <- original$transform(x)
    valid <- is.finite(x) & x >= original$domain[1L] & x <= original$domain[2L]
    valid <- which(valid)
    if (length(valid)) {
      index <- valid[order(x[valid])]
      index <- index[!duplicated(x[index])]
      values <- result[index]
      restored <- suppressWarnings(original$inverse(values))
      direction <- sign(diff(values))
      if (length(result) != length(x) || any(!is.finite(values)) ||
          length(restored) != length(index) || any(!is.finite(restored)) ||
          any(abs(restored - x[index]) > sqrt(.Machine$double.eps) * pmax(1, abs(x[index]))) ||
          (length(direction) && (!all(direction > 0) && !all(direction < 0))))
        ggchord_stop("Track '", track, "' requires a monotone, invertible value transform")
    }
    result
  }
  transformation
}

#' Report shared numerical limits across independent plots
#'
#' @param plots A ggchord plot or a list of independent ggchord plots.
#' @param tracks IDs of numerical tracks present in every plot.
#' @return A data frame of current limits, units and available measure evidence.
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
#' p1 <- ggchord(data = seqs) +
#'   geom_track_line(data = signal, aes(x = pos, y = value)) + coord_collinear()
#' signal2 <- transform(signal, value = value * 2)
#' p2 <- ggchord(data = seqs) +
#'   geom_track_line(data = signal2, aes(x = pos, y = value)) + coord_collinear()
#' rng <- track_limits(list(first = p1, second = p2))
#' p1 + scale_track_value(limits = rng$limits[[1]])
#' p2 + scale_track_value(limits = rng$limits[[1]])
track_limits <- function(plots, tracks = "default") {
  ggchord_check_argument_names(sys.call(), sys.function())
  if (inherits(plots, "ggchord")) plots <- list(plots)
  if (!is.list(plots) || !length(plots) || !all(vapply(plots, inherits, logical(1), "ggchord")))
    ggchord_stop("track_limits(): supply independent ggchord plots")
  ggchord_id_vector(tracks, "tracks")
  labels <- names(plots) %||% as.character(seq_along(plots))
  labels[!nzchar(labels)] <- as.character(which(!nzchar(labels)))
  reports <- lapply(plots, function(plot) {
    plot <- ggchord_render_inputs(plot)
    plot <- ggchord_register_layers(plot)
    plot <- ggchord_prepare_facet_sources(plot)
    ggchord_prepare_track_layers(plot)
  })
  compatible <- function(values, field) {
    known <- unique(values[!is.na(values)])
    if (length(known) > 1L) ggchord_stop("Incompatible track ", field)
    if (anyNA(values)) NA_character_ else known[1L]
  }
  do.call(rbind, lapply(tracks, function(id) {
    ranges <- lapply(reports, function(plot) {
      owners <- Filter(function(l) identical(l$ggchord_params$track_id, id) &&
        !identical(l$ggchord_params$track_kind, "tile"), plot$layers)
      if (!length(owners)) ggchord_stop("No numerical layer in track '", id, "'")
      values <- unlist(lapply(owners, function(l)
        ggchord_track_domain_values(ggchord_resolve_layer_input(l), l$ggchord_params)), use.names = FALSE)
      scale <- ggchord_train_track_scale(plot, id, values)
      if (!scale$has_domain) ggchord_stop("No numerical domain in track '", id, "'")
      sc <- scale$scale
      measures <- vapply(owners, function(l) l$ggchord_params$measure %||% NA_character_, character(1))
      list(limits = sort(sc$get_transformation()$inverse(sc$get_limits())),
        transform = sc$get_transformation()$name,
        unit = sc$ggchord_track_unit %||% NA_character_, measure = compatible(measures, "measure"))
    })
    collect <- function(key) vapply(ranges, `[[`, character(1), key)
    compatible(collect("transform"), "transform")
    unit <- compatible(collect("unit"), "unit")
    measure <- compatible(collect("measure"), "measure")
    data.frame(track = id, limits = I(list(range(unlist(lapply(ranges, `[[`, "limits"))))),
      unit = unit, measure = measure, plots = I(list(labels)),
      status = if (is.na(measure) || is.na(unit)) "unknown" else "available")
  }))
}

# Axis selection changes decorations, never source sequence membership.
ggchord_filter_sequence_axes <- function(layout, selection = NULL) {
  seqs <- as.character(layout$seqs %||% character())
  if (!is.null(selection)) {
    if (!is.character(selection) || anyNA(selection) || any(!nzchar(selection)) ||
        anyDuplicated(selection)) ggchord_stop(
      "guide_seq_axis(): sequences must name distinct displayed sequence IDs")
    unknown <- setdiff(selection, seqs)
    if (length(unknown)) ggchord_stop(
      "guide_seq_axis(): unknown or undisplayed axis sequence ID(s): ",
      paste(unknown, collapse = ", "), "; focused/segmented views require window IDs")
    for (component in c("axis_lines", "axis_ticks")) {
      data <- layout[[component]]
      if (is.data.frame(data) && nrow(data) && "accver" %in% names(data)) {
        layout[[component]] <- data[data$accver %in% selection, , drop = FALSE]
      }
    }
  }
  layout$axis_selection <- data.frame(accver = seqs,
    selected = if (is.null(selection)) rep(TRUE, length(seqs)) else seqs %in% selection,
    stringsAsFactors = FALSE)
  layout
}
