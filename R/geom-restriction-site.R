GeomSequenceSite <- ggproto("GeomSequenceSite", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = aes(colour = "black", fill = "grey70", size = 2, alpha = 1,
    linewidth = .4, linetype = 1, shape = 21, stroke = .4),
  optional_aes = c(".component", "group"),
  extra_params = c("na.rm", "site_key_shape"),
  draw_key = function(data, params, size) {
    shape <- as.character(data$shape)[1L]
    if (!shape %in% c("point", "cut", "tick")) shape <- params$site_key_shape
    if (identical(shape, "point")) {
      data$shape <- 21
      return(ggplot2::draw_key_point(data, params, size))
    }
    if (identical(shape, "cut")) {
      return(grid::polylineGrob(unit(c(.3, .3, .7, .7), "npc"),
        unit(c(.78, .5, .5, .22), "npc"),
        gp = gpar(col = scales::alpha(data$colour, data$alpha),
          lwd = data$linewidth * ggplot2::.pt, lty = data$linetype)))
    }
    ggplot2::draw_key_vline(data, params, size)
  },
  draw_panel = function(data, panel_params, coord, na.rm = FALSE) {
    paths <- data[data$.component == "path", , drop = FALSE]
    points <- data[data$.component == "point", , drop = FALSE]
    grobs <- list()
    if (nrow(paths)) grobs[[1L]] <- ggplot2::GeomPath$draw_panel(
      paths, panel_params, coord, lineend = "butt", na.rm = na.rm)
    if (nrow(points)) {
      points$shape <- 21
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPoint$draw_panel(
        points, panel_params, coord, na.rm = na.rm)
    }
    do.call(grobTree, grobs)
  }
)

ggchord_site_layer <- function(mapping, data, shape, fixed_shape, size, fixed_size,
    stat, position, dots, na.rm, show.legend, inherit.aes, id) {
  shape <- ggchord_match_arg(shape, c("tick", "point", "cut"))
  size <- ggchord_physical_mm(size, "geom_site(size)", TRUE)
  if (length(size) != 1L || size <= 0) ggchord_stop("geom_site(): size must be positive")
  position <- ggchord_as_feature_position(position, "geom_site()")
  if (ggchord_local_position(position) && !identical(position$placement, "offset"))
    ggchord_stop("geom_site(): use identity or position_offset()")
  dots <- ggchord_colour_dots(dots, "geom_site()")
  mapped_shape <- !fixed_shape && !is.null(mapping$shape)
  if ((fixed_shape || !mapped_shape) && !inherit.aes) {
    if ("fill" %in% c(names(dots), names(mapping)) && shape != "point")
      ggchord_stop("geom_site(): fill is only available for point symbols")
    if ("unknown_cut" %in% names(dots) && shape != "cut")
      ggchord_stop("geom_site(): unknown_cut is only available for cut symbols")
  }
  unknown_cut <- ggchord_match_arg(dots$unknown_cut %||% "site", c("site", "error"))
  explicit_cut <- "unknown_cut" %in% names(dots)
  dots$unknown_cut <- NULL
  if (fixed_shape) mapping$shape <- NULL
  if (fixed_size) {
    mapping$size <- NULL
    dots$size <- size
  }
  ggchord_check_geom_params(dots, GeomSequenceSite, "geom_site()", stat = stat)
  parameters <- ggchord_split_source_stat(stat, dots, GeomSequenceSite)
  layer <- ggplot2::layer(
    data = ggchord_site_empty(), mapping = aes(x = x, y = y, group = group,
      .component = I(.component)),
    geom = GeomSequenceSite, stat = "identity", position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes, check.aes = FALSE,
    params = c(list(na.rm = na.rm, site_key_shape = shape), parameters$geom)
  )
  layer$ggchord_type <- "restriction_site"
  layer$ggchord_params <- list(id = id, type = "restriction_site",
    site_shape = shape, fixed_shape = fixed_shape, site_size = size, fixed_size = fixed_size,
    unknown_cut = unknown_cut, explicit_cut = explicit_cut, feature_position = position,
    physical_layout = TRUE, source_stat = parameters$stat)
  ggchord_capture_layer_input(layer, data, mapping, c("accver", "pos", "cut_top", "cut_bottom"))
}

#' Draw source point annotations
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param shape Default shape. An explicitly supplied fixed shape takes priority over its mapping.
#' @param size Fixed symbol size in millimetres. To map size to data, use
#'   `aes(size = ...)` with a shared or [scale_for()] size scale.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details Site data supplies `accver` and a one-based `pos`. A double-strand cut uses `cut_top` and `cut_bottom` as interbase positions. Add site labels separately or bind them to the site layer.
#' In a discrete ggchord legend, the default Site tick uses a diamond-headed
#' stem. Explicit or mapped Site shapes keep their own symbols.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' sites <- data.frame(accver = "A", pos = 300, site_id = "s1", label = "EcoRI")
#'
#' cut_sites <- transform(sites, cut_top = 300, cut_bottom = 304)
#' ggchord(one_seq) +
#'   geom_seq() +
#'   geom_site(data = cut_sites, id = "sites") +
#'   geom_site_label(aes(label = label), anchor_layer = "sites") +
#'   coord_circular()
geom_site <- function(mapping = NULL, data = NULL, shape = "tick", size = 2,
    stat = "identity", position = "identity", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_annotation_layer(mapping, data, "site", shape, !missing(shape), size,
    stat, position, list(...), na.rm, show.legend, inherit.aes, id,
    fixed_size = !missing(size), strict_shape = FALSE)
}

#' Draw documented restriction cuts
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param shape Default shape. An explicitly supplied fixed shape takes priority over its mapping.
#' @param size Fixed symbol size in millimetres. To map size to data, use
#'   `aes(size = ...)` with a shared or [scale_for()] size scale.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details Site data supplies `accver` and a one-based `pos`. A double-strand cut uses `cut_top` and `cut_bottom` as interbase positions. Add site labels separately or bind them to the site layer.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' sites <- data.frame(accver = "A", pos = 300, site_id = "s1", label = "EcoRI")
#'
#' cut_sites <- transform(sites, cut_top = 300, cut_bottom = 304)
#' ggchord(one_seq) +
#'   geom_seq() +
#'   geom_restriction_site(data = cut_sites, id = "sites") +
#'   geom_site_label(aes(label = label), anchor_layer = "sites") +
#'   coord_circular()
geom_restriction_site <- function(mapping = NULL, data = NULL, shape = "cut", size = 2,
    stat = "identity", position = "identity", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_annotation_layer(mapping, data, "cut", shape, !missing(shape), size,
    stat, position, list(...), na.rm, show.legend, inherit.aes, id,
    fixed_size = !missing(size), strict_shape = FALSE)
}

#' Label source point annotations
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional label table or function. With a bound anchor, the default data consists of the target entities.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param anchor_layer Use [auto_anchor()] to follow the unique compatible layer, a layer ID to select a target, or `NULL` for independent labels.
#' @param position Position adjustment. For bound labels, leave this unset to follow the target layer.
#' @param label_layout A fixed or repelled label layout declaration.
#' @param leader How to route the optional label leader: `"auto"`, `"none"`,
#'   `"straight"` or `"elbow"`. Set fixed `segment_colour`, `segment_linewidth`,
#'   `segment_linetype` or `segment_alpha` through `...`.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details Site data supplies `accver` and a one-based `pos`. A double-strand cut uses `cut_top` and `cut_bottom` as interbase positions. Add site labels separately or bind them to the site layer.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' sites <- data.frame(accver = "A", pos = 300, site_id = "s1", label = "EcoRI")
#'
#' ggchord(one_seq) +
#'   geom_site(data = sites, id = "sites") +
#'   geom_site_label(aes(label = label), anchor_layer = "sites") +
#'   coord_circular()
geom_site_label <- function(mapping = NULL, data = NULL, stat = "identity",
    anchor_layer = auto_anchor(), position = NULL, label_layout = layout_label_fixed(),
    leader = "auto", ..., na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_domain_label(mapping, data, stat, anchor_layer, "site", position,
    label_layout, leader, "site", list(...), na.rm, show.legend, inherit.aes, id)
}

ggchord_site_empty <- function() data.frame(
  x = numeric(), y = numeric(), group = integer(), accver = character(),
  source_row = integer(), .component = character(), normal_offset = numeric()
)

# Site dimensions enter the geometry before ggplot2 draws the marks. Train
# each participating size domain once across its layers, before faceting.
ggchord_prepare_site_size_scales <- function(plot) {
  layers <- plot$layers
  mapping <- lapply(layers, function(layer) {
    if (!is.null(layer$aes_params$size)) return(NULL)
    if (!is.null(layer$ggchord_type))
      return(ggchord_start_mapping(layer$ggchord_input_mapping$size))
    if (!"size" %in% layer$geom$aesthetics()) return(NULL)
    ggchord_start_mapping(layer$mapping$size %||% if (isTRUE(layer$inherit.aes)) plot$mapping$size)
  })
  sites <- which(vapply(layers, function(layer) identical(layer$ggchord_type, "restriction_site") &&
    !isTRUE(layer$ggchord_params$fixed_size), logical(1)) & lengths(mapping) > 0L)
  if (!length(sites)) return(plot)
  previous_plot <- last_plot()
  on.exit(set_last_plot(previous_plot), add = TRUE)

  scopes <- Filter(function(scope) "size" %in% scope$aesthetics, ggchord_scale_scopes(plot))
  mapped <- which(lengths(mapping) > 0L)
  independent <- ggchord_independent_layers(plot, "size", mapped, scopes)
  scopes <- c(scopes, lapply(independent, function(i) list(selected = i)))
  global <- setdiff(mapped, unlist(lapply(scopes, `[[`, "selected")))
  domains <- c(scopes, list(list(selected = global, scale = plot$scales$get_scales("size"))))
  for (domain in domains) {
    owners <- intersect(sites, domain$selected)
    if (!length(owners)) next
    values <- lapply(domain$selected[lengths(mapping[domain$selected]) > 0L], function(i) {
      layer <- layers[[i]]
      data <- if (is.null(layer$ggchord_type)) layer$layer_data(plot$data) else
        ggchord_resolve_layer_input(layer, ggchord_layer_fallback(layer, plot$ggchord$data))
      if (!NROW(data)) return(NULL)
      cached <- ggchord_visual_value(layer$ggchord_params$visual_source, "size",
        mapping[[i]], data, seq_len(nrow(data)))
      if (!is.null(cached)) return(cached[[1L]])
      ggchord_eval_input_mapping("size", mapping[[i]], data, ggchord_layer_input_context(layer))
    })
    values <- Filter(length, values)
    if (!length(values)) next
    scale <- (domain$scale %||% ggplot2::scale_size())$clone()
    # Use ggplot2's active theme and default palette, which can be deferred
    # until its scale is initialized.
    carrier <- ggplot(mapping = aes(x = 0, y = 0, size = value)) + plot$theme
    for (value in values) carrier <- carrier + ggplot2::geom_point(data = data.frame(value = unique(value)))
    scale$aesthetics <- "size"
    trained <- ggplot_build(carrier + scale)$plot$scales$get_scales("size")
    for (i in owners) plot$layers[[i]]$ggchord_params$size_scale <- trained
  }
  plot
}

ggchord_site_geometry <- function(layer, input, layout, plot) {
  if (!NROW(input)) return(ggchord_site_empty())
  ggchord_require_columns(input, c("accver", "pos"), "geom_site()")
  p <- layer$ggchord_params
  visual <- layer$ggchord_input_mapping
  ids <- as.character(input$accver)
  lens <- unname(layout$sequence_reference$lens[ids])
  if (anyNA(lens) || !is.numeric(input$pos) || any(!is.finite(input$pos) |
      input$pos != trunc(input$pos) | input$pos < 1 | input$pos > lens))
    ggchord_stop("geom_site(): positions must be integer residues within a known reference")
  offsets <- ggchord_resolve_position_offset(p$feature_position$offset %||% 0,
    layout$seqs, "offset", "geom_site()")
  shape <- rep(p$site_shape, nrow(input))
  size <- rep(p$site_size, nrow(input))
  for (channel in c("shape", "size")) {
    if (isTRUE(p[[paste0("fixed_", channel)]])) next
    start <- ggchord_start_mapping(visual[[channel]])
    if (is.null(start)) next
    cached <- ggchord_visual_value(p$visual_source, channel, start, input, seq_len(nrow(input)))
    mapped <- if (is.null(cached)) ggchord_eval_input_mapping(channel, start, input, "geom_site()") else cached[[1L]]
    if (!length(mapped) %in% c(1L, nrow(input)))
      ggchord_stop("geom_site(): invalid length for ", channel, " mapping")
    if (channel == "shape") {
      scale <- p$shape_scale %||% plot$scales$get_scales("shape")
      if (is.null(scale)) ggchord_stop("geom_site(): shape mapping requires an explicit visual scale")
      scale <- scale$clone()
      transformed <- scale$transform(mapped)
      scale$train(transformed)
    } else {
      scale <- p$size_scale
      if (is.null(scale)) next
      transformed <- scale$transform(mapped)
    }
    mapped <- rep_len(scale$map(transformed), nrow(input))
    if (channel == "shape") shape <- as.character(mapped) else size <- mapped
  }
  allowed <- if (!isTRUE(p$strict_annotation_shape)) c("tick", "point", "cut") else
    if (p$annotation_location == "cut") "cut" else c("tick", "point")
  if (anyNA(shape) || any(!shape %in% allowed))
    ggchord_stop("geom_site(): shape scale must produce ", paste(allowed, collapse = ", "))
  if (any(shape != "point") && "fill" %in% c(names(layer$aes_params), names(visual)))
    ggchord_stop("geom_site(): fill applies only to point symbols")
  if (any(shape != "cut") && p$explicit_cut)
    ggchord_stop("geom_site(): unknown_cut applies only to cut symbols")
  if (any(shape == "cut") && identical(layout$sequence_unit, "aa"))
    ggchord_stop("geom_site(): cleavage boundaries require bp coordinates")
  boundary <- lapply(c("cut_top", "cut_bottom"), function(column) {
    values <- input[[column]] %||% rep(NA_real_, nrow(input))
    if (!ggchord_interbase_numeric(values) ||
        any(ggchord_invalid_interbase(values) | (!is.na(values) & values > lens)))
      ggchord_stop("geom_site(): cut positions must be interbase boundaries in 0..length")
    values
  })
  pieces <- list()
  for (row in seq_len(nrow(input))) {
    if (is.na(size[row])) next
    id <- ids[row]
    mm <- ggchord_units_per_mm()
    offset <- offsets[[id]][[1L]]
    add <- function(pos, extent, component) {
      xy <- ggchord_track_point(rep(pos, length(extent)), id, 0, offset + extent, layout)
      pieces[[length(pieces) + 1L]] <<- data.frame(x = xy[, 1L], y = xy[, 2L],
        group = length(pieces) + 1L, accver = id, source_row = row,
        .component = component, normal_offset = offset * mm)
    }
    half <- size[row] / 2
    if (shape[row] == "point") {
      add(input$pos[row], 0, "point")
      next
    }
    cut <- c(boundary[[1L]][row], boundary[[2L]][row])
    if (shape[row] == "tick" || all(is.na(cut))) {
      if (shape[row] == "cut" && p$unknown_cut == "error")
        ggchord_stop("geom_site(): no documented cut for source row ", row)
      add(input$pos[row], c(-half, half), "path")
      next
    }
    if (!is.na(cut[1L])) add(cut[1L], c(0, half), "path")
    if (!is.na(cut[2L])) add(cut[2L], c(-half, 0), "path")
  }
  ggchord_rbind_fill(pieces) %||% ggchord_site_empty()
}
