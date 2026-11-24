ggchord_site_shapes <- c("tick", "point", "diamond", "lollipop", "block")
ggchord_site_filled_shapes <- c("point", "diamond", "lollipop", "block")

ggchord_site_location_shapes <- function(location) {
  if (identical(location, "cut")) c("cut", "block", "diamond", "point", "lollipop") else
    ggchord_site_shapes
}

ggchord_site_key_spec <- function(data, params) {
  shape <- as.character(data$shape)[1L]
  if (is.na(shape)) shape <- params$site_key_shape
  spec <- ggchord_parse_feature_shape(shape)
  allowed <- ggchord_site_location_shapes(params$site_key_location)
  if (spec$name %in% allowed) spec else
    ggchord_parse_feature_shape(params$site_key_shape)
}

ggchord_site_key_shape <- function(data, params) {
  ggchord_site_key_spec(data, params)$name
}

GeomSequenceSite <- ggproto("GeomSequenceSite", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = aes(colour = "black", fill = "grey70", size = 2, alpha = 1,
    linewidth = .4, linetype = 1, shape = 21, stroke = .4),
  optional_aes = c(".component", ".shape_spec", ".shape_part", ".block_size", "group"),
  extra_params = c("na.rm", "site_key_shape", "site_key_location", "body_style_explicit"),
  draw_key = function(data, params, size) {
    spec <- ggchord_site_key_spec(data, params)
    location <- params$site_key_location %||% "site"
    width <- data$size %||% 2
    cut <- identical(location, "cut")
    lollipop <- identical(spec$name, "lollipop")
    marks <- ggchord_anchor_mark(spec, location, .5, c(.28, .72), width)
    span <- max(abs(unlist(lapply(marks, `[[`, "normal"))), width / 2)
    centre <- if (lollipop && !cut) .20 else .50
    extent <- if (cut) switch(spec$name, block = .20, lollipop = .32, .26) else
      if (lollipop) .48 else .30
    grobs <- lapply(marks, function(mark) {
      style <- ggchord_shape_part_styles(data, spec$value, mark$part, params$body_style_explicit)
      gp <- gpar(fill = scales::alpha(style$fill, style$alpha),
        col = scales::alpha(style$colour, style$alpha),
        lwd = style$linewidth * ggplot2::.pt, lty = style$linetype)
      y <- centre + extent * mark$normal / span
      if (cut && spec$name %in% c("point", "diamond"))
        y <- y + if (mark$side == "top") .14 else -.14
      # Start both chain strokes at the centre so dash patterns rotate too.
      if (cut && spec$name == "cut" && mark$side == "bottom") y <- rev(y)
      if (mark$component == "block") {
        side <- unit(min(if (cut) .36 else .58, mark$width / min(size)), "snpc")
        return(grid::rectGrob(mark$pos, y, width = side, height = side, gp = gp))
      }
      if (mark$component %in% c("point", "diamond")) {
        style$shape <- if (mark$component == "point") 21 else 23
        style$size <- if (lollipop) min(.9 * mark$width, .31 * min(size)) else
          min(mark$width, .51 * min(size))
        point <- ggplot2::draw_key_point(style, params, size)
        point$x <- unit(mark$pos, "npc")
        point$y <- unit(y, "npc")
        return(point)
      }
      if (lollipop && !is.null(spec$body_width))
        gp$lwd <- spec$body_width * ggplot2::.pt
      segmentsGrob(mark$pos, y[1L], mark$pos, y[2L], gp = gp)
    })
    if (length(grobs) == 1L) grobs[[1L]] else do.call(grobTree, grobs)
  },
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        body_style_explicit = character()) {
    data <- ggchord_shape_part_styles(data, data$.shape_spec, data$.shape_part, body_style_explicit)
    paths <- data[data$.component == "path", , drop = FALSE]
    for (value in unique(paths$.shape_spec)) {
      spec <- ggchord_parse_feature_shape(value)
      if (spec$name == "lollipop" && !is.null(spec$body_width))
        paths$linewidth[paths$.shape_spec == value] <- spec$body_width
    }
    points <- data[data$.component %in% c("point", "diamond"), , drop = FALSE]
    blocks <- data[data$.component == "block", , drop = FALSE]
    grobs <- list()
    if (nrow(paths)) grobs[[1L]] <- ggplot2::GeomPath$draw_panel(
      paths, panel_params, coord, lineend = "butt", na.rm = na.rm)
    if (nrow(points)) {
      points$shape <- ifelse(points$.component == "point", 21, 23)
      points$size <- points$.block_size
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPoint$draw_panel(
        points, panel_params, coord, na.rm = na.rm)
    }
    if (nrow(blocks)) {
      blocks <- coord$transform(blocks, panel_params)
      side <- unit(blocks$.block_size, "mm")
      grobs[[length(grobs) + 1L]] <- grid::rectGrob(blocks$x, blocks$y,
        width = side, height = side, default.units = "native", gp = gpar(
          fill = scales::alpha(blocks$fill, blocks$alpha),
          col = scales::alpha(blocks$colour, blocks$alpha),
          lwd = blocks$linewidth * ggplot2::.pt, lty = blocks$linetype))
    }
    do.call(grobTree, grobs)
  }
)

ggchord_site_layer <- function(mapping, data, shape, fixed_shape, size, fixed_size,
    stat, position, dots, na.rm, show.legend, inherit.aes, id) {
  preset <- ggchord_parse_feature_shape(shape)$name
  size <- ggchord_physical_mm(size, "geom_site(size)", TRUE)
  if (length(size) != 1L || size <= 0) ggchord_stop("geom_site(): size must be positive")
  position <- ggchord_as_feature_position(position, "geom_site()")
  if (ggchord_local_position(position) && !identical(position$placement, "offset"))
    ggchord_stop("geom_site(): use identity or position_offset()")
  dots <- ggchord_colour_dots(dots, "geom_site()")
  mapped_shape <- !fixed_shape && !is.null(mapping$shape)
  if ((fixed_shape || !mapped_shape) && !inherit.aes) {
    if ("fill" %in% c(names(dots), names(mapping)) && !preset %in% ggchord_site_filled_shapes)
      ggchord_stop("geom_site(): fill is available for point, diamond, block and lollipop symbols")
    if ("unknown_cut" %in% names(dots) && preset != "cut")
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
      .component = I(.component), .shape_spec = I(.shape_spec),
      .shape_part = I(.shape_part), .block_size = I(.block_size)),
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
#' @param shape Preset name or [shape_feature()] specification. A fixed shape takes priority over its mapping.
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
#' @details Provide `accver` and one-based `pos`. Use `"tick"` (default)
#' for a normal stroke, `"point"` for a circle, `"diamond"` for a diamond,
#' `"block"` for a square, or `"lollipop"` for a stem and circular head.
#' The lollipop extends from `pos` in the positive normal direction. Its
#' stem length is 1.5 times `size`. Configure its head size or independently
#' style the stem and head with [shape_feature()]. `position_offset()`
#' shifts the entire mark. Map shapes with `aes(shape = ...)` and
#' [ggplot2::scale_shape_manual()].
#' @examples
#' library(ggplot2)
#'
#' # Site shapes are anchored at recorded residue positions.
#' site_reference <- data.frame(
#'   accver = "demo", length = 2200L, sequence_unit = "bp", topology = "circular"
#' )
#' site_kinds <- c("tick", "point", "diamond", "block", "lollipop")
#' position_sites <- data.frame(
#'   accver = "demo",
#'   pos = c(300L, 700L, 1100L, 1500L, 1900L),
#'   kind = factor(site_kinds, levels = site_kinds),
#'   label = c("300", "700", "1100", "1500", "1900")
#' )
#' site_comparison <- ggchord(site_reference) +
#'   geom_seq(geometry = "line", colour = "#89979A", show.legend = FALSE) +
#'   geom_site(aes(shape = kind), data = position_sites, id = "position_marks", size = 3,
#'     colour = "#24474A") +
#'   geom_site_label(aes(label = label), anchor_layer = "position_marks",
#'     label_layout = layout_label_fixed(orientation = "horizontal", offset = 13),
#'     leader = "none", size = 3.2) +
#'   scale_shape_manual(name = "Site shape", values = stats::setNames(site_kinds, site_kinds),
#'     breaks = site_kinds, guide = guide_legend(nrow = 1)) +
#'   scale_seq_continuous(guide = "none") +
#'   coord_circular(fit = "labels") +
#'   theme_ggchord_minimal() +
#'   theme(legend.position = "bottom", legend.key.size = grid::unit(7, "mm"),
#'     plot.margin = margin(9, 12, 9, 12)) +
#'   labs(title = "Site shapes anchor at recorded positions")
#' site_comparison
geom_site <- function(mapping = NULL, data = NULL, shape = "tick", size = 2,
    stat = "identity", position = "identity", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_annotation_layer(mapping, data, "site", shape, !missing(shape), size,
    stat, position, list(...), na.rm, show.legend, inherit.aes, id,
    fixed_size = !missing(size))
}

#' Draw documented restriction cuts
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param shape Preset name or [shape_feature()] specification. A fixed shape takes priority over its mapping.
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
#' @details The default `"cut"` draws a separate normal stroke for each
#' known boundary in `cut_top` or `cut_bottom`. The top stroke extends in
#' the positive normal direction and the bottom stroke in the negative
#' direction. Cut legend keys place the top boundary mark at upper left
#' and the bottom boundary mark at lower right.
#' With neither boundary known, `unknown_cut = "site"` (default) draws a
#' tick at `pos`, and `unknown_cut = "error"` reports an error.
#' `"block"` draws a square on each known boundary's side. `"point"` and
#' `"diamond"` centre a glyph at each known boundary. `"lollipop"` extends
#' a stem and head outward from each known boundary. These four presets
#' omit cuts whose boundaries are both missing.
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
    fixed_size = !missing(size))
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
  source_row = integer(), .component = character(), normal_offset = numeric(),
  .shape_spec = character(), .shape_part = character(), .block_size = numeric(), cut_side = character()
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
  allowed <- ggchord_site_location_shapes(p$annotation_location)
  values <- unique(shape)
  specs <- lapply(values, ggchord_parse_feature_shape)
  by_row <- specs[match(shape, values)]
  shape_values <- shape
  shape <- vapply(by_row, `[[`, character(1), "name")
  if (anyNA(shape) || any(!shape %in% allowed))
    ggchord_stop("geom_site(): shape scale must produce ", paste(allowed, collapse = ", "))
  if (any(!shape %in% ggchord_site_filled_shapes) &&
      "fill" %in% c(names(layer$aes_params), names(visual)))
    ggchord_stop("geom_site(): fill applies to point, diamond, block and lollipop symbols")
  if (any(shape != "cut") && p$explicit_cut)
    ggchord_stop("geom_site(): unknown_cut applies only to cut symbols")
  at_cut <- identical(p$annotation_location, "cut")
  if ((at_cut || any(shape == "cut")) && identical(layout$sequence_unit, "aa"))
    ggchord_stop("geom_site(): cleavage boundaries require bp coordinates")
  boundary <- lapply(c("cut_top", "cut_bottom"), function(column) {
    values <- input[[column]] %||% rep(NA_real_, nrow(input))
    if (!ggchord_interbase_numeric(values) ||
        any(ggchord_invalid_interbase(values) | (!is.na(values) & values > lens)))
      ggchord_stop("geom_site(): cut positions must be interbase boundaries in 0..length")
    values
  })
  pieces <- list()
  mm <- ggchord_units_per_mm()
  for (row in seq_len(nrow(input))) {
    profile <- by_row[[row]]
    width <- size[row]
    if (profile$name %in% c("block", "tick", "cut")) width <- profile$body_width %||% width
    if (profile$name %in% c("point", "diamond")) width <- profile$head_width %||% width
    if (is.na(width)) next
    id <- ids[row]
    offset <- offsets[[id]][[1L]]
    cut <- c(boundary[[1L]][row], boundary[[2L]][row])
    if (profile$name == "cut" && p$unknown_cut == "error" && all(is.na(cut)))
      ggchord_stop("geom_site(): no documented cut for source row ", row)
    marks <- ggchord_anchor_mark(profile, p$annotation_location, input$pos[row], cut, width)
    for (mark in marks) {
      xy <- ggchord_track_point(rep(mark$pos, length(mark$normal)), id, 0,
        offset + mark$normal, layout)
      pieces[[length(pieces) + 1L]] <- data.frame(x = xy[, 1L], y = xy[, 2L],
        group = length(pieces) + 1L, accver = id, source_row = row,
        .component = mark$component, normal_offset = offset * mm,
        .shape_spec = shape_values[row], .shape_part = mark$part,
        .block_size = mark$width, cut_side = mark$side)
    }
  }
  ggchord_rbind_fill(pieces) %||% ggchord_site_empty()
}
