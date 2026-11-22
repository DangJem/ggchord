ggchord_label_integer <- function(value, name, minimum = 1, infinite = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < minimum ||
      (!is.finite(value) && !infinite) || value != floor(value))
    ggchord_stop(name, " must be a whole number >= ", minimum)
  value
}

#' Declare fixed label placement
#' @export
#' @param orientation Text orientation relative to the reference path or the canvas. Use
#'   `"tangent"` along the path or `"normal"` perpendicular to it.
#' @param side Side of the reference path on which to lay out labels.
#' @param offset Physical distance from the anchor in millimetres.
#' @param along Signed displacement in millimetres along the reference path.
#' @param overlap Policy for overlapping fixed labels.
#' @param fit Whether and how to fit text to available space.
#' @param wrap Optional text wrapping specification.
#' @return A placement declaration for compatible ggchord layers.
#' @details Layout policies arrange text in physical units while preserving the entity and its anchor. Supply the resulting declaration to a label geometry.
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
#' ggchord(seqs) +
#'   geom_gene(data = genes, id = "genes") +
#'   geom_feature_label(
#'     aes(label = label), anchor_layer = "genes",
#'     label_layout = label_layout_fixed(
#'       orientation = "horizontal", offset = 2, overlap = "allow"
#'     )
#'   ) +
#'   coord_collinear()
label_layout_fixed <- function(orientation = "feature", side = "positive",
    offset = 1, along = 0, overlap = "hide", fit = "none", wrap = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  orientation <- ggchord_match_arg(orientation, c("feature", "tangent", "normal", "horizontal", "vertical"))
  side <- ggchord_match_arg(side, c("positive", "negative"))
  overlap <- ggchord_match_arg(overlap, c("hide", "allow"))
  fit <- ggchord_match_arg(fit, c("none", "wrap"))
  if (!is.null(wrap)) {
    if (fit != "wrap") ggchord_stop("wrap requires fit = 'wrap'")
    wrap <- ggchord_label_integer(wrap, "wrap")
  }
  structure(list(orientation = orientation, side = side, overlap = overlap,
    offset = ggchord_physical_mm(offset, "label offset", TRUE),
    along = ggchord_physical_mm(along, "label along"), fit = fit, wrap = wrap),
    class = c("LabelLayoutFixed", "ggchord_label_layout"))
}

#' Declare collision-aware label placement
#' @export
#' @param arrangement Repel arrangement appropriate to the panel projection.
#' @param orientation Text orientation relative to the reference path or the canvas.
#' @param side Side of the reference path on which to lay out labels.
#' @param offset Physical distance from the anchor in millimetres.
#' @param padding Minimum clearance around a label in millimetres.
#' @param fit Whether and how to fit text to available space.
#'   `"wrap"` (default) and `"none"` keep the label size as supplied.
#'   `"shrink"` and `"wrap_shrink"` allow labels to shrink individually.
#' @param wrap Optional text wrapping specification.
#' @param max_lines Maximum number of lines after text wrapping.
#' @param min_size Minimum text size in policies that can shrink text.
#' @param max_overlaps Maximum number of conflicts permitted before hiding a label.
#' @param leader Policy for drawing a leader between the entity and text.
#' @param outer_size Optional text size in millimetres for labels outside a
#'   circular backbone. Use an absolute grid unit or a positive number.
#'   `NULL` uses the label layer's size. Mapped size and the inward fit policy
#'   keep their usual behavior.
#' @return A placement declaration for compatible ggchord layers.
#' @details Layout policies arrange text in physical units while preserving the entity and its anchor. Supply the resulting declaration to a label geometry.
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
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' ggchord(one_seq) +
#'   geom_feature(data = subset(genes, accver == "A"), id = "features") +
#'   geom_feature_label_repel(
#'     aes(label = label), anchor_layer = "features",
#'     label_layout = label_layout_repel(
#'       arrangement = "columns", side = "positive", padding = 1.5,
#'       outer_size = 3.3
#'     )
#'   ) +
#'   coord_circular()
label_layout_repel <- function(arrangement = "auto", orientation = "auto",
    side = "both", offset = 1, padding = 1, fit = "wrap", wrap = NULL,
    max_lines = 2, min_size = 1.5, max_overlaps = Inf, leader = "auto",
    outer_size = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  arrangement <- ggchord_match_arg(arrangement, c("auto", "local", "columns"))
  orientation <- ggchord_match_arg(orientation, c("auto", "feature", "tangent", "horizontal", "vertical"))
  side <- ggchord_match_arg(side, c("both", "positive", "negative"))
  fit <- ggchord_match_arg(fit, c("none", "wrap", "shrink", "wrap_shrink"))
  leader <- ggchord_match_arg(leader, c("auto", "straight", "elbow", "none"))
  if (!is.null(wrap)) {
    if (!fit %in% c("wrap", "wrap_shrink")) ggchord_stop("wrap requires a wrapping fit")
    wrap <- ggchord_label_integer(wrap, "wrap")
  }
  if (!missing(min_size) && !fit %in% c("shrink", "wrap_shrink"))
    ggchord_stop("min_size requires a shrinking fit")
  padding <- ggchord_physical_mm(padding, "label padding", TRUE)
  min_size <- ggchord_physical_mm(min_size, "minimum label size", TRUE)
  if (length(padding) != 1L || length(min_size) != 1L || min_size <= 0)
    ggchord_stop("Label padding and minimum size must be physical scalars")
  if (!is.null(outer_size)) {
    outer_size <- ggchord_physical_mm(outer_size, "outer label size")
    if (!is.numeric(outer_size) || length(outer_size) != 1L || outer_size <= 0)
      ggchord_stop("outer_size must be one positive physical size")
  }
  structure(list(arrangement = arrangement, orientation = orientation, side = side,
    offset = ggchord_physical_mm(offset, "label offset", TRUE), padding = padding,
    fit = fit, wrap = wrap, max_lines = ggchord_label_integer(max_lines, "max_lines"),
    min_size = min_size, max_overlaps = ggchord_label_integer(max_overlaps,
      "max_overlaps", minimum = 0, infinite = TRUE), leader = leader,
    outer_size = outer_size),
    class = c("LabelLayoutRepel", "ggchord_label_layout"))
}

GeomDomainLabel <- ggproto("GeomDomainLabel", ggplot2::GeomText,
  default_aes = utils::modifyList(ggplot2::GeomText$default_aes, aes(fill = NA)),
  optional_aes = c("xend", "yend", ".component", ".perimeter_callout",
    ".label_over_feature", ".arc_units", "enzyme_label", "coordinate_label",
    "label_order", "annotation_side", ".ggchord_composite_label", ".leader_crosses_backbone",
    ".leader_crosses_text", "group"),
  extra_params = c("na.rm", "leader_style", "auto_contrast", "bold_callout"),
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
      leader_style = list(), parse = FALSE, check_overlap = FALSE,
      auto_contrast = FALSE, bold_callout = FALSE) {
    leader <- data[data$.component == "leader", , drop = FALSE]
    labels <- data[data$.component == "text", , drop = FALSE]
    grobs <- list()
    if (nrow(leader)) {
      leader$linewidth <- leader_style$linewidth %||% .3
      type <- leader_style$linetype %||% "auto"
      leader$linetype <- if (identical(type, "auto"))
        ifelse(leader$.leader_crosses_backbone %in% TRUE |
          leader$.leader_crosses_text %in% TRUE, "dashed", "solid") else type
      leader$colour <- leader_style$colour %||% leader$colour
      grobs[[1L]] <- ggplot2::GeomSegment$draw_panel(leader, panel_params, coord,
        na.rm = na.rm)
    }
    boxed <- (labels$annotation_side %in% "outer" |
      labels$.perimeter_callout %in% TRUE) & !is.na(labels$fill)
    if (bold_callout) {
      bold <- labels$annotation_side %in% "outer"
      labels$fontface[bold] <- if (is.numeric(labels$fontface)) 2 else "bold"
    }
    if (auto_contrast) {
      over_feature <- (labels$.label_over_feature %||% FALSE) %in% TRUE &
        !is.na(labels$fill)
      labels$colour[over_feature] <- ggchord_contrast_colour(labels$fill[over_feature])
      labels$fill[boxed] <- ggchord_feature_callout_fill(labels$fill[boxed])
      labels$colour[boxed] <- ggchord_contrast_colour(labels$fill[boxed])
    }
    if (any(boxed)) {
      callouts <- labels[boxed, , drop = FALSE]
      callouts$linewidth <- .34
      callouts$linetype <- 1
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomLabel$draw_panel(
        callouts, panel_params, coord, parse = parse, na.rm = na.rm,
        label.padding = unit(.12, "lines"), label.r = unit(.08, "lines"),
        border.colour = "#514A43")
    }
    composite <- !boxed & !parse & labels$.ggchord_composite_label %in% TRUE
    if (any(composite)) {
      sites <- coord$transform(labels[composite, , drop = FALSE], panel_params)
      for (i in seq_len(nrow(sites))) {
        site <- sites[i, ]
        colour <- scales::alpha(site$colour, site$alpha)
        parts <- ggchord_site_label_parts(site$enzyme_label, site$coordinate_label,
          site$size, site$family, "bold", site$lineheight, colour)
        order <- if (site$label_order %in% "position_enzyme")
          c("coordinate", "enzyme") else c("enzyme", "coordinate")
        first <- parts[[order[1L]]]
        second <- parts[[order[2L]]]
        width <- grobWidth(first) + parts$gap + grobWidth(second)
        first$x <- unit(site$x, "native") - site$hjust * width
        second$x <- first$x + grobWidth(first) + parts$gap
        first$y <- second$y <- unit(site$y, "native")
        first$hjust <- second$hjust <- 0
        first$vjust <- second$vjust <- site$vjust
        grobs[[length(grobs) + 1L]] <- grobTree(first, second)
      }
    }
    arc <- !boxed & is.finite(labels$.arc_units) & !parse
    if (any(arc)) {
      curved <- labels[arc, , drop = FALSE]
      glyphs <- ggchord_arc_text_layout(
        transform(curved, text_x = x, text_y = y, text_angle = angle, text = label),
        units_per_inch = curved$.arc_units[1L]
      )$glyphs
      if (nrow(glyphs)) grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
        transform(glyphs, x = text_x, y = text_y, angle = text_angle, label = text),
        panel_params, coord, na.rm = na.rm)
    }
    plain <- !boxed & !arc & !composite
    if (any(plain)) grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
      labels[plain, , drop = FALSE], panel_params, coord, na.rm = na.rm,
      parse = parse, check_overlap = check_overlap)
    do.call(grobTree, grobs)
  }
)

ggchord_label_empty <- function() data.frame(
  x = numeric(), y = numeric(), xend = numeric(), yend = numeric(), label = character(),
  size = numeric(), angle = numeric(), source_row = integer(), accver = character(),
  .component = character(), .perimeter_callout = logical(), .label_over_feature = logical(),
  .arc_units = numeric(), enzyme_label = character(), coordinate_label = character(),
  label_order = character(), .leader_crosses_backbone = logical(),
  .ggchord_composite_label = logical(),
  annotation_side = character(), .leader_crosses_text = logical(),
  anchor_x = numeric(), anchor_y = numeric())

ggchord_domain_bold_callout <- function(layer, plot) {
  params <- layer$ggchord_params
  mapped <- c(names(ggchord_effective_mapping(layer)),
    if (isTRUE(layer$inherit.aes)) names(plot$mapping))
  identical(params$label_kind, "feature") && is.null(params$text_style$fontface) &&
    !"fontface" %in% mapped
}

ggchord_domain_label <- function(mapping, data, stat, anchor_layer, anchor, position,
    label_layout, kind, dots, na.rm, show.legend, inherit.aes, id) {
  if (!inherits(label_layout, "ggchord_label_layout"))
    ggchord_stop("label_layout must be a label layout strategy")
  if (!is.null(anchor_layer)) {
    ggchord_check_layer_id(anchor_layer)
    if (!is.null(position) || !ggchord_identity_stat(stat))
      ggchord_stop("Bound labels require an identity Stat and an unspecified Position")
  }
  position <- ggchord_as_feature_position(position %||% "identity", "label position")
  if (ggchord_local_position(position) && !identical(position$placement, "offset"))
    ggchord_stop("An independent label may use identity or position_offset()")
  if (identical(anchor, "center") &&
      (inherits(label_layout, "LabelLayoutRepel") || ggchord_local_position(position)))
    ggchord_stop("A centre label requires a fixed strategy and identity Position")
  dots <- ggchord_colour_dots(dots, "label layer")
  leader_names <- c("segment_colour", "segment_linewidth", "segment_linetype")
  leader <- dots[intersect(names(dots), leader_names)]
  if (length(leader) && !inherits(label_layout, "LabelLayoutRepel"))
    ggchord_stop("Leader styles require a repel label layout")
  names(leader) <- sub("^segment_", "", names(leader))
  dots[leader_names] <- NULL
  ggchord_check_geom_params(dots, GeomDomainLabel, "label layer", stat = stat)
  parameters <- ggchord_split_source_stat(stat, dots, GeomDomainLabel)
  initial_size <- dots$size %||% if (kind == "sequence") 3.5 else 3
  if (label_layout$fit %in% c("shrink", "wrap_shrink") && any(initial_size < label_layout$min_size))
    ggchord_stop("Minimum label size exceeds initial size")
  if (!is.null(label_layout$outer_size)) {
    if (!is.null(dots$size)) mapping$size <- NULL
    else if (!is.null(mapping$size)) label_layout$outer_size <- NULL
  }
  geom_parameters <- parameters$geom
  if (label_layout$fit %in% c("shrink", "wrap_shrink") || !is.null(label_layout$outer_size))
    geom_parameters$size <- NULL
  layer <- ggplot2::layer(
    data = ggchord_label_empty(), mapping = aes(
      x = x, y = y, xend = xend, yend = yend, label = label,
      size = I(size), angle = angle, hjust = hjust, vjust = vjust,
      .component = I(.component), .perimeter_callout = I(.perimeter_callout),
      .label_over_feature = I(.label_over_feature), .arc_units = I(.arc_units),
      enzyme_label = I(enzyme_label), coordinate_label = I(coordinate_label),
      label_order = I(label_order), annotation_side = I(annotation_side),
      .ggchord_composite_label = I(.ggchord_composite_label),
      .leader_crosses_backbone = I(.leader_crosses_backbone),
      .leader_crosses_text = I(.leader_crosses_text)),
    geom = GeomDomainLabel, stat = "identity", position = "identity", inherit.aes = inherit.aes,
    show.legend = show.legend, check.aes = FALSE, check.param = FALSE,
    params = c(list(na.rm = na.rm, leader_style = leader), geom_parameters)
  )
  layer$ggchord_type <- "domain_label"
  layer$ggchord_params <- list(id = id, type = "domain_label", label_kind = kind, anchor = anchor,
    anchor_layer = anchor_layer, label_layout = label_layout, feature_position = position,
    source_stat = parameters$stat, physical_layout = TRUE, label_size = initial_size,
    text_style = dots[intersect(names(dots), c("angle", "hjust", "vjust",
      "family", "fontface", "lineheight"))])
  roles <- switch(kind,
    sequence = "label", feature = c("accver", "start", "end", "strand", "label"),
    site = c("accver", "pos", "label"))
  ggchord_capture_layer_input(layer, data, mapping, roles)
}

#' Label source feature intervals
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional label table or function. With a named anchor, the default data consists of the target entities.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param anchor_layer ID of the target layer. Labels use the target entities after its Stat and follow their displayed placement.
#' @param anchor Position on the entity or reference path that anchors the text.
#' @param position Position adjustment. For bound labels, leave this unset to follow the target layer.
#' @param label_layout A fixed or repelled label layout declaration.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details The source entities have `accver`, `start`, `end` and a label or label mapping. A bound label follows the target layer by identity. Text size, collision handling and leader placement are independent of the underlying feature geometry. A segmented feature is anchored to its visible runs. For a circular feature label, mapping `fill` sets the background of an exterior callout. If `colour` is omitted, text over a feature or on a callout uses a contrasting colour.
#' Filled circular exterior labels have a visible outline. Circular
#' exterior feature labels use bold text by default. Set or
#' map `fontface` to override. Local inward labels can slide along an arc
#' or use a deeper text track. Automatic leaders stay on the feature's
#' radial line when it reaches a visible glyph, including off-centre glyphs.
#' Nearby inward callouts can adjust their arc positions and rows together
#' to keep text clear of other leaders.
#' Inward leaders that cross other label text are dashed. Automatic
#' leaders also use dashes when they cross their own backbone. Set
#' `segment_linetype` to override the automatic style.
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
#' ggchord(seqs) +
#'   geom_gene(data = genes, id = "genes") +
#'   geom_feature_label(aes(label = label), anchor_layer = "genes") +
#'   coord_collinear()
geom_feature_label <- function(mapping = NULL, data = NULL, stat = "identity",
    anchor_layer = NULL, anchor = "middle", position = NULL,
    label_layout = label_layout_fixed(), ...,
    na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  anchor <- ggchord_match_arg(anchor, c("middle", "start", "end", "head", "tail"))
  ggchord_domain_label(mapping, data, stat, anchor_layer, anchor, position,
    label_layout, "feature", list(...), na.rm, show.legend, inherit.aes, id)
}

#' Label source feature intervals with collision avoidance
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional label table or function. With a named anchor, the default data consists of the target entities.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param anchor_layer ID of the target layer. Labels use the target entities after its Stat and follow their displayed placement.
#' @param anchor Position on the entity or reference path that anchors the text.
#' @param position Position adjustment. For bound labels, leave this unset to follow the target layer.
#' @param label_layout A fixed or repelled label layout declaration.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details The source entities have `accver`, `start`, `end` and a label or label mapping. A bound label follows the target layer by identity. Text size, collision handling and leader placement are independent of the underlying feature geometry. A segmented feature is anchored to its visible runs. For a circular feature label, mapping `fill` sets the background of an exterior callout. If `colour` is omitted, text over a feature or on a callout uses a contrasting colour.
#' Filled circular exterior labels have a visible outline. Circular
#' exterior feature labels use bold text by default. Set or
#' map `fontface` to override. Local inward labels can slide along an arc
#' or use a deeper text track. Automatic leaders stay on the feature's
#' radial line when it reaches a visible glyph, including off-centre glyphs.
#' Nearby inward callouts can adjust their arc positions and rows together
#' to keep text clear of other leaders.
#' Inward leaders that cross other label text are dashed. Automatic
#' leaders also use dashes when they cross their own backbone. Set
#' `segment_linetype` to override the automatic style.
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
#' ggchord(seqs) +
#'   geom_gene(data = genes, id = "genes") +
#'   geom_feature_label_repel(aes(label = label), anchor_layer = "genes") +
#'   coord_collinear()
geom_feature_label_repel <- function(mapping = NULL, data = NULL, stat = "identity",
    anchor_layer = NULL, anchor = "middle", position = NULL,
    label_layout = label_layout_repel(), ...,
    na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_feature_label(mapping, data, stat, anchor_layer, anchor, position, label_layout,
    ..., na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes, id = id)
}

# Bind to the requested layer's source/statistic before expanding either layer
# into display instances. The target geometry is joined after its own layout.
ggchord_bind_domain_labels <- function(plot) {
  ids <- vapply(plot$layers, function(layer) layer$ggchord_layer_id %||% "", character(1))
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    params <- layer$ggchord_params
    anchor <- params$anchor_layer
    if (!identical(params$type, "domain_label") || is.null(anchor)) next
    target <- match(anchor, ids)
    expected <- if (params$label_kind == "feature") "gene" else "restriction_site"
    if (is.na(target) || !identical(plot$layers[[target]]$ggchord_params$type, expected))
      ggchord_stop("Label layer '", ids[i], "' requires a compatible anchor layer: ", anchor)
    original <- plot$layers[[target]]
    owner <- reconstruct_layer(original, original$data, mapping = original$mapping)
    owner$ggchord_input_transform <- NULL
    raw <- ggchord_layer_source(owner, ggchord_layer_fallback(owner, plot$ggchord$data))
    source <- ggchord_resolve_layer_input(owner, data = raw)
    layer$ggchord_params$anchor_data <- source
    layer$ggchord_params$anchor_position <- owner$ggchord_params$feature_position
    layer$ggchord_params$anchor_same_source <- is.null(layer$ggchord_input_data) ||
      (!is.null(source) && is_reference(layer$ggchord_input_data, raw))
  }
  plot
}

ggchord_label_match <- function(input, target, kind, allow_row_match) {
  key <- if (kind == "feature") "feature_id" else "site_id"
  if (key %in% names(input) && key %in% names(target)) {
    id <- function(x) paste(nchar(x$accver), x$accver, x[[key]], sep = ":")
    needed <- id(input); available <- id(target)
    if (anyNA(input[[key]]) || anyDuplicated(available[available %in% needed]))
      ggchord_stop("Ambiguous or missing ", key, " for label anchor")
    matched <- match(needed, available)
    if (anyNA(matched)) ggchord_stop("A label has no matching ", key, " in its target layer")
    return(matched)
  }
  coordinates <- if (kind == "feature") c("accver", "start", "end") else c("accver", "pos")
  if (allow_row_match && nrow(input) == nrow(target) &&
      identical(as.data.frame(input[coordinates]), as.data.frame(target[coordinates])))
    return(seq_len(nrow(input)))
  ggchord_stop("Supply a stable ", key, " to bind these labels")
}

ggchord_feature_label_position <- function(source, lengths, anchor, strand) {
  start <- source$start
  end <- source$end
  middle <- (start + end) / 2
  wraps <- (source$cross_origin %||% source$crosses_origin) %in% TRUE
  for (i in seq_len(nrow(source))) {
    segments <- source$segments[[i]]
    if (!is.data.frame(segments) || !nrow(segments)) next
    parts <- segments[!(segments$segment_type %||% "") %in% "gap", , drop = FALSE]
    if (!nrow(parts)) next
    before <- which(parts$end == lengths[i])
    after <- which(parts$start == 1)
    circular_seam <- length(before) == 1L && length(after) == 1L && before != after &&
      (isTRUE(wraps[i]) || source$start[i] > source$end[i] ||
        ggchord_segments_cross_seam(parts, lengths[i]))
    order <- if (circular_seam) order(parts$segment_index %||% seq_len(nrow(parts))) else
      order(parts$start, parts$end)
    if (circular_seam && match(before, order) > match(after, order)) order <- rev(order)
    parts <- parts[order, , drop = FALSE]
    start[i] <- parts$start[1L]
    end[i] <- tail(parts$end, 1L)
    widths <- parts$end - parts$start + 1
    edges <- c(0, cumsum(widths))
    centre <- (tail(edges, 1L) - 1) / 2
    at <- findInterval(centre, edges)
    middle[i] <- min(parts$end[at], parts$start[at] + centre - edges[at])
    # A half-position between the final and first base lies on the seam.
    if (circular_seam && middle[i] > lengths[i]) middle[i] <- lengths[i]
  }
  switch(anchor, start = start, end = end, middle = middle,
    head = ifelse(strand == "-", start, end),
    tail = ifelse(strand == "-", end, start))
}

# Reserve some space for the next label on an occupied feature band. The
# exact glyph and feature checks still decide whether a shifted slot can fit.
ggchord_arc_neighbours <- function(order, angles, bands, sizes) {
  before <- after <- vector("list", length(order))
  for (j in seq_along(order)) {
    i <- order[j]
    if (!is.finite(angles[i]) || is.na(bands[i]) || j == length(order)) next
    later <- order[seq.int(j + 1L, length(order))]
    same_band <- later[!is.na(bands[later]) & bands[later] == bands[i]]
    difference <- angles[same_band] - angles[i]
    difference <- atan2(sin(difference), cos(difference))
    for (side in c(-1, 1)) {
      available <- which(side * difference > 1e-8 & side * difference <= .45)
      if (!length(available)) next
      nearest <- available[which.min(abs(difference[available]))]
      entry <- c(angle = abs(difference[nearest]), space = sizes[same_band[nearest]])
      if (side < 0) before[[i]] <- entry else after[[i]] <- entry
    }
  }
  list(before = before, after = after)
}

ggchord_arc_label_slots <- function(slots, radius, delta, metric, width, before, after) {
  # Centred, shallow candidates may share an arc with a nearby smaller label
  # when the first label gives that neighbour room on its own source ray.
  shifts <- slots[slots$along == 0 & slots$side < 0 & delta + slots$extra > 0, ]
  r <- radius - (delta + shifts$extra) * metric
  free_before <- if (is.null(before)) Inf else r * before[["angle"]] / metric - before[["space"]]
  free_after <- if (is.null(after)) Inf else r * after[["angle"]] / metric - after[["space"]]
  low <- width / 2 - free_before
  high <- free_after - width / 2
  fits <- low <= high
  preferred <- ifelse(fits, pmax(low, pmin(0, high)), (low + high) / 2)
  preferred <- pmax(-.44 * width, pmin(.44 * width, preferred))
  shallow <- shifts$extra <= min(shifts$extra) + 1e-8
  at <- which(shallow & is.finite(preferred) & abs(preferred) > .2)
  if (length(at)) {
    adjusted <- shifts[at, , drop = FALSE]
    adjusted$along <- preferred[at]
    slots <- rbind(slots, adjusted)
  }
  cost <- slots$extra + .30 * abs(slots$along)
  if (length(at)) {
    local <- slots$side < 0 & slots$extra == min(shifts$extra)
    # On the nearest track favour the smallest shift that frees space.
    cost[local] <- slots$extra[local] + .30 * abs(slots$along[local] - preferred[at][1L])
  }
  slots[order(cost, abs(slots$along)), , drop = FALSE]
}

ggchord_domain_label_geometry <- function(layer, input, layout, target = NULL,
    target_geometry = NULL, instances = NULL, radial = FALSE, vertical = FALSE,
    shared = NULL, bold_callout = FALSE) {
  if (!NROW(input)) return(ggchord_label_empty())
  p <- layer$ggchord_params
  kind <- p$label_kind
  strategy <- p$label_layout
  repel <- inherits(strategy, "LabelLayoutRepel")
  arrangement <- strategy$arrangement %||% "local"
  if (arrangement == "auto") arrangement <- if (!isTRUE(layout$collinear) || radial) "columns" else "local"
  if (repel && arrangement == "columns" && isTRUE(layout$collinear) && !radial)
    ggchord_stop("Column label layout requires circular, chord, or radial coordinates")
  orientation <- strategy$orientation
  if (orientation == "auto") orientation <- if (arrangement == "columns") "horizontal" else "tangent"
  share_perimeter <- isTRUE(layout$circular) && !is.null(shared) && repel &&
    strategy$side != "negative" &&
    strategy$leader == "auto" &&
    !any(c("angle", "hjust", "vjust") %in% names(p$text_style)) &&
    (arrangement != "columns" || orientation == "horizontal")
  ref <- layout$sequence_reference
  ids <- as.character(input$accver)
  if (any(!ids %in% layout$seqs)) ggchord_stop("Label data contains an unknown reference")
  anchor <- p$anchor
  source <- input
  matched <- seq_len(nrow(input))
  if (!is.null(target)) {
    matched <- ggchord_label_match(input, target, kind, isTRUE(p$anchor_same_source))
    source <- target[matched, , drop = FALSE]
  }
  lens <- unname(ref$lens[ids])
  if (kind == "sequence") {
    if (anchor == "center" && !isTRUE(layout$circular))
      ggchord_stop("The centre label requires coord_circular()")
    pos <- switch(anchor, start = rep(0, nrow(input)), end = lens, lens / 2)
  } else if (kind == "site") {
    ggchord_require_columns(source, c("accver", "pos"), "Site label")
    pos <- source$pos
  } else {
    ggchord_require_columns(source, c("accver", "start", "end"), "Feature label")
    strand <- source$.feature_biological_strand %||% source$strand %||% rep(".", nrow(source))
    if (anchor %in% c("head", "tail") && anyNA(match(strand, c("+", "-"))))
      ggchord_stop("head and tail anchors require a known biological strand")
    pos <- ggchord_feature_label_position(source, lens, anchor, strand)
  }
  if (!is.numeric(pos) || any(!is.finite(pos) | pos < 0 | pos > lens))
    ggchord_stop("Label anchor positions are outside their reference")
  if (is.null(input$label) && kind != "sequence")
    ggchord_stop("Provide a label column or map aes(label = ...)")
  label <- as.character(input$label %||% input$accver)
  if (kind == "sequence" && is.null(layer$ggchord_input_mapping$label)) {
    label <- if (is.null(instances)) ids else instances$accver[match(ids, instances$key)]
  }
  offsets <- function(value, label, placement = "offset") {
    if (!is.null(instances)) value <- ggchord_instance_param(value, instances)
    ggchord_resolve_position_offset(value, layout$seqs, placement, label)
  }
  distance <- offsets(strategy$offset, "label offset")
  along <- offsets(strategy$along %||% 0, "label along")
  source_position <- p$anchor_position %||% p$feature_position
  position <- offsets(source_position$offset %||% 0, "label position",
    source_position$placement %||% "offset")
  placed <- list()
  occupied <- shared$boxes %||% data.frame()
  endpoints <- if (repel) shared$endpoints else NULL
  metric <- shared$metric %||% ggchord_annotation_metric(layout)
  # Cartesian projection can have different layout units per device inch on
  # each axis. Test rotated text envelopes in physical coordinates.
  circular <- isTRUE(layout$circular)
  outer_size <- if (circular) strategy$outer_size
  box_scale <- if (circular) c(1, 1) else unname(metric) * 25.4
  box_units <- if (isTRUE(layout$circular)) metric[2L] * 25.4 else 1
  use_arc <- circular && repel && kind == "feature" && strategy$side != "positive" &&
    orientation == "tangent" && arrangement == "local" &&
    !any(c("angle", "hjust", "vjust") %in% names(p$text_style)) &&
    !isTRUE(layer$geom_params$parse) &&
    NROW(target_geometry) > 0L
  polygons <- if (repel) shared$polygons %||% list() else list()
  own <- vapply(polygons, function(d) identical(d$.label_polygon_owner[1L], p$anchor_layer),
    logical(1L))
  own_polygons <- polygons[own]
  for (j in which(!own)) polygons[[j]]$source_row <- NA_integer_
  polygon_bounds <- if (repel) shared$polygons_bounds else NULL
  inner_leaders <- shared$leaders %||% data.frame()
  backbones <- if (repel && !circular) shared$backbones %||% data.frame() else data.frame()
  # Record both visible edges of each occupied band. A gap between those
  # edges is a real text corridor, including when a shoulder extends farther
  # than its arrow body. Other sequence rings keep their own band ordering.
  bands <- if (use_arc && length(own_polygons)) {
    outlines <- lapply(own_polygons, function(d) {
      radius <- sqrt(d$x^2 + d$y^2)
      data.frame(accver = d$accver[1L], offset = d$normal_offset[1L],
        inner = min(radius), outer = max(radius))
    })
    outlines <- do.call(rbind, outlines)
    lapply(split(outlines, outlines$accver), function(ring) {
      edges <- lapply(split(ring, ring$offset), function(track)
        data.frame(offset = track$offset[1L], inner = min(track$inner), outer = max(track$outer)))
      edges <- do.call(rbind, edges)
      edges[order(-edges$offset), , drop = FALSE]
    })
  } else list()
  bounds <- c(layout$extremes$x_min, layout$extremes$x_max)
  if (length(bounds) != 2L || any(!is.finite(bounds))) bounds <- c(-1, 1)
  label_order <- order(ids, pos, label, na.last = TRUE)
  label_sizes <- rep_len(p$label_size, nrow(input))
  radial_leader <- use_arc && strategy$leader == "auto"
  if (radial_leader) {
    initial_boxes <- occupied
    initial_leaders <- inner_leaders
    arc_choices <- vector("list", nrow(input))
    angles <- rep(NA_real_, nrow(input))
    frames <- lapply(split(seq_len(nrow(input)), ids), function(rows) {
      frame <- ggchord_panel_backbone_frame(layout, ids[rows[1L]], pos[rows])
      angles[rows] <<- atan2(frame$xy[, 2L], frame$xy[, 1L])
      list(frame = frame, rows = rows)
    })
    frame_rows <- integer(nrow(input))
    for (frame in frames) frame_rows[frame$rows] <- seq_along(frame$rows)
    offset <- target_geometry$normal_offset[match(matched, target_geometry$source_row)]
    band <- ifelse(is.finite(offset) & !is.na(label) & nzchar(label), paste(ids, offset), NA_character_)
    neighbours <- ggchord_arc_neighbours(label_order, angles, band, label_sizes)
  }
  for (i in label_order) {
    if (is.na(label[i]) || !nzchar(label[i])) next
    id <- ids[i]
    frame <- if (radial_leader) lapply(frames[[id]]$frame, function(values)
      values[frame_rows[i], , drop = FALSE]) else ggchord_panel_backbone_frame(layout, id, pos[i])
    normal <- drop(frame$normal); tangent <- c(normal[2L], -normal[1L])
    anchor_xy <- drop(frame$xy)
    backbone_outer <- if (circular) {
      outer <- layout$backbone_bounds$outer[match(id, layout$backbone_bounds$accver)]
      if (length(outer) && is.finite(outer)) outer else sqrt(sum(anchor_xy^2))
    } else Inf
    normal_mm <- 1 / sqrt(sum((normal / metric)^2))
    tangent_mm <- 1 / sqrt(sum((tangent / metric)^2))
    feature_strand <- if (kind == "feature" && strand[i] %in% c("+", "-")) strand[i] else "."
    normal_offset <- position[[id]][[feature_strand]] * normal_mm
    if (!is.null(target_geometry) && "normal_offset" %in% names(target_geometry)) {
      index <- match(matched[i], target_geometry$source_row)
      if (!is.na(index)) normal_offset <- target_geometry$normal_offset[index]
    }
    if (identical(anchor, "center")) {
      anchor_xy <- c(0, 0); normal <- c(0, 1); tangent <- c(1, 0)
      normal_offset <- 0
      normal_mm <- metric[2L]; tangent_mm <- metric[1L]
    }
    anchor_xy <- anchor_xy + normal * normal_offset
    if (circular && kind == "feature" && anchor %in% c("head", "tail") &&
        !is.null(target_geometry$.feature_tip)) {
      candidates <- which(target_geometry$source_row == matched[i] &
        target_geometry$.feature_tip %in% anchor)
      if (length(candidates)) {
        tips <- as.matrix(target_geometry[candidates, c("x", "y")])
        near <- which.min(rowSums(sweep(tips, 2L, anchor_xy)^2))
        anchor_xy <- drop(tips[near, ])
      }
    }
    angle <- switch(orientation, horizontal = 0,
      vertical = 90, normal = atan2(tangent[2L], tangent[1L]) * 180 / pi + 90,
      atan2(tangent[2L], tangent[1L]) * 180 / pi)
    if (angle > 90) angle <- angle - 180
    if (angle < -90) angle <- angle + 180
    angle <- p$text_style$angle %||% angle
    size <- label_sizes[i]
    if (strategy$fit %in% c("shrink", "wrap_shrink") && strategy$min_size > size)
      ggchord_stop("Minimum label size exceeds the initial size")
    base_row <- data.frame(x = anchor_xy[1L], y = anchor_xy[2L],
      xend = NA_real_, yend = NA_real_, label = label[i], size = size, angle = angle,
      source_row = i, accver = id, group = i, hjust = .5, vjust = .5,
      anchor_x = anchor_xy[1L], anchor_y = anchor_xy[2L],
      anchor_layer = p$anchor_layer %||% NA_character_, anchor_position = pos[i],
      .perimeter_callout = FALSE, .label_over_feature = FALSE,
      .arc_units = NA_real_,
      .ggchord_composite_label = FALSE,
      .leader_crosses_backbone = FALSE,
      .leader_crosses_text = FALSE, annotation_side = NA_character_,
      enzyme_label = if (kind == "site") input$enzyme_label[i] %||% NA_character_ else NA_character_,
      coordinate_label = if (kind == "site") input$coordinate_label[i] %||% NA_character_ else NA_character_,
      label_order = NA_character_,
      measurement = ggchord_render_context()$status %||% "estimated")
    hidden <- function(reason) {
      row <- base_row
      row$.component <- "hidden_label"
      row$label_status <- reason
      row
    }
    text <- label[i]
    if (strategy$fit %in% c("wrap", "wrap_shrink")) {
      wrap <- strategy$wrap %||% 24
      lines <- strwrap(text, width = wrap)
      if (!length(lines) || length(lines) > (strategy$max_lines %||% Inf)) {
        placed[[length(placed) + 1L]] <- hidden("wrap_limit")
        next
      }
      text <- paste(lines, collapse = "\n")
    }
    allowed_sides <- switch(strategy$side, positive = 1, negative = -1, both = c(1, -1))
    reach <- if (repel) seq(0, 8 * max(size, 1), by = max(size / 2, 0.5)) else 0
    perimeter <- function() {
      row <- base_row
      xy <- anchor_xy + normal * distance[[id]][[1L]] * normal_mm
      row$x <- xy[1L]
      row$y <- xy[2L]
      row$label <- text
      row$size <- outer_size %||% row$size
      row$angle <- if (strategy$orientation == "auto") 0 else angle
      row$.component <- "text"
      row$.perimeter_callout <- TRUE
      row$annotation_side <- "outer"
      row$label_status <- "visible"
      row$annotation_class <- if (kind == "site") "site" else if ("primer_id" %in% names(source))
        "primer" else "feature"
      if (circular && kind == "feature" && !anchor %in% c("head", "tail")) {
        attached <- ggchord_feature_label_attachment(
          transform(row, source_row = matched[i]), own_polygons, TRUE, box_units)
        row$anchor_x <- attached$anchor_x; row$anchor_y <- attached$anchor_y
      } else if (circular && kind == "site" && NROW(layout$backbone_bounds)) {
        radius <- sqrt(sum(anchor_xy^2))
        outer <- max(layout$backbone_bounds$outer, na.rm = TRUE)
        if (is.finite(outer) && radius > 1e-8) {
          row$anchor_x <- outer * row$anchor_x / radius
          row$anchor_y <- outer * row$anchor_y / radius
        }
      }
      row
    }
    if (share_perimeter && (arrangement == "columns" || strategy$side == "positive")) {
      placed[[length(placed) + 1L]] <- perimeter()
      next
    }
    if (share_perimeter || (use_arc && strategy$side == "negative")) {
      allowed_sides <- -1
      reach <- reach[reach <= if (use_arc) 18 else 6]
      if (kind == "feature" && !use_arc) {
        span <- layout$seq_arcs[[id]]
        segments <- source$segments[[i]]
        if (is.data.frame(segments) && nrow(segments)) {
          parts <- segments[!(segments$segment_type %||% "") %in% "gap", , drop = FALSE]
          feature_length <- sum(parts$end - parts$start + 1)
        } else feature_length <- source$end[i] - source$start[i] + 1
        length_on_ring <- feature_length / lens[i] *
          sum(sqrt(diff(span$x)^2 + diff(span$y)^2)) *
          sqrt(sum(anchor_xy^2)) / mean(sqrt(span$x^2 + span$y^2))
        horizontal <- ggchord_text_boxes(data.frame(x = 0, y = 0, text = text,
          size = size, angle = 0, family = p$text_style$family %||% "",
          fontface = p$text_style$fontface %||% 1), x_col = "x", y_col = "y",
          angle_col = "angle", units_per_inch = metric[2L] * 25.4)
        minimum <- if (strategy$fit %in% c("shrink", "wrap_shrink")) strategy$min_size else size
        fits <- min(size, size * length_on_ring / horizontal$bw[1L])
        if (fits < minimum) {
          placed[[length(placed) + 1L]] <- perimeter()
          next
        }
        size <- fits
      }
    }
    slots <- if (arrangement == "columns") expand.grid(side = allowed_sides,
      extra = c(0, rep(reach[-1L], each = 2L) * rep(c(1, -1), length(reach) - 1L))) else
      expand.grid(side = allowed_sides, extra = reach)
    fit_sizes <- if (strategy$fit %in% c("shrink", "wrap_shrink"))
      unique(c(size, seq(size, strategy$min_size, length.out = 8L), strategy$min_size)) else size
    slots$along <- 0
    if (use_arc) {
      radius <- sqrt(sum(anchor_xy^2))
      # With an exterior fallback, keep inward callouts associated with
      # their source neighbourhood. An explicitly inward-only layout can
      # also use the deeper free disk.
      auto_side <- share_perimeter && strategy$side == "both"
      max_depth <- if (auto_side) min(22, .43 * radius / normal_mm) else Inf
      available_depth <- max_depth - distance[[id]][[1L]]
      band <- bands[[id]]
      source_track <- match(normal_offset, band$offset)
      preferred <- numeric()
      if (!is.na(source_track)) {
        # Open enough physical text rows for each measured corridor. The
        # inner disk also provides rows, so local congestion can move labels
        # radially rather than packing them along the same circle.
        inward <- seq.int(source_track, nrow(band))
        clearance <- (.55 * size + .45) * normal_mm
        upper <- band$inner[inward] - clearance
        lower <- c(band$outer[-1L], .20)[inward] + clearance
        pitch <- (size + .45 + strategy$padding) * normal_mm
        candidates <- Map(function(a, b) {
          if (a < b) return(numeric())
          start <- a - min((a - b) / 2, .35 * size * normal_mm)
          start - seq.int(0, min(9L, floor((start - b) / pitch))) * pitch
        }, upper, lower)
        preferred <- (radius - unlist(candidates)) / normal_mm - distance[[id]][[1L]]
        preferred <- preferred[is.finite(preferred) & preferred > 0 &
          preferred <= available_depth]
      }
      steps <- if (auto_side) c(0, -.6, .6, -1.25, 1.25) else
        c(0, -.6, .6, -1.25, 1.25, -2.1, 2.1, -3, 3)
      if (radial_leader) steps <- c(steps, -.3, .3)
      along_steps <- size * steps
      reach <- reach[reach <= available_depth]
      if (length(preferred)) reach <- c(0, reach[reach >= min(preferred)])
      extra <- sort(unique(c(0, reach, preferred)))
      slots <- expand.grid(side = allowed_sides, extra = extra, along = along_steps)
      same_radius <- slots$extra == 0 & slots$along != 0
      slots <- slots[!same_radius | !length(preferred), , drop = FALSE]
      cost <- slots$extra + abs(slots$along) * 1.25
      if (length(preferred)) {
        corridor <- slots$extra %in% preferred
        cost <- cost + ifelse(corridor | slots$extra == 0, 0, .75 * size)
      }
      slots <- slots[order(slots$along != 0, cost), , drop = FALSE]
    } else if (repel && arrangement == "local" && !circular) {
      slots <- expand.grid(side = allowed_sides, extra = reach, along = size * c(0, -.6, .6, -1.3, 1.3))
      slots <- slots[order(slots$extra + abs(slots$along) * 1.2), , drop = FALSE]
    }
    inner_attachment <- base_row
    if (circular && kind == "feature" && !anchor %in% c("head", "tail")) {
      inner_attachment <- ggchord_feature_label_attachment(
        transform(base_row, source_row = matched[i]), own_polygons, FALSE, box_units)
    }
    best <- NULL
    avoid <- !identical(strategy$overlap, "allow")
    draw_leader <- repel && !identical(strategy$leader, "none")
    obstacles <- ggchord_rbind_fill(list(occupied, endpoints)) %||% data.frame()
    padding <- (strategy$padding %||% 0) / (2 * 25.4)
    for (test_size in fit_sizes) {
      arc_metrics <- if (use_arc) new.env(parent = emptyenv()) else NULL
      bold_arc_metrics <- if (use_arc && bold_callout && length(allowed_sides) > 1L)
        new.env(parent = emptyenv()) else arc_metrics
      sample <- data.frame(text_x = 0, text_y = 0, text = text, size = test_size,
        family = p$text_style$family %||% "", fontface = p$text_style$fontface %||% 1,
        lineheight = p$text_style$lineheight %||% 1.2)
      if (!use_arc) measured <- ggchord_text_boxes(sample, units_per_inch = box_units)
      candidates <- slots
      fallback <- NULL
      if (radial_leader && any(slots$side < 0 & distance[[id]][[1L]] + slots$extra > 0)) {
        width <- ggchord_text_boxes(sample, units_per_inch = box_units)$w / tangent_mm
        candidates <- ggchord_arc_label_slots(slots, sqrt(sum(anchor_xy^2)),
          distance[[id]][[1L]], tangent_mm, width, neighbours$before[[i]], neighbours$after[[i]])
      }
      measured_outer <- NULL
      for (j in seq_len(nrow(candidates))) {
        side <- candidates$side[j]
        delta <- distance[[id]][[1L]]
        xy <- anchor_xy + normal * (delta + abs(candidates$extra[j])) * normal_mm * side +
          tangent * along[[id]][[1L]] * tangent_mm
        if (candidates$along[j] != 0) {
          if (use_arc) {
            radius <- sqrt(sum(xy^2))
            if (radius <= 1e-8) next
            theta <- atan2(xy[2L], xy[1L]) + candidates$along[j] * tangent_mm / radius
            xy <- radius * c(cos(theta), sin(theta))
          } else {
            xy <- xy + tangent * candidates$along[j] * tangent_mm
          }
        }
        # Reject candidates beyond the centre rather than allowing them to
        # reappear on the opposite side of the map.
        if (use_arc && side < 0 && sum(xy * normal) < .20) next
        exterior <- circular && sqrt(sum(xy^2)) > backbone_outer + 1e-8
        if (arrangement == "columns") {
          left <- normal[1L] * side < 0 || (abs(normal[1L]) < 1e-8 &&
            anchor_xy[1L] < mean(bounds))
          xy <- c(if (left) bounds[1L] - delta * metric[1L] else bounds[2L] + delta * metric[1L],
            anchor_xy[2L] + candidates$extra[j] * metric[2L])
          justify <- if (left) 1 else 0
          exterior <- circular && sqrt(sum(xy^2)) > backbone_outer + 1e-8
        } else justify <- .5
        candidate_size <- if (exterior) outer_size %||% test_size else test_size
        hjust <- p$text_style$hjust %||% justify
        vjust <- p$text_style$vjust %||% .5
        test_angle <- if (use_arc)
          ggchord_normalize_text_orientation(atan2(xy[2L], xy[1L]) * 180 / pi + 90,
            .5)$angle else angle
        over <- use_arc && delta == 0 && candidates$extra[j] == 0 && candidates$along[j] == 0
        if (use_arc) {
          if (sqrt(sum(xy^2)) < .20) next
          sample$text_x <- xy[1L]; sample$text_y <- xy[2L]
          sample$text_angle <- test_angle; sample$hjust <- hjust; sample$vjust <- vjust
          if (!is.null(outer_size)) sample$size <- candidate_size
          if (bold_callout) sample$fontface <- if (exterior) 2 else 1
          glyph_layout <- ggchord_arc_text_layout(
            sample,
            units_per_inch = box_units, box_padding = if (over) 0 else padding,
            line_metrics = if (exterior) bold_arc_metrics else arc_metrics,
            boxes_only = TRUE)
          box <- glyph_layout$boxes
          if (over && !ggchord_feature_contains_label(box, own_polygons, matched[i])) next
        } else {
          metrics <- measured
          if (exterior && (bold_callout || !is.null(outer_size))) {
            if (is.null(measured_outer)) measured_outer <- ggchord_text_boxes(
              transform(sample, size = candidate_size,
                fontface = if (bold_callout) 2 else fontface), units_per_inch = box_units)
            metrics <- measured_outer
          }
          box <- ggchord_project_text_boxes(xy[1L] / box_scale[1L], xy[2L] / box_scale[2L],
            metrics$w, metrics$h, test_angle * pi / 180, hjust, vjust, box_units, padding)
        }
        count <- if (!avoid) 0L else if (use_arc)
          as.integer(ggchord_boxes_overlap_any(box, obstacles)) else
          sum(ggchord_oriented_box_overlaps(box, obstacles), na.rm = TRUE)
        if (share_perimeter && count > 0L) next
        if (!is.null(best) && count >= best$count) next
        allow_own <- over || (!use_arc && delta == 0 && candidates$extra[j] == 0)
        if (length(polygons) && ggchord_boxes_hit_features(box, polygons, matched[i],
          allow_own = allow_own, tol = if (circular) .004 else padding, bounds = polygon_bounds)) next
        inward_route <- use_arc && side < 0 && !exterior
        if (repel && (ggchord_inner_leader_hits(backbones, box) ||
            (!inward_route && ggchord_inner_leader_hits(inner_leaders, box)))) next
        line <- data.frame()
        normal_hit <- TRUE
        if (draw_leader && !over && any(xy != anchor_xy)) {
          attachment <- if (circular && side < 0) inner_attachment else base_row
          attachment$anchor_x <- attachment$anchor_x / box_scale[1L]
          attachment$anchor_y <- attachment$anchor_y / box_scale[2L]
          attachment$text_x <- xy[1L] / box_scale[1L]
          attachment$text_y <- xy[2L] / box_scale[2L]
          elbow <- strategy$leader == "elbow" || (strategy$leader == "auto" && arrangement == "columns")
          end_clearance <- if (inward_route) max(0, padding - .25 / 25.4) * box_units
            else .03 * box_units
          if (radial_leader && inward_route) {
            line <- ggchord_radial_label_leader(attachment, box, glyph_layout$visible, .25 * normal_mm)
            normal_hit <- !is.null(line)
          }
          if (!radial_leader || !inward_route || !normal_hit)
            line <- ggchord_domain_leader_route(attachment, box, box_units, elbow,
              shrink = end_clearance)
          crossing <- (!inward_route && ggchord_inner_leader_hits(line, obstacles)) ||
            ggchord_inner_leaders_cross(line, inner_leaders)
          if (crossing && !circular && !elbow && strategy$leader == "auto") {
            line <- ggchord_domain_leader_route(attachment, box, box_units, TRUE)
            crossing <- ggchord_inner_leader_hits(line, obstacles) ||
              ggchord_inner_leaders_cross(line, inner_leaders)
          }
          if (crossing) next
        }
        if (!normal_hit || is.null(best) || count < best$count) {
          candidate <- list(xy = xy, box = box, count = count, size = candidate_size,
            hjust = hjust, vjust = vjust, angle = test_angle,
            line = line, exterior = exterior, inward_route = inward_route)
          if (!normal_hit) {
            if (is.null(fallback) || count < fallback$count) fallback <- candidate
            next
          }
          best <- candidate
        }
        if (count == 0L) break
      }
      if (!is.null(fallback) && (is.null(best) || fallback$count < best$count)) best <- fallback
      if (!is.null(best) && best$count == 0L) break
    }
    if (share_perimeter && (is.null(best) || best$count > 0)) {
      placed[[length(placed) + 1L]] <- perimeter()
      next
    }
    if (is.null(best) || best$count > (strategy$max_overlaps %||% 0L)) {
      placed[[length(placed) + 1L]] <- hidden("overlap")
      next
    }
    occupied <- ggchord_rbind_fill(list(occupied, best$box))
    row <- base_row
    row$x <- best$xy[1L]; row$y <- best$xy[2L]
    row$label <- text; row$size <- best$size; row$hjust <- best$hjust; row$vjust <- best$vjust
    row$angle <- best$angle
    if (use_arc) row$.arc_units <- box_units
    row$.component <- "text"
    row$.label_over_feature <- (share_perimeter || use_arc) && kind == "feature" &&
      !is.null(target_geometry) && all(best$xy == anchor_xy)
    if (best$exterior) row$annotation_side <- "outer"
    if (best$inward_route) row$annotation_side <- "inner"
    row$label_status <- if (best$count) "overlap_unresolved" else "visible"
    for (key in c("xmin", "xmax", "ymin", "ymax")) {
      scale <- if (startsWith(key, "x")) box_scale[1L] else box_scale[2L]
      edge <- if (endsWith(key, "min")) min else max
      row[[paste0("bbox_", key)]] <- edge(best$box[[key]]) * scale
    }
    placed[[length(placed) + 1L]] <- row
    if (radial_leader)
      arc_choices[[i]] <- list(box = best$box, line = best$line,
        text_row = length(placed), feature = matched[i],
        radius_limits = sqrt(sum(anchor_xy^2)) - c(max_depth, distance[[id]][[1L]]) * normal_mm,
        attachment = c(inner_attachment$anchor_x, inner_attachment$anchor_y))
    line <- best$line
    if (nrow(line)) {
      inner_leaders <- rbind(inner_leaders, line)
      leader <- row[rep(1L, nrow(line)), , drop = FALSE]
      leader$x <- line$x0 * box_scale[1L]; leader$y <- line$y0 * box_scale[2L]
      leader$xend <- line$x1 * box_scale[1L]; leader$yend <- line$y1 * box_scale[2L]
      leader$.component <- "leader"
      placed[[length(placed) + 1L]] <- leader
    }

  }
  if (radial_leader) {
    revised <- ggchord_refine_circular_callouts(placed, arc_choices, initial_boxes,
      initial_leaders, endpoints, polygons, polygon_bounds, box_units, p$text_style,
      strategy$padding)
    if (!is.null(revised)) {
      placed <- revised$placed
      occupied <- revised$boxes
      inner_leaders <- revised$leaders
    }
  }
  if (!is.null(shared)) {
    shared$boxes <- occupied
    shared$leaders <- inner_leaders
  }
  ggchord_rbind_fill(placed) %||% ggchord_label_empty()
}

ggchord_domain_leader_route <- function(attachment, box, units, elbow, shrink = .03 * units) {
  if (!elbow) return(ggchord_inner_label_leader(attachment, box, units, shrink = shrink))
  corner <- (attachment$anchor_x + attachment$text_x) / 2
  first <- data.frame(x0 = attachment$anchor_x, y0 = attachment$anchor_y,
    x1 = corner, y1 = attachment$text_y)
  attachment$anchor_x <- corner
  attachment$anchor_y <- attachment$text_y
  rbind(first, ggchord_inner_label_leader(attachment, box, units, shrink = shrink))
}
