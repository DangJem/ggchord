ggchord_label_integer <- function(value, name, minimum = 1, infinite = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < minimum ||
      (!is.finite(value) && !infinite) || value != floor(value))
    ggchord_stop(name, " must be a whole number >= ", minimum)
  value
}

#' Declare fixed label placement
#' @export
#' @param orientation Text orientation relative to the reference path or the canvas.
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
  orientation <- match.arg(orientation, c("feature", "tangent", "horizontal", "vertical"))
  side <- match.arg(side, c("positive", "negative"))
  overlap <- match.arg(overlap, c("hide", "allow"))
  fit <- match.arg(fit, c("none", "wrap"))
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
#' @param wrap Optional text wrapping specification.
#' @param max_lines Maximum number of lines after text wrapping.
#' @param min_size Minimum text size in policies that can shrink text.
#' @param max_overlaps Maximum number of conflicts permitted before hiding a label.
#' @param leader Policy for drawing a leader between the entity and text.
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
#'       arrangement = "columns", side = "positive", padding = 1.5
#'     )
#'   ) +
#'   coord_circular()
label_layout_repel <- function(arrangement = "auto", orientation = "auto",
    side = "both", offset = 1, padding = 1, fit = "wrap", wrap = NULL,
    max_lines = 2, min_size = 1.5, max_overlaps = Inf, leader = "auto") {
  ggchord_check_argument_names(sys.call(), sys.function())
  arrangement <- match.arg(arrangement, c("auto", "local", "columns"))
  orientation <- match.arg(orientation, c("auto", "feature", "tangent", "horizontal", "vertical"))
  side <- match.arg(side, c("both", "positive", "negative"))
  fit <- match.arg(fit, c("none", "wrap", "shrink", "wrap_shrink"))
  leader <- match.arg(leader, c("auto", "straight", "elbow", "none"))
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
  structure(list(arrangement = arrangement, orientation = orientation, side = side,
    offset = ggchord_physical_mm(offset, "label offset", TRUE), padding = padding,
    fit = fit, wrap = wrap, max_lines = ggchord_label_integer(max_lines, "max_lines"),
    min_size = min_size, max_overlaps = ggchord_label_integer(max_overlaps,
      "max_overlaps", minimum = 0, infinite = TRUE), leader = leader),
    class = c("LabelLayoutRepel", "ggchord_label_layout"))
}

GeomDomainLabel <- ggproto("GeomDomainLabel", ggplot2::GeomText,
  default_aes = utils::modifyList(ggplot2::GeomText$default_aes, aes(fill = NA)),
  optional_aes = c("xend", "yend", ".component", ".perimeter_callout",
    ".label_over_feature", ".arc_units", "enzyme_label", "coordinate_label",
    "label_order", "group"),
  extra_params = c("na.rm", "leader_style", "auto_contrast"),
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
      leader_style = list(), parse = FALSE, check_overlap = FALSE,
      auto_contrast = FALSE) {
    leader <- data[data$.component == "leader", , drop = FALSE]
    labels <- data[data$.component == "text", , drop = FALSE]
    grobs <- list()
    if (nrow(leader)) {
      leader$linewidth <- leader_style$linewidth %||% .3
      leader$linetype <- leader_style$linetype %||% 1
      leader$colour <- leader_style$colour %||% leader$colour
      grobs[[1L]] <- ggplot2::GeomSegment$draw_panel(leader, panel_params, coord,
        na.rm = na.rm)
    }
    boxed <- (labels$.perimeter_callout %||% FALSE) %in% TRUE & !is.na(labels$fill)
    if (auto_contrast) {
      over_feature <- (labels$.label_over_feature %||% FALSE) %in% TRUE &
        !is.na(labels$fill)
      labels$colour[over_feature] <- ggchord_contrast_colour(labels$fill[over_feature])
      labels$fill[boxed] <- ggchord_feature_callout_fill(labels$fill[boxed])
      labels$colour[boxed] <- ggchord_contrast_colour(labels$fill[boxed])
    }
    if (any(boxed)) {
      callouts <- labels[boxed, , drop = FALSE]
      callouts$linewidth <- .14
      callouts$linetype <- 1
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomLabel$draw_panel(
        callouts, panel_params, coord, parse = parse, na.rm = na.rm,
        label.padding = unit(.12, "lines"), label.r = unit(.08, "lines"),
        border.colour = "#76716A")
    }
    composite <- !boxed & !parse & !is.na(labels$enzyme_label) &
      !is.na(labels$coordinate_label) & nzchar(labels$enzyme_label) &
      nzchar(labels$coordinate_label) & labels$.perimeter_callout %in% TRUE
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
  label_order = character(), anchor_x = numeric(), anchor_y = numeric())

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
  geom_parameters <- parameters$geom
  if (label_layout$fit %in% c("shrink", "wrap_shrink")) geom_parameters$size <- NULL
  layer <- ggplot2::layer(
    data = ggchord_label_empty(), mapping = aes(
      x = x, y = y, xend = xend, yend = yend, label = label,
      size = I(size), angle = angle, hjust = hjust, vjust = vjust,
      .component = I(.component), .perimeter_callout = I(.perimeter_callout),
      .label_over_feature = I(.label_over_feature), .arc_units = I(.arc_units),
      enzyme_label = I(enzyme_label), coordinate_label = I(coordinate_label),
      label_order = I(label_order)),
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
  anchor <- match.arg(anchor, c("middle", "start", "end", "head", "tail"))
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
    circular_seam <- isTRUE(wraps[i]) && any(parts$start == 1) &&
      any(parts$end == lengths[i])
    order <- if (circular_seam) order(parts$segment_index %||% seq_len(nrow(parts))) else
      order(parts$start, parts$end)
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

ggchord_domain_label_geometry <- function(layer, input, layout, target = NULL,
    target_geometry = NULL, instances = NULL, radial = FALSE, vertical = FALSE,
    shared = NULL) {
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
  mm <- ggchord_units_per_mm(if (isTRUE(layout$collinear) && !radial)
    if (vertical) "x" else "y" else NULL)
  placed <- list()
  occupied <- shared$boxes
  metric <- ggchord_render_context()$units_per_mm
  if (length(metric) != 2L || any(!is.finite(metric) | metric <= 0))
    metric <- rep(mm, 2L)
  # Cartesian projection can have different layout units per device inch on
  # each axis. Test rotated text envelopes in physical coordinates.
  box_scale <- if (isTRUE(layout$circular)) c(1, 1) else metric * 25.4
  box_units <- if (isTRUE(layout$circular)) metric[2L] * 25.4 else 1
  use_arc <- share_perimeter && kind == "feature" && orientation == "tangent" &&
    arrangement == "local" && !isTRUE(layer$geom_params$parse) &&
    NROW(target_geometry) > 0L
  polygons <- list()
  if (use_arc || (repel && !isTRUE(layout$circular) && kind == "feature")) {
    feature_rows <- lapply(layout$layer_geometry, function(layer) layer$gene_poly)
    for (owner in names(feature_rows)) {
      d <- feature_rows[[owner]]
      if (!NROW(d)) next
      d <- d[d$.component %in% "polygon", , drop = FALSE]
      if (nrow(d)) {
        d$x <- d$x / box_scale[1L]
        d$y <- d$y / box_scale[2L]
        if (!identical(owner, p$anchor_layer)) d$source_row <- rep(NA_integer_, nrow(d))
        polygons <- c(polygons, unname(split(d, d$group)))
      }
    }
  }
  polygon_bounds <- ggchord_polygon_bounds(polygons)
  inner_leaders <- data.frame()
  bounds <- c(layout$extremes$x_min, layout$extremes$x_max)
  if (length(bounds) != 2L || any(!is.finite(bounds))) bounds <- c(-1, 1)
  for (i in order(ids, pos, label, na.last = TRUE)) {
    if (is.na(label[i]) || !nzchar(label[i])) next
    id <- ids[i]
    frame <- ggchord_panel_backbone_frame(layout, id, pos[i])
    normal <- drop(frame$normal); tangent <- c(normal[2L], -normal[1L])
    anchor_xy <- drop(frame$xy)
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
    angle <- switch(orientation, horizontal = 0,
      vertical = 90,
      atan2(tangent[2L], tangent[1L]) * 180 / pi)
    if (angle > 90) angle <- angle - 180
    if (angle < -90) angle <- angle + 180
    angle <- p$text_style$angle %||% angle
    size <- rep_len(p$label_size, nrow(input))[i]
    if (strategy$fit %in% c("shrink", "wrap_shrink") && strategy$min_size > size)
      ggchord_stop("Minimum label size exceeds the initial size")
    base_row <- data.frame(x = anchor_xy[1L], y = anchor_xy[2L],
      xend = NA_real_, yend = NA_real_, label = label[i], size = size, angle = angle,
      source_row = i, accver = id, group = i, hjust = .5, vjust = .5,
      anchor_x = anchor_xy[1L], anchor_y = anchor_xy[2L],
      anchor_layer = p$anchor_layer %||% NA_character_, anchor_position = pos[i],
      .perimeter_callout = FALSE, .label_over_feature = FALSE,
      .arc_units = NA_real_,
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
      row$angle <- if (strategy$orientation == "auto") 0 else angle
      row$.component <- "text"
      row$.perimeter_callout <- TRUE
      row$annotation_side <- "outer"
      row$label_status <- "visible"
      row$annotation_class <- if (kind == "site") "site" else if ("primer_id" %in% names(source))
        "primer" else "feature"
      row
    }
    if (share_perimeter && (arrangement == "columns" || strategy$side == "positive")) {
      placed[[length(placed) + 1L]] <- perimeter()
      next
    }
    if (share_perimeter) {
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
      slots <- expand.grid(side = allowed_sides, extra = reach, along = size * (-3:3))
      slots <- slots[order(slots$extra + abs(slots$along) * 1.3), , drop = FALSE]
    }
    inner_attachment <- if (use_arc) ggchord_feature_label_attachment(
      transform(base_row, source_row = matched[i]), polygons, FALSE, box_units)
    best <- NULL
    for (test_size in fit_sizes) {
      arc_metrics <- if (use_arc) new.env(parent = emptyenv()) else NULL
      for (j in seq_len(nrow(slots))) {
        side <- slots$side[j]
        delta <- distance[[id]][[1L]]
        xy <- anchor_xy + normal * (delta + abs(slots$extra[j])) * normal_mm * side +
          tangent * along[[id]][[1L]] * tangent_mm
        if (use_arc && slots$along[j] != 0) {
          radius <- sqrt(sum(xy^2))
          theta <- atan2(xy[2L], xy[1L]) + slots$along[j] * tangent_mm / radius
          xy <- radius * c(cos(theta), sin(theta))
        }
        if (arrangement == "columns") {
          left <- normal[1L] * side < 0 || (abs(normal[1L]) < 1e-8 &&
            anchor_xy[1L] < mean(bounds))
          xy <- c(if (left) bounds[1L] - delta * metric[1L] else bounds[2L] + delta * metric[1L],
            anchor_xy[2L] + slots$extra[j] * metric[2L])
          justify <- if (left) 1 else 0
        } else justify <- .5
        hjust <- p$text_style$hjust %||% justify
        vjust <- p$text_style$vjust %||% .5
        test_angle <- if (use_arc)
          ggchord_normalize_text_orientation(atan2(xy[2L], xy[1L]) * 180 / pi + 90,
            .5)$angle else angle
        text_at <- data.frame(x = xy[1L] / box_scale[1L],
          y = xy[2L] / box_scale[2L], text = text,
          angle = test_angle, size = test_size, hjust = hjust, vjust = vjust,
          family = p$text_style$family %||% "", fontface = p$text_style$fontface %||% 1,
          lineheight = p$text_style$lineheight %||% 1.2)
        padding <- (strategy$padding %||% 0) / (2 * 25.4)
        over <- use_arc && delta == 0 && slots$extra[j] == 0 && slots$along[j] == 0
        if (use_arc) {
          if (sqrt(sum(xy^2)) < .20) next
          box <- ggchord_arc_text_layout(
            transform(text_at, text_x = x, text_y = y, text_angle = angle),
            units_per_inch = box_units, box_padding = if (over) 0 else padding,
            line_metrics = arc_metrics, boxes_only = TRUE)$boxes
          if (over && !ggchord_feature_contains_label(box, polygons, matched[i])) next
          if (ggchord_boxes_hit_features(box, polygons, matched[i], allow_own = over,
            bounds = polygon_bounds)) next
          if (ggchord_inner_leader_hits(inner_leaders, box)) next
          line <- data.frame()
          if (!over) {
            attached <- inner_attachment
            attached$text_x <- xy[1L]
            attached$text_y <- xy[2L]
            line <- ggchord_inner_label_leader(attached, box, box_units)
          }
          if (ggchord_inner_leader_hits(line, occupied %||% data.frame()) ||
              ggchord_inner_leaders_cross(line, inner_leaders)) next
        } else {
          box <- ggchord_text_boxes(text_at, x_col = "x", y_col = "y",
            angle_col = "angle", units_per_inch = box_units, box_padding = padding)
          if (length(polygons) && ggchord_boxes_hit_features(box, polygons,
            matched[i], allow_own = delta == 0 && slots$extra[j] == 0,
            tol = padding, bounds = polygon_bounds)) next
        }
        count <- if (identical(strategy$overlap, "allow")) 0L else if (use_arc)
          as.integer(ggchord_boxes_overlap_any(box, occupied)) else
          sum(ggchord_oriented_box_overlaps(box, occupied), na.rm = TRUE)
        if (is.null(best) || count < best$count)
          best <- list(xy = xy, box = box, count = count, size = test_size,
            hjust = hjust, vjust = vjust, angle = test_angle,
            line = if (use_arc) line else NULL)
        if (count == 0L) break
      }
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
    row$.label_over_feature <- share_perimeter && kind == "feature" &&
      !is.null(target_geometry) && all(best$xy == anchor_xy)
    row$label_status <- if (best$count) "overlap_unresolved" else "visible"
    for (key in c("xmin", "xmax", "ymin", "ymax")) {
      scale <- if (startsWith(key, "x")) box_scale[1L] else box_scale[2L]
      edge <- if (endsWith(key, "min")) min else max
      row[[paste0("bbox_", key)]] <- edge(best$box[[key]]) * scale
    }
    placed[[length(placed) + 1L]] <- row
    leader <- strategy$leader %||% "none"
    if (repel && leader != "none" && any(best$xy != anchor_xy)) {
      if (use_arc) {
        line <- best$line
        if (nrow(line)) {
          inner_leaders <- rbind(inner_leaders, line)
          row$x <- line$x0; row$y <- line$y0
          row$xend <- line$x1; row$yend <- line$y1
          row$.component <- "leader"
          placed[[length(placed) + 1L]] <- row
        }
        next
      }
      row$x <- anchor_xy[1L]; row$y <- anchor_xy[2L]
      row$xend <- best$xy[1L]; row$yend <- best$xy[2L]; row$.component <- "leader"
      if (leader == "elbow" || (leader == "auto" && arrangement == "columns")) {
        corner <- c((anchor_xy[1L] + best$xy[1L]) / 2, best$xy[2L])
        first <- row
        first$xend <- corner[1L]; first$yend <- corner[2L]
        placed[[length(placed) + 1L]] <- first
        row$x <- corner[1L]; row$y <- corner[2L]
      }
      placed[[length(placed) + 1L]] <- row
    }

  }
  if (!is.null(shared)) shared$boxes <- occupied
  ggchord_rbind_fill(placed) %||% ggchord_label_empty()
}
