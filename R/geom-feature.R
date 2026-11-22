# geom-feature.R - generic feature layer (v0.9.0)

GeomChordFeatureBase <- ggproto(
  "GeomChordFeatureBase", ggplot2::GeomPolygon
)
feature_polygon_geom <- GeomChordFeatureBase

ggchord_check_feature_bounds <- function(data, lengths) {
  limit <- unname(lengths[as.character(data$accver)])
  if (any(data$start > limit | data$end > limit, na.rm = TRUE))
    ggchord_stop("geom_feature(): intervals extend outside their reference")
}

feature_geom_defaults <- feature_polygon_geom$default_aes
feature_geom_defaults$colour <- "#353A3E"
feature_geom_defaults$linewidth <- .25
feature_geom_defaults$feature_shape <- "arrow"
feature_geom <- ggproto(
  "GeomChordFeature", feature_polygon_geom,
  default_aes = feature_geom_defaults,
  optional_aes = "shape",
  extra_params = c(
    "na.rm", "segment_boundaries", "boundary_colour",
    "boundary_linewidth", "boundary_linetype", "adaptive_outline", "arrow_key"
  ),
  draw_key = function(data, params, size) {
    key_glyph_feature(data, params, size)
  },
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        segment_boundaries = TRUE,
                        boundary_colour = "#666666",
                        boundary_linewidth = .25,
                        boundary_linetype = "22",
                        adaptive_outline = TRUE, arrow_key = NULL) {
    boundary_colour <- boundary_colour %||% (data$colour[1L] %||% "black")
    component <- data$.component %||% rep("polygon", nrow(data))
    polygons <- data[is.na(component) | component == "polygon", , drop = FALSE]
    boundaries <- data[component == "boundary", , drop = FALSE]
    points <- data[component == "point", , drop = FALSE]
    grobs <- list()
    if (nrow(polygons)) {
      if (isTRUE(adaptive_outline) && "feature_shape" %in% names(polygons)) {
        small <- polygons$feature_shape %in% c(
          "compact_arrow", "promoter_arrow", "primer_arrow", "marker",
          "chevron", "lollipop"
        )
        polygons$linewidth[small] <- pmin(polygons$linewidth[small], .18)
      }
      if ("feature_shape" %in% names(polygons)) {
        line_only <- polygons$feature_shape %in% "capped_line"
        polygons$fill[line_only] <- NA
      }
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPolygon$draw_panel(
        polygons, panel_params, coord
      )
    }
    if (isTRUE(segment_boundaries) && nrow(boundaries)) {
      boundaries$colour <- boundary_colour
      boundaries$linewidth <- boundary_linewidth
      boundaries$linetype <- boundary_linetype
      if ("boundary_linetype" %in% names(boundaries)) {
        explicit <- !is.na(boundaries$boundary_linetype) &
          nzchar(as.character(boundaries$boundary_linetype))
        boundaries$linetype[explicit] <-
          as.character(boundaries$boundary_linetype[explicit])
      }
      if ("boundary_draw_linetype" %in% names(boundaries)) {
        draw_style <- !is.na(boundaries$boundary_draw_linetype) &
          nzchar(as.character(boundaries$boundary_draw_linetype))
        boundaries$linetype[draw_style] <-
          as.character(boundaries$boundary_draw_linetype[draw_style])
      }
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        boundaries, panel_params, coord, lineend = "butt",
        linejoin = "round", na.rm = na.rm
      )
    }
    if (nrow(points)) {
      if (isTRUE(adaptive_outline)) points$linewidth <- pmax(points$linewidth, .55)
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        points, panel_params, coord, lineend = "butt", linejoin = "round",
        na.rm = na.rm
      )
    }
    do.call(grobTree, grobs)
  }
)

#' Normalize common and mapped feature roles once at build time
#' @noRd
ggchord_feature_data <- function(data, fixed_shape = "block",
                                 expand_segments = TRUE) {
  ggchord_require_columns(data, c("accver", "start", "end"), "geom_feature()")
  if (!"strand" %in% names(data)) data$strand <- rep(".", nrow(data))
  for (coordinate in c("start", "end")) {
    if (is.factor(data[[coordinate]]))
      data[[coordinate]] <- as.character(data[[coordinate]])
    data[[coordinate]] <- suppressWarnings(as.numeric(data[[coordinate]]))
  }
  if (any(!is.finite(data$start) | !is.finite(data$end) |
      data$start < 1 | data$end < data$start)) {
    ggchord_stop("geom_feature(): start and end must be ordered, positive finite coordinates")
  }
  if (isTRUE(expand_segments)) data <- ggchord_expand_feature_segments(data)
  out <- as.data.frame(data, stringsAsFactors = FALSE)
  out$accver <- as.character(out$accver)
  out$strand <- as.character(out$strand)
  out$strand[is.na(out$strand)] <- "."
  if (any(!out$strand %in% c("+", "-", "+/-", ".")))
    ggchord_stop("geom_feature(): strand must be '+', '-', '+/-' or '.'")
  out$.feature_biological_strand <- out$strand
  # The renderer walks both neutral and bidirectional intervals in source
  # order. Placement and arrow heads use the separate biological strand.
  out$strand[out$strand %in% c(".", "+/-")] <- "+"
  out$anno <- rep(NA_character_, nrow(out))
  out$.feature_shape_raw <- if (is.null(fixed_shape)) {
    out$feature_shape
  } else rep(fixed_shape, nrow(out))
  out
}

ggchord_feature_shapes <- c(
  "block", "arrow", "compact_arrow", "promoter_arrow", "primer_arrow",
  "marker", "capped_line", "primer_arc", "chevron", "lollipop"
)

ggchord_feature_heads <- c("arrow", "compact_arrow", "promoter_arrow",
  "primer_arrow", "primer_arc")

# Convert absolute device units without opening a graphics device. Keep the
# names in a flexible parameter so layout can resolve them against instances.
ggchord_physical_mm <- function(value, caller, nonnegative = FALSE) {
  if (inherits(value, "unit")) {
    conversion <- c(mm = 1, cm = 10, inches = 25.4, points = 25.4 / 72.27,
      picas = 25.4 / 6.0225, bigpts = 25.4 / 72)
    unit_type <- unitType(value)
    if (any(!unit_type %in% names(conversion)))
      ggchord_stop(caller, ": use absolute physical units")
    value <- as.numeric(value) * unname(conversion[unit_type])
  }
  if (is.list(value)) return(lapply(value, ggchord_physical_mm,
    caller = caller, nonnegative = nonnegative))
  if (!is.numeric(value) || !length(value) || any(!is.finite(value)) ||
      (nonnegative && any(value < 0))) ggchord_stop(caller, ": invalid physical size")
  value
}

ggchord_units_per_mm <- function(axis = NULL) {
  units <- ggchord_render_context()$units_per_mm
  if (!is.null(axis)) units <- units[axis]
  units <- units[is.finite(units) & units > 0]
  if (length(units)) max(units) else .35 / 25.4
}

ggchord_identity_stat <- function(stat) {
  identical(stat, "identity") || inherits(stat, "StatIdentity")
}

# Split source Stat parameters from those used by the geometry renderer.
ggchord_split_source_stat <- function(stat, dots, geom) {
  if (ggchord_identity_stat(stat)) return(list(geom = dots))
  if (!inherits(stat, "Stat")) {
    stat <- layer(geom = "blank", stat = stat, position = "identity")$stat
  }
  stat_keys <- intersect(names(dots), stat$parameters(TRUE))
  list(stat = list(stat = stat, params = dots[stat_keys]),
    geom = dots[setdiff(names(dots), stat_keys)])
}

# A domain Stat receives source roles and source-row membership.
ggchord_run_source_stat <- function(data, specification, roles, caller, facet, layer) {
  if (is.null(specification) || !nrow(data)) return(data)
  stat <- specification$stat
  data$source_row <- seq_len(nrow(data))
  data$group <- data$group %||% rep(-1L, nrow(data))
  params <- c(specification$params, list(na.rm = isTRUE(layer$geom_params$na.rm)))
  members <- ggchord_source_stat_members(data, facet, layer)
  if (is.null(members)) {
    data$PANEL <- rep.int(1L, nrow(data))
  } else {
    data <- data[members$.ggchord_input_index, , drop = FALSE]
    data$PANEL <- as.integer(members$PANEL)
    for (key in facet$keys) data[[key]] <- facet$panels[[key]][data$PANEL]
  }
  params <- stat$setup_params(data, params)
  data <- stat$setup_data(data, params)
  needed <- Filter(function(role) !any(strsplit(role, "|", fixed = TRUE)[[1L]] %in% names(data)),
    stat$required_aes)
  if (length(needed)) ggchord_stop(caller, ": Stat requires ", paste(needed, collapse = ", "))
  param_names <- stat$parameters()
  compute <- function(part) {
    if (!nrow(part)) return(part)
    axes <- list(x = scale_x_continuous(), y = scale_y_continuous())
    for (axis in intersect(c("x", "y"), names(part))) axes[[axis]]$train(part[[axis]])
    computed <- do.call(stat$compute_panel, c(list(data = part, scales = axes),
      params[intersect(names(params), param_names)]))
    if (!is.data.frame(computed)) ggchord_stop(caller, ": Stat must return a data frame")
    ggchord_require_columns(computed, setdiff(roles, "group"), caller)
    if (is.null(computed$source_members)) {
      source <- computed$source_rows
      if (is.null(source) && !is.null(computed$source_row))
        source <- lapply(computed$source_row, function(i) if (!is.na(i)) i)
      if (!is.null(source))
        computed$source_members <- I(source)
    }
    computed
  }
  if (is.null(members)) return(compute(data))
  panels <- facet$panels
  ggchord_rbind_fill(lapply(seq_len(nrow(panels)), function(i) {
    out <- compute(data[data$PANEL == i, , drop = FALSE])
    out$PANEL <- out$.ggchord_stat_panel <- rep.int(i, nrow(out))
    for (key in facet$keys) out[[key]] <- rep(panels[[key]][i], nrow(out))
    out
  }))
}

ggchord_feature_layer <- function(mapping, data, shape, fixed_shape, thickness,
    stat, position, dots, na.rm, show.legend, inherit.aes, id,
    default_head_style = "shouldered", feature_role = "feature") {
  caller <- "geom_feature()"
  if (!is.character(shape) || length(shape) != 1L || !shape %in% ggchord_feature_shapes)
    ggchord_stop(caller, ": unknown shape")
  thickness <- ggchord_physical_mm(thickness, "thickness", nonnegative = TRUE)
  if (any(unlist(thickness, use.names = FALSE) <= 0))
    ggchord_stop(caller, ": thickness must be positive")
  position <- ggchord_as_feature_position(position, caller)
  dots <- ggchord_colour_dots(dots, caller)
  head_args <- c("head_length", "head_width", "head_style", "short_feature")
  boundary_args <- c("boundary_colour", "boundary_linewidth", "boundary_linetype")
  custom <- c(head_args, "unknown_strand", "segment_boundaries", boundary_args)
  explicit_heads <- intersect(head_args, names(dots))
  if (fixed_shape && length(explicit_heads) && !shape %in% ggchord_feature_heads)
    ggchord_stop(caller, ": ", explicit_heads[1L], " requires a headed geometry")
  if (length(intersect(boundary_args, names(dots))) &&
      !isTRUE(dots$segment_boundaries))
    ggchord_stop(caller, ": boundary styles require segment_boundaries = TRUE")
  head_length <- ggchord_physical_mm(dots$head_length %||% 2, "head_length", TRUE)
  head_width <- dots$head_width
  if (!is.null(head_width)) head_width <- ggchord_physical_mm(head_width, "head_width", TRUE)
  if (length(head_length) != 1L || head_length <= 0 ||
      (!is.null(head_width) && (length(head_width) != 1L || head_width <= 0)))
    ggchord_stop(caller, ": head dimensions must be positive scalars")
  head_style <- ggchord_match_arg(dots$head_style %||% default_head_style,
    c("shouldered", "flush", "triangle"))
  short <- ggchord_match_arg(dots$short_feature %||% "shrink", c("shrink", "block", "error"))
  unknown <- ggchord_match_arg(dots$unknown_strand %||% "undirected", c("undirected", "error"))
  geom_dots <- dots[setdiff(names(dots), custom)]
  ggchord_check_geom_params(geom_dots, feature_geom, caller, position, stat)
  geom_params <- ggchord_split_source_stat(stat, geom_dots, feature_geom)
  mapping <- ggchord_normalize_mapping(mapping)
  if (fixed_shape) mapping$shape <- NULL
  shape_mapped <- "shape" %in% names(mapping)
  if (shape_mapped) names(mapping)[names(mapping) == "shape"] <- "feature_shape"
  layer <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), group = integer()),
    mapping = aes(x = x, y = y, group = group),
    stat = "identity", geom = feature_geom, position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, check.param = FALSE,
    params = c(list(na.rm = na.rm,
      segment_boundaries = isTRUE(dots$segment_boundaries),
      boundary_colour = dots$boundary_colour,
      boundary_linewidth = dots$boundary_linewidth %||% 0.3,
      boundary_linetype = dots$boundary_linetype %||% 1,
      adaptive_outline = is.null(dots$linewidth),
      arrow_key = list(style = head_style,
        head_ratio = if (length(head_width) && is.numeric(thickness) && length(thickness) == 1L)
          head_width / thickness)), geom_params$geom)
  )
  layer$ggchord_type <- "gene_poly"
  layer$ggchord_params <- list(
    type = "gene", id = id, is_feature = TRUE,
    style_defaults = c(feature_shape = !fixed_shape),
    feature_mm = thickness, physical_layout = TRUE,
    head_length_mm = head_length, head_width_mm = head_width,
    head_args = explicit_heads, arrow_head_style = head_style, short_feature = short,
    unknown_strand = unknown, feature_position = position,
    gene_color_scheme = "manual", feature_role = feature_role,
    feature_shape_mapped = shape_mapped,
    feature_shape = shape, user_feature_fill = TRUE,
    feature_fill_mapping = mapping$fill, feature_fill_fixed = dots$fill,
    source_stat = geom_params$stat
  )
  layer <- ggchord_capture_layer_input(layer, data, mapping,
    c("accver", "start", "end", "strand", "feature_shape"))
  if (shape_mapped) layer$mapping$shape <- as.name(".feature_shape_raw")
  layer$ggchord_input_transform <- function(x)
    ggchord_feature_data(x, fixed_shape = if (shape_mapped) NULL else shape)
  layer$ggchord_obstacle_provider <- ggchord_entity_obstacles
  layer
}

#' Draw sequence features
#'
#' Features have a dark 0.25 mm outline by default. Compact symbols use a
#' finer outline unless a linewidth is set or mapped.
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param shape Default shape. An explicitly supplied fixed shape takes priority over its mapping.
#' @param thickness Physical thickness in mm or absolute grid units. Named values select references.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @details Provide one-based closed intervals with `accver`, `start` and `end`.
#' Supply strand `+`, `-`, `+/-` (both directions), or `.` (undirected),
#' with disjoint segments when applicable. Bidirectional arrow shapes retain
#' one source entity and use heads at both ends. The gene shortcut defaults to
#' an arrow and strand placement. The primer shortcut defaults to a directional
#' narrow path, and the region shortcut defaults to a translucent block.
#' In circular coordinates, shouldered arrows and their compact and promoter
#' forms retain the full body thickness. With no explicit `head_width`, their
#' head is 1.42 times the body thickness.
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
#' gene_colours <- c(CDS = "#187B77", tRNA = "#CC7D3A")
#'
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_feature(
#'     aes(shape = type, fill = type), data = genes,
#'     thickness = c(A = 2.5, B = 1.5)
#'   ) +
#'   scale_shape_manual(values = c(CDS = "arrow", tRNA = "block")) +
#'   scale_fill_manual(values = gene_colours) +
#'   coord_collinear()
geom_feature <- function(mapping = NULL, data = NULL, shape = "block",
    thickness = 2, stat = "identity", position = "identity", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_feature_layer(mapping, data, shape, !missing(shape), thickness,
    stat, position, list(...), na.rm, show.legend, inherit.aes, id)
}

# A plot mapping becomes known after layer assembly. Project the standard
# shape mapping once as a source role, and train its standard Scale below.
ggchord_inherit_feature_shape <- function(layer) {
  params <- layer$ggchord_params
  if (!identical(layer$ggchord_type, "gene_poly")) return(layer)
  mapping <- layer$ggchord_input_mapping
  if (!isTRUE(params$style_defaults["feature_shape"])) {
    layer$ggchord_input_mapping$shape <- NULL
    return(layer)
  }
  if (isTRUE(params$feature_shape_mapped) || is.null(mapping$shape)) return(layer)
  layer$ggchord_input_mapping$feature_shape <- mapping$shape
  layer$ggchord_input_mapping$shape <- NULL
  layer$ggchord_params$feature_shape_mapped <- TRUE
  layer$ggchord_input_transform <- function(data) ggchord_feature_data(data, fixed_shape = NULL)
  layer$mapping$shape <- as.name(".feature_shape_raw")
  layer
}
