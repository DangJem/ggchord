# Composition owns device-local panels; source plots remain independent.

#' Compose independent plots
#' @param panels List of ggchord plots, Cartesian ggplots or
#'   nested compositions. Unnamed entries receive stable view IDs. Source
#'   endpoints require a registered source map.
#' @param layout Panel allocation from [layout_fusion()] or [layout_fusion_grid()].
#'   NULL selects an automatic grid. A layout added with `+` replaces it.
#' @param guides Guide policy: `"auto"` allows collection by a parent,
#'   `"keep"` preserves a local boundary, and `"collect"` collects outside guides.
#' @return A composition accepted by print and [ggplot2::ggsave()].
#' @export
#' @examples
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
#' p_over <- ggchord(seqs) + geom_seq() +
#'   geom_link_ribbon(data = links) + coord_chord()
#' p_detail <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#' ggchord_fusion(list(overview = p_over, detail = p_detail))
ggchord_fusion <- function(panels, layout = NULL, guides = "auto") {
  ggchord_local_console_errors()
  if (!is.list(panels) || !length(panels))
    ggchord_stop("ggchord_fusion(): panels must be a nonempty list")
  if (!all(vapply(panels, inherits, logical(1), "ggplot")))
    ggchord_stop("ggchord_fusion(): each panel must be an independent plot or Fusion")
  ids <- names(panels) %||% rep("", length(panels))
  if (anyNA(ids) || anyDuplicated(ids[nzchar(ids)]) ||
      any(nzchar(ids) & !nzchar(trimws(ids))) || any(grepl("[[:cntrl:]]", ids)))
    ggchord_stop("ggchord_fusion(): explicit panel IDs must be unique, nonmissing strings")
  number <- 0L
  for (i in which(!nzchar(ids))) {
    repeat {
      number <- number + 1L
      candidate <- paste0("view-", number)
      if (!candidate %in% ids) break
    }
    ids[i] <- candidate
  }
  layout <- layout %||% layout_fusion_grid()
  if (!inherits(layout, c("ggchord_fusion_layout", "ggchord_fusion_grid")))
    ggchord_stop("ggchord_fusion(): layout must come from layout_fusion() or layout_fusion_grid()")
  keys <- vapply(ids, ggchord_fusion_path_key, character(1))
  plot <- ggplot()
  plot$panels <- setNames(lapply(panels, function(child) list(plot = child)), keys)
  plot$fusion_areas <- layout
  plot$ribbons <- list()
  plot$zooms <- list()
  plot$styles <- list()
  plot$scope_scales <- list()
  plot$allow_overlap <- layout$allow_overlap
  plot$fusion_guides <- ggchord_match_arg(guides, c("auto", "keep", "collect"))
  class(plot) <- c("ggchord_fusion", class(plot))
  set_last_plot(plot)
  plot
}

#' Define rectangular plot allocations
#' @param areas Data frame with panel, xmin, xmax, ymin, ymax and optional z.
#' @param units Coordinates in relative canvas fractions, millimetres or inches.
#' @param allow_overlap Whether panel rectangles may overlap.
#' @return A Fusion layout declaration.
#' @export
#' @examples
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
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' areas <- data.frame(panel = c("left", "right"),
#'   xmin = c(0, .65), xmax = c(.6, 1), ymin = 0, ymax = 1)
#' ggchord_fusion(list(left = p_basic, right = p_basic),
#'   layout = layout_fusion(areas))
layout_fusion <- function(areas, units = "relative", allow_overlap = FALSE) {
  ggchord_local_console_errors()
  caller <- "layout_fusion()"
  fields <- c("xmin", "xmax", "ymin", "ymax")
  ggchord_require_columns(areas, c("panel", fields), caller)
  units <- ggchord_match_arg(units, c("relative", "mm", "in"))
  if (!is.logical(allow_overlap) || length(allow_overlap) != 1L || is.na(allow_overlap))
    ggchord_stop(caller, ": allow_overlap must be TRUE or FALSE")
  if (!nrow(areas) || !is.character(areas$panel) || anyNA(areas$panel) ||
      any(!nzchar(areas$panel)) || anyDuplicated(areas$panel))
    ggchord_stop(caller, ": each panel must be a unique nonempty string")
  if (any(!vapply(areas[fields], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(areas[fields]))) || any(as.matrix(areas[fields]) < 0) ||
      any(areas$xmin >= areas$xmax | areas$ymin >= areas$ymax) ||
      (units == "relative" && any(as.matrix(areas[fields]) > 1)))
    ggchord_stop(caller, ": areas must have finite ordered bounds within the canvas")
  if ("z" %in% names(areas) && (!is.numeric(areas$z) || any(!is.finite(areas$z))))
    ggchord_stop(caller, ": z must contain finite numbers")
  if (!allow_overlap) for (i in seq_len(nrow(areas))) {
    j <- seq_len(i - 1L)
    overlap <- pmin(areas$xmax[i], areas$xmax[j]) > pmax(areas$xmin[i], areas$xmin[j]) &
      pmin(areas$ymax[i], areas$ymax[j]) > pmax(areas$ymin[i], areas$ymin[j])
    if (any(overlap)) ggchord_stop(caller, ": overlapping areas require allow_overlap = TRUE")
  }
  structure(list(areas = as.data.frame(areas), units = units, allow_overlap = allow_overlap),
    class = "ggchord_fusion_layout")
}

#' Allocate plots on a rectangular grid
#' @param design Character matrix of panel IDs with NA for empty cells.
#' @param widths,heights Positive column and row weights.
#' @param nrow,ncol Optional positive integer grid dimensions.
#' @param gap Absolute space between cells, in millimetres or as an absolute grid unit.
#' @return A Fusion layout declaration.
#' @export
#' @examples
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
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' ggchord_fusion(list(a = p_basic, b = p_basic),
#'   layout = layout_fusion_grid(matrix(c("a", "b"), nrow = 1),
#'                              widths = c(2, 1)))
layout_fusion_grid <- function(design = NULL, widths = 1, heights = 1,
    nrow = NULL, ncol = NULL, gap = 0) {
  ggchord_local_console_errors()
  if (!is.null(design) && (!is.matrix(design) || !is.character(design) ||
      !length(design) || all(is.na(design)) || any(!is.na(design) & !nzchar(design))))
    ggchord_stop("layout_fusion_grid(): design must be a nonempty character matrix")
  for (weights in list(widths, heights)) {
    if (!is.numeric(weights) || !length(weights) || any(!is.finite(weights) | weights <= 0))
      ggchord_stop("layout_fusion_grid(): widths and heights must be positive finite weights")
  }
  for (size in list(nrow, ncol)) {
    if (!is.null(size) && (!is.numeric(size) || length(size) != 1L ||
        !is.finite(size) || size < 1 || size != floor(size)))
      ggchord_stop("layout_fusion_grid(): nrow and ncol must be positive integers")
  }
  gap <- ggchord_physical_mm(gap, "layout_fusion_grid(): gap", TRUE)
  if (!is.numeric(gap) || length(gap) != 1L)
    ggchord_stop("layout_fusion_grid(): gap must be a single absolute distance")
  if (!is.null(design) && (!is.null(nrow) && nrow != base::nrow(design) ||
      !is.null(ncol) && ncol != base::ncol(design)))
    ggchord_stop("layout_fusion_grid(): dimensions must agree with design")
  spec <- structure(list(design = design, widths = widths, heights = heights,
    nrow = nrow, ncol = ncol, gap = gap), class = "ggchord_fusion_grid")
  if (!is.null(design))
    ggchord_fusion_grid_areas(spec, unique(design[!is.na(design)]), c(1, 1), gap = 0)
  spec
}

ggchord_fusion_grid_areas <- function(spec, ids, canvas, gap = spec$gap) {
  design <- spec$design
  if (is.null(design)) {
    n <- length(ids)
    columns <- spec$ncol %||% if (is.null(spec$nrow)) ceiling(sqrt(n)) else ceiling(n / spec$nrow)
    rows <- spec$nrow %||% ceiling(n / columns)
    if (rows * columns < n)
      ggchord_stop("layout_fusion_grid(): the grid has fewer cells than views")
    design <- matrix(c(ids, rep(NA_character_, rows * columns - n)),
      nrow = rows, ncol = columns, byrow = TRUE)
  }
  if (!setequal(ids, design[!is.na(design)]))
    ggchord_stop("layout_fusion_grid(): design must allocate exactly the registered views",
      "; missing: ", paste(setdiff(ids, design), collapse = ", "),
      "; unknown: ", paste(setdiff(design[!is.na(design)], ids), collapse = ", "))
  edges <- function(weights, n, space) {
    if (!length(weights) %in% c(1L, n))
      ggchord_stop("layout_fusion_grid(): weights must match the grid dimensions")
    available <- space - (n - 1) * gap
    if (available <= 0)
      ggchord_stop("layout_fusion_grid(): gaps exceed the available canvas")
    weights <- rep_len(weights, n) / max(weights)
    sizes <- available * weights / sum(weights)
    start <- c(0, head(cumsum(sizes + gap), -1L))
    list(start = start / space, end = pmin(1, (start + sizes) / space))
  }
  x <- edges(spec$widths, ncol(design), canvas[1L])
  y <- edges(spec$heights, nrow(design), canvas[2L])
  rows <- lapply(ids, function(id) {
    cells <- which(design == id, arr.ind = TRUE)
    if (!nrow(cells)) ggchord_stop("layout_fusion_grid(): missing panel '", id, "'")
    left <- min(cells[, 2]); right <- max(cells[, 2])
    top <- min(cells[, 1]); bottom <- max(cells[, 1])
    if (nrow(cells) != (right - left + 1L) * (bottom - top + 1L))
      ggchord_stop("layout_fusion_grid(): panel '", id, "' must occupy a rectangle")
    data.frame(panel = id, xmin = x$start[left], xmax = x$end[right],
      ymin = 1 - y$end[bottom], ymax = 1 - y$start[top])
  })
  layout_fusion(do.call(rbind, rows))
}

ggchord_fusion_path_key <- function(path) {
  paste(gsub("/", "%2F", gsub("%", "%25", path, fixed = TRUE), fixed = TRUE), collapse = "/")
}

ggchord_fusion_path_ids <- function(key) {
  if (!nzchar(key)) return(character())
  utils::URLdecode(strsplit(key, "/", fixed = TRUE)[[1L]])
}

ggchord_fusion_check_path <- function(path, label) {
  if (!is.character(path) || !length(path) || anyNA(path) ||
      any(!nzchar(trimws(path))) || any(grepl("[[:cntrl:]]", path)))
    ggchord_stop(label, ": panel paths must contain nonempty string IDs")
  ggchord_fusion_path_key(path)
}

ggchord_fusion_link_roles <- c("qpanel_path", "spanel_path", "qinstance",
  "sinstance", "qfacet", "sfacet")

ggchord_fusion_link_spec <- function(layer) {
  p <- layer$ggchord_params
  kind <- layer$ggchord_type
  caller <- if (kind == "ribbon") "geom_link_ribbon()" else "geom_link_line()"
  data <- layer$ggchord_input_data
  mapping <- layer$ggchord_input_mapping
  if (!is.data.frame(data)) ggchord_stop(caller, " on Fusion: provide a relationship data frame")
  if (!is.null(p$source_stat) || !is.null(p$link_stat) ||
      ggchord_local_position(p$feature_position))
    ggchord_stop(caller, " on Fusion: use identity Stat and Position")
  unsupported <- c(p$fusion_controls, if (p$route$explicit_bulge) "bulge",
    if (p$route$route == "shared") "shared route")
  styles <- c(if (kind == "ribbon") "fill", "colour", "alpha", "linewidth", "linetype")
  required <- c("qaccver", "saccver", if (kind == "ribbon")
    c("qstart", "qend", "sstart", "send") else c("qpos", "spos"))
  roles <- union(c(required, "qstrand", "sstrand", "direction", "relative_strand"),
    ggchord_fusion_link_roles)
  fixed <- layer$aes_params
  for (name in names(fixed)) mapping[[name]] <- NULL
  params <- layer$geom_params
  options <- intersect(names(params), c("arrow", "lineend", "linejoin", "linemitre", "rule"))
  unsupported <- union(unsupported, c(setdiff(names(mapping), c(roles, styles)),
    setdiff(names(fixed), styles)))
  if (length(unsupported))
    ggchord_stop(caller, " on Fusion: unsupported parameters or aesthetics: ",
      paste(unsupported, collapse = ", "))
  if (any(lengths(fixed) != 1L)) ggchord_stop(caller, " on Fusion: fixed styles must be scalar")
  ggchord_require_columns(data, setdiff(required, names(mapping)), caller)
  structure(list(data = data, mapping = mapping, roles = roles, required = required,
    kind = kind, route = p$route$route, curvature = p$route$curvature,
    geom = layer$geom, params = params, fixed = c(fixed, params[options]),
    na.rm = params$na.rm, show.legend = layer$show.legend, id = p$id),
    class = "ggchord_fusion_link")
}

ggchord_fusion_input <- function(layer) {
  data <- as.data.frame(layer$data)
  for (role in intersect(names(layer$mapping), layer$roles)) {
    value <- eval_tidy(layer$mapping[[role]], data = layer$data)
    if (role %in% c("qfacet", "sfacet") && is.list(value) && length(names(value)))
      value <- list(value)
    if (length(value) == 1L) value <- rep(value, nrow(data))
    if (length(value) != nrow(data))
      ggchord_stop("Fusion layer '", layer$id, "': ", role, " must have one value per row")
    data[[role]] <- value
  }
  data
}

ggchord_fusion_incomplete <- function(input, required) {
  Reduce(`|`, lapply(input[required], function(value) {
    if (is.list(value)) vapply(value, function(v) !length(v) || anyNA(v), logical(1)) else is.na(value)
  }), init = rep(FALSE, nrow(input)))
}

ggchord_fusion_layer_error <- function(layer, stage, error, on_error) {
  if (is.null(on_error)) ggchord_rethrow(error)
  if (!is.null(layer$capture)) {
    for (field in if (stage == "mapped") c("mapped", "statistic", "geometry") else stage)
      layer$capture[[paste0(field, "_error")]] <- conditionMessage(error)
  }
  on_error(character(), error, layer$id)
  NULL
}

ggchord_fusion_subset_ribbons <- function(plot, ids) {
  plot$ribbons <- plot$ribbons[ids]
  plot$scope_scales <- Filter(function(scope) length(scope$layers), lapply(plot$scope_scales, function(scope) {
    scope$layers <- intersect(scope$layers, ids)
    scope
  }))
  plot
}

#' @export
`+.ggchord_fusion` <- function(e1, e2) {
  ggchord_local_console_errors()
  if (is.list(e2) && !is.object(e2)) {
    for (component in e2) e1 <- e1 + component
    return(e1)
  }
  if (inherits(e2, "Layer") && isTRUE(e2$ggchord_type %in% c("link", "ribbon")))
    e2 <- ggchord_fusion_link_spec(e2)
  if (inherits(e2, c("ggchord_fusion_layout", "ggchord_fusion_grid"))) {
    e1$fusion_areas <- e2
    e1$allow_overlap <- e2$allow_overlap
  } else if (inherits(e2, c("ggchord_fusion_link", "ggchord_fusion_zoom"))) {
    family <- if (inherits(e2, "ggchord_fusion_link")) "ribbons" else "zooms"
    prefix <- if (family == "ribbons") e2$kind else "magnify"
    ids <- c(names(e1$ribbons), names(e1$zooms))
    if (family == "ribbons") {
      e2$scale_id <- e2$id
      e2$id <- ggchord_declared_id(e2$id)
    }
    if (is.null(e2$id)) {
      number <- length(e1[[family]]) + 1L
      repeat {
        e2$id <- paste0("fusion-", prefix, "-", number)
        if (!e2$id %in% ids) break
        number <- number + 1L
      }
    }
    if (!is.character(e2$id) || length(e2$id) != 1L || is.na(e2$id) ||
        !nzchar(e2$id) || e2$id %in% ids)
      ggchord_stop("Fusion: each root layer must have a unique nonempty ID")
    if (family == "zooms") {
      registered <- ggchord_magnify_register(e1, e2)
      e1 <- registered$plot
      e2 <- registered$spec
    }
    e1[[family]][[e2$id]] <- unclass(e2)
    e1$active_layer <- e2$id
  } else if (inherits(e2, "ggchord_palette_preset")) {
    e1$palette_preset <- e2
  } else if (inherits(e2, "ggchord_scoped_scale")) {
    if (isTRUE(e2$y) || !is.null(e2$tracks))
      ggchord_stop("Fusion: scale_for() accepts visual scopes on Link or child layers")
    if (is.data.frame(e2$layers)) {
      e2$targets <- ggchord_output_selection(e2$layers, "Fusion scale_for()")
      if (!length(e2$targets)) ggchord_stop("Fusion scale_for(): select at least one layer")
      e1$fusion_shared_scales <- c(e1$fusion_shared_scales, list(e2))
    } else {
      e2$layers <- e2$layers %||% e1$active_layer
      if (is.null(e2$layers) || any(!e2$layers %in% names(e1$ribbons)))
        ggchord_stop("Fusion: add a Link or explicitly select a visual layer before scale_for()")
      e1$scope_scales <- c(e1$scope_scales, list(e2))
    }
  } else if (inherits(e2, "ggchord_layout_policy")) {
    e1$fusion_policy <- unclass(e2)
  } else if (inherits(e2, c("Scale", "theme", "Guides", "ggplot2::labels"))) {
    if (inherits(e2, "theme")) ggchord_check_theme_elements(e2)
    if (inherits(e2, "Scale") && any(e2$aesthetics %in% c("x", "y")))
      ggchord_stop("Fusion: root scales map visual aesthetics")
    e1$styles <- c(e1$styles, list(e2))
  } else ggchord_stop("Fusion accepts Link Geoms, Magnify, Fusion layouts, render policies and visual plot components")
  set_last_plot(e1)
  e1
}

# Compare native graphical contents; generated grob names are not semantics.
ggchord_fusion_guide_signature <- function(x) {
  if (inherits(x,"unit")) return(list(unit=as.character(x)))
  if (is.data.frame(x)) {
    x$name <- NULL
    return(x)
  }
  if (!is.list(x)) return(x)
  x$name <- NULL; x$childrenOrder <- NULL
  if (!is.null(x$children)) names(x$children) <- NULL
  if (!is.null(x$grobs)) names(x$grobs) <- NULL
  x <- lapply(x, ggchord_fusion_guide_signature)
  x
}

ggchord_fusion_guide_semantics <- function(built, guide, source) {
  aesthetics <- guide$aesthetic
  details <- lapply(built$plot$scales$scales, function(scale) {
    kept <- intersect(aesthetics, scale$aesthetics)
    if (!length(kept)) return(NULL)
    breaks <- scale$get_breaks()
    list(aesthetics = sub("^ggchord_scope_[0-9]+_", "", kept),
      type = if (scale$is_discrete()) "discrete" else class(scale)[1L],
      limits = scale$get_limits(), mapped_breaks = scale$map(breaks),
      transform = if (!scale$is_discrete()) scale$get_transformation()$name,
      breaks = breaks, labels = scale$get_labels(breaks))
  })
  variables <- lapply(built$plot$layers, function(layer) {
    mapping <- layer$mapping
    if (isTRUE(layer$inherit.aes))
      mapping <- c(mapping, built$plot$mapping[setdiff(names(built$plot$mapping), names(mapping))])
    include <- intersect(names(mapping), aesthetics)
    if (!length(include) || identical(layer$show.legend, FALSE)) return(NULL)
    data <- layer$ggchord_resolved_input %||% layer$ggchord_input_data %||% layer$data
    if (!is.data.frame(data)) data <- built$plot$data
    values <- lapply(mapping[include], function(value) {
      fields <- intersect(all.vars(get_expr(value)), names(data))
      units <- lapply(data[fields], function(column) list(units = attr(column, "units"),
        unit = attr(column, "unit"), measure = attr(column, "measure")))
      list(expression = as_label(value), columns = units)
    })
    list(mapping = values, unit = layer$ggchord_params$unit,
      measure = layer$ggchord_params$measure)
  })
  variables <- unique(Filter(Negate(is.null), variables))
  scopes <- Filter(Negate(is.null), details)
  list(title = guide$title, variables = variables, scales = scopes,
    unresolved = if (!length(variables) || !length(scopes) ||
      any(vapply(variables, function(v) any(vapply(v$mapping, function(x)
        length(x$columns) != 1L, logical(1))), logical(1)))) source)
}

ggchord_fusion_outside_sides <- c("top", "bottom", "left", "right")

ggchord_fusion_guide_positions <- function(built) {
  default <- built$plot$theme$legend.position %||% "right"
  if (length(default) != 1L) default <- "inside"
  vapply(built$plot$guides$params, function(guide) guide$position %||% default, character(1))
}

ggchord_fusion_hide_guides <- function(built, sides = built$plot$ggchord_fusion_omit_guides) {
  if (!length(sides)) return(built)
  guides <- built$plot$guides
  keep <- !ggchord_fusion_guide_positions(built) %in% sides
  built$plot$guides <- ggproto(NULL, guides, guides = guides$guides[keep],
    params = guides$params[keep])
  built
}

ggchord_fusion_collect_guides <- function(x, partial = FALSE) {
  boxes <- list()
  retained <- list()
  root_theme <- Reduce(`+`, Filter(function(s) inherits(s, "theme"), x$styles), init = theme())
  theme_names <- grep("^legend[.]", names(root_theme), value = TRUE)
  theme_names <- setdiff(theme_names, c("legend.position", "legend.position.inside",
    "legend.justification.inside"))
  legend_theme <- if (length(theme_names)) do.call(theme, as.list(root_theme)[theme_names]) else theme()
  collect <- function(built, path) {
    # The child owns guide training. Root theme settings only style the
    # collected copy, while Guide-level overrides retain native precedence.
    changed <- built
    class(changed) <- setdiff(class(changed), "ggchord_built")
    child_theme <- changed$plot$theme
    changed$plot$theme <- child_theme + legend_theme
    table <- ggplot_gtable(changed)
    params <- built$plot$guides$params
    positions <- ggchord_fusion_guide_positions(built)
    at <- which(grepl("^guide-box-", table$layout$name) &
      !vapply(table$grobs, inherits, logical(1), "zeroGrob"))
    sides <- character()
    for (i in at) {
      side <- sub("^guide-box-", "", table$layout$name[i])
      if (!side %in% ggchord_fusion_outside_sides) {
        retained[[length(retained) + 1L]] <<- list(source = path, side = side, reason = "local_inside")
        next
      }
      box <- table$grobs[[i]]
      ids <- which(box$layout$name == "guides")
      candidates <- params[positions == side]
      if (length(candidates) != length(ids) || !length(ids)) {
        retained[[length(retained) + 1L]] <<- list(source = path, side = side, reason = "unknown_guide_ownership")
        next
      }
      for (j in seq_along(ids)) {
        semantics <- ggchord_fusion_guide_semantics(built, candidates[[j]], path)
        boxes[[length(boxes) + 1L]] <<- list(side = side, grob = box$grobs[[ids[j]]],
          template = box, source = path, semantics = semantics,
          position_lock = candidates[[j]]$position)
      }
      sides <- c(sides, side)
    }
    sides
  }
  attempt <- function(expr) tryCatch(expr, error = function(e) {
    if (!partial) ggchord_rethrow(e)
    NULL
  })
  walk <- function(plot, path) {
    if (inherits(plot, "ggchord_fusion")) {
      if (identical(plot$fusion_guides, "keep")) {
        retained[[length(retained) + 1L]] <<- list(source = path, reason = "nested_keep")
        return(plot)
      }
      plot$fusion_guides <- "auto"
      if (length(plot$ribbons))
        plot$ggchord_fusion_omit_guides <- attempt(collect(ggchord_fusion_build_styles(plot)$built, path))
      for (id in names(plot$panels)) plot$panels[[id]]$plot <-
        walk(plot$panels[[id]]$plot, paste(path, id, sep = "/"))
      return(plot)
    }
    plot$ggchord_fusion_omit_guides <- attempt(collect(ggchord_fusion_build_leaf(plot), path))
    plot
  }
  for (id in names(x$panels)) x$panels[[id]]$plot <- walk(x$panels[[id]]$plot, id)
  list(plot = x, boxes = boxes, retained = retained)
}

# Native guides keep their physical text/key sizes. Only boxes are arranged.
ggchord_fusion_shared_guides <- function(styles, boxes, x, canvas) {
  if (!length(boxes)) return(styles)
  position <- NULL
  for (style in x$styles) if (inherits(style, "theme") && !is.null(style$legend.position))
    position <- style$legend.position
  if (identical(position, "none")) return(styles)
  if (length(position) == 1L && position %in% ggchord_fusion_outside_sides) {
    moved <- Filter(function(box) !is.null(box$position_lock) && box$side != position, boxes)
    styles$report$guide_position_overrides <- lapply(moved, function(box)
      list(source = box$source, guide_position = box$side, composition_position = position))
    boxes <- lapply(boxes, function(box) {
      box$side <- position
      box
    })
  }
  unique_boxes <- list(); signatures <- list()
  for (box in boxes) {
    signature <- list(semantics = box$semantics, appearance = ggchord_fusion_guide_signature(box$grob))
    duplicate <- which(vapply(seq_along(signatures), function(i)
      identical(unique_boxes[[i]]$side,box$side) && identical(signatures[[i]],signature), logical(1)))
    if (!length(duplicate)) {
      unique_boxes[[length(unique_boxes)+1L]] <- box
      signatures[[length(signatures)+1L]] <- signature
    }
  }
  for (side in unique(vapply(unique_boxes,`[[`,"","side"))) {
    items <- Filter(function(b) b$side==side,unique_boxes)
    native <- styles$guides[[side]]
    if (!is.null(native)) for (i in which(native$layout$name=="guides"))
      items[[length(items)+1L]] <- list(grob=native$grobs[[i]])
    horizontal <- side %in% c("top","bottom")
    sizes <- lapply(items,function(item) ggchord_guide_size(item$grob))
    limit <- canvas[if(horizontal) 1L else 2L]-8
    if (any(vapply(sizes,function(s) s[if(horizontal) 1L else 2L]>limit,logical(1))))
      ggchord_stop("Fusion: a collected legend exceeds the canvas; revise its native guide layout")
    positions <- list(); row <- 0; used <- 0; depth <- 0; extent <- 0
    for (i in seq_along(items)) {
      s <- sizes[[i]]; along <- s[if(horizontal) 1L else 2L]; across <- s[if(horizontal) 2L else 1L]
      if (used>0 && used+along>limit) {row<-row+depth+4;used<-0;depth<-0}
      positions[[i]] <- c(used+along/2,row+across/2)
      used <- used+along+5; depth<-max(depth,across);extent<-max(extent,used-5)
    }
    dims <- if(horizontal) c(extent,row+depth) else c(row+depth,extent)
    grobs <- lapply(seq_along(items),function(i) {
      pos<-positions[[i]];if(!horizontal) pos<-rev(pos)
      grobTree(items[[i]]$grob,vp=viewport(x=pos[1],y=pos[2],
        width=sizes[[i]][1],height=sizes[[i]][2],default.units="mm"))
    })
    box <- items[[1]]$template %||% unique_boxes[[1]]$template
    box$grobs <- grobs
    box$layout <- data.frame(t=1,l=1,b=1,r=1,z=seq_along(grobs),clip="off",name=paste0("shared-",seq_along(grobs)))
    box$widths <- unit(dims[1],"mm");box$heights<-unit(dims[2],"mm")
    styles$guides[[side]] <- box;styles$sizes[side,]<-dims
  }
  styles$side <- if(length(styles$guides)==1L) names(styles$guides) else "multiple"
  styles$report$collected <- list(input=length(boxes),shared=length(unique_boxes),
    sources=vapply(boxes,`[[`,"","source"),mode="collect")
  styles
}

ggchord_fusion_style_layer <- function(layer) {
  d <- layer$data
  # The carrier's x is a stable relationship row key. It has no source or
  # display coordinate meaning and never leaves the style build.
  row_key <- make.unique(c(names(d), ".fusion_row"))[ncol(d) + 1L]
  d[[row_key]] <- seq_len(nrow(d))
  mapping <- layer$mapping %||% aes()
  for (role in intersect(layer$roles, names(mapping))) mapping[[role]] <- NULL
  mapping$x <- aes(x = .data[[row_key]])$x
  mapping$y <- aes(y = 0)$y
  params <- c(layer$params, layer$fixed[setdiff(names(layer$fixed), names(layer$params))])
  geom <- ggproto(NULL, layer$geom, draw_panel = function(...) nullGrob())
  carrier_layer <- ggplot2::layer(data = d, mapping = mapping, geom = geom,
    stat = "identity", position = "identity", inherit.aes = FALSE,
    show.legend = layer$show.legend, params = params)
  carrier_layer$ggchord_layer_id <- layer$id
  carrier_layer$ggchord_params <- list(id = layer$scale_id)
  carrier_layer$ggchord_type <- layer$kind

  carrier_layer
}

ggchord_fusion_check_style_mapping <- function(layer) {
  previous_plot <- last_plot()
  on.exit(set_last_plot(previous_plot), add = TRUE)
  carrier <- ggchord_fusion_style_layer(layer)
  plot <- ggplot()
  data <- carrier$setup_layer(carrier$data, plot)
  carrier$compute_aesthetics(data, plot)
  invisible(NULL)
}

# Train cross-panel styles through ggplot2; no second colour mapping engine.
ggchord_fusion_build_styles <- function(x) {
  p <- ggplot() + theme_void(base_size = 9) +
    theme(legend.position = "right", plot.margin = margin(0,0,0,0))
  for (layer in x$ribbons) p <- p + ggchord_fusion_style_layer(layer)
  explicit <- theme()
  for (style in x$styles) {
    p <- suppressMessages(p + style)
    if (inherits(style, "theme")) explicit <- explicit + style
  }
  # Capture user locks separately from the carrier's default right-hand guide.
  p$ggchord <- list(theme_explicit = explicit, palette_preset = x$palette_preset)
  if (!is.null(x$fusion_policy)) {
    p$ggchord$layout_policy <- x$fusion_policy
    p <- ggchord_apply_adaptive_styles(p)
  }
  scoped <- Filter(function(scope) any(scope$layers %in% names(x$ribbons)), x$scope_scales)
  for (scope in x$scope_scales) {
    if (any(!scope$layers %in% names(x$ribbons)))
      ggchord_stop("Fusion: scale_for() must select root ribbon layers")
  }
  specs <- ggchord_drawing_scopes(p, scoped)
  aliases <- lapply(specs, function(scope) list(selected = scope$selected,
    inverse = setNames(names(scope$aliases), scope$aliases)))
  p <- ggchord_apply_scoped_scales(p, specs)
  p <- ggchord_add_visual_scales(p)
  list(plot = p, built = ggplot_build(p), aliases = aliases)
}

ggchord_fusion_measure_styles <- function(x, styles) {
  p <- styles$plot
  built <- styles$built
  aliases <- styles$aliases
  report <- list(status = "fixed", render_policy = p$ggchord$adaptive_report)
  if (identical(x$fusion_policy$mode, "adaptive")) {
    canvas <- c(convertWidth(unit(1, "npc"), "mm", TRUE),
      convertHeight(unit(1, "npc"), "mm", TRUE))
    evaluate <- function(sizes) {
      budget <- ggchord_fusion_guide_budget(canvas, sizes)
      rects <- tryCatch(ggchord_fusion_rectangles(x, budget$area), error = function(e) NULL)
      if (min(budget$area) < 15 || is.null(rects)) return(list(penalty = 1e6, feasible = FALSE))
      dimensions <- lapply(rects, function(r) c(diff(r[1:2]), diff(r[3:4]))*budget$area)
      automatic <- vapply(x$panels, function(panel) is.null(panel$size) && is.null(panel$box), logical(1))
      small <- sum(vapply(dimensions[automatic], function(d) sum(pmax(15-d, 0)), numeric(1)))
      # Side boxes must coexist without meeting at the reserved corners.
      small <- small + sum(pmax(sizes[c("left", "right"), "height"]-budget$area[2], 0)) +
        sum(pmax(sizes[c("top", "bottom"), "width"]-budget$area[1], 0))
      list(penalty = small*1e4, feasible = small == 0)
    }
    joint <- ggchord_with_render_context(list(available_canvas_mm = canvas),
      ggchord_guide_joint_layout(built, p, evaluate = evaluate,
        positions = c("left", "right", "top", "bottom")))
    report <- joint$report
    if (!report$status %in% c("feasible", "no_guides")) ggchord_stop(
      "Fusion: cross-panel guides and panel allocations cannot fit the composition under the explicit constraints; enlarge the canvas or revise guide layout")
    built <- joint$built; table <- joint$table %||% ggplot_gtable(built)
  } else table <- ggplot_gtable(ggchord_fusion_hide_guides(built, x$ggchord_fusion_omit_guides))
  at <- which(grepl("^guide-box", table$layout$name) & !vapply(table$grobs, inherits, logical(1), "zeroGrob"))
  sides <- sub("guide-box-", "", table$layout$name[at])
  if (any(!sides %in% c("left", "right", "top", "bottom")))
    ggchord_stop("Fusion: cross-panel legends must be outside the panels")
  guides <- setNames(table$grobs[at], sides)
  sizes <- matrix(0, 4, 2, dimnames = list(c("left", "right", "top", "bottom"), c("width", "height")))
  for (side in sides) sizes[side, ] <- ggchord_guide_size(guides[[side]])
  styles <- lapply(seq_along(x$ribbons), function(i) {
    d <- built$data[[i]]
    for (alias in aliases) if (i %in% alias$selected) d <- ggchord_alias_aes(d, alias$inverse)
    d[match(seq_len(nrow(x$ribbons[[i]]$data)), d$x), , drop = FALSE]
  })
  list(data = styles, guides = guides, sizes = sizes, report = report,
    guide = if (length(at) == 1L) guides[[1L]] else nullGrob(),
    side = if (!length(at)) "none" else if (length(at) == 1L) sides else "multiple")
}


# Reserve all native outside guide boxes before measuring any child panel.
ggchord_fusion_guide_budget <- function(canvas, sizes) {
  horizontal <- c("left", "right"); vertical <- c("top", "bottom")
  widths <- ifelse(sizes[horizontal, "width"] > 0, sizes[horizontal, "width"]+4, 0)
  heights <- ifelse(sizes[vertical, "height"] > 0, sizes[vertical, "height"]+4, 0)
  list(area = canvas-c(sum(widths), sum(heights)), shift = c(widths[1L], heights[2L]))
}

ggchord_fusion_rectangles <- function(x, canvas) {
  layout <- x$fusion_areas
  ids <- vapply(names(x$panels), function(key) ggchord_fusion_path_ids(key)[1L], character(1))
  if (inherits(layout, "ggchord_fusion_grid"))
    layout <- ggchord_fusion_grid_areas(layout, ids, canvas)
  if (!setequal(ids, layout$areas$panel))
    ggchord_stop("Fusion: allocate exactly one area to each registered view",
      "; missing: ", paste(setdiff(ids, layout$areas$panel), collapse = ", "),
      "; unknown: ", paste(setdiff(layout$areas$panel, ids), collapse = ", "))
  areas <- layout$areas[match(ids, layout$areas$panel), , drop = FALSE]
  bounds <- as.matrix(areas[c("xmin", "xmax", "ymin", "ymax")])
  if (layout$units != "relative") {
    multiplier <- if (layout$units == "mm") 1 else 25.4
    bounds <- sweep(bounds * multiplier, 2L, rep(canvas, each = 2L), "/")
  }
  if (any(bounds > 1 + 1e-10)) ggchord_stop("Fusion: panel allocation lies outside the canvas")
  bounds[bounds > 1] <- 1
  z <- areas$z %||% seq_len(nrow(areas))
  rects <- setNames(lapply(seq_len(nrow(areas)), function(i) bounds[i, ]), names(x$panels))
  rects[order(z)]
}

# Partition backgrounds without deleting explicit fills or repainting foreground.
ggchord_fusion_partition <- function(g, background = FALSE, inherited = FALSE) {
  is_background <- inherited || grepl("^(plot[.]background|panel[.]background|background)([.-]|$)", g$name %||% "")
  if (!is.null(g$grobs)) {
    g$grobs <- lapply(g$grobs, ggchord_fusion_partition, background = background, inherited = is_background)
  } else if (!is.null(g$children) && length(g$children)) {
    children <- lapply(g$children, ggchord_fusion_partition, background = background, inherited = is_background)
    g$children <- do.call(gList, children)
  } else if (background != is_background) return(nullGrob(name = g$name))
  g
}

# Optional roles distinguish an omitted selector from an explicit invalid ID.
ggchord_fusion_selector <- function(row, name) {
  value <- row[[name]][[1L]]
  if (is.null(value) || (is.atomic(value) && length(value) == 1L && is.na(value)))
    return(NULL)
  value
}

ggchord_fusion_facet <- function(value, where) {
  if (is.null(value)) return(NULL)
  if (!is.list(value) || !length(value) || is.null(names(value)) ||
      anyNA(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value)) ||
      !all(vapply(value, function(x) is.atomic(x) && length(x) == 1L, logical(1))))
    ggchord_stop(where, "facet selectors must be named lists of facet values")
  value
}

ggchord_fusion_facet_match <- function(key, selector) {
  all(names(selector) %in% names(key)) && all(vapply(names(selector), function(name)
    identical(as.character(key[[name]]), as.character(selector[[name]])), logical(1)))
}

ggchord_fusion_endpoint_selectors <- function(row, prefix, where) {
  path <- ggchord_fusion_selector(row, paste0(prefix, "panel_path"))
  panel <- if (!is.null(path)) ggchord_fusion_check_path(path, where)
  instance <- ggchord_fusion_selector(row, paste0(prefix, "instance"))
  if (!is.null(instance) && (!is.character(instance) || length(instance) != 1L || !nzchar(instance)))
    ggchord_stop(where, "instance selector must be a nonempty string")
  facet <- ggchord_fusion_facet(ggchord_fusion_selector(row, paste0(prefix, "facet")), where)
  list(panel = panel, instance = instance, facet = facet)
}

ggchord_fusion_endpoint <- function(panels, row, prefix, layer_id, source_row, kind) {
  where <- paste0("Fusion ", layer_id, " row ", source_row, " ", prefix, ": ")
  selected <- ggchord_fusion_endpoint_selectors(row, prefix, where)
  point <- kind == "link"
  entity <- as.character(row[[paste0(prefix, "accver")]])
  start <- row[[paste0(prefix, if (point) "pos" else "start")]]
  end <- if (point) start else row[[paste0(prefix, "end")]]
  strand <- if (!point) ggchord_link_endpoint_strand(start, end,
    row[[paste0(prefix, "strand")]], where, paste0(prefix, "strand"))
  if (identical(strand, "-")) {
    lower <- min(start, end)
    start <- max(start, end)
    end <- lower
  }
  resolved <- ggchord_resolve_panel_interval(panels, entity, start, end,
    selected$panel, selected$instance, where, selected$facet, point)
  resolved$strand <- strand
  resolved
}

ggchord_resolve_panel_interval <- function(panels, entity, start, end,
    selected_panel, selected_instance, where, selected_facet = NULL, point = FALSE) {
  if (length(entity)!=1 || is.na(entity) || !nzchar(entity) ||
      !is.numeric(start) || !is.numeric(end) || any(!is.finite(c(start,end))) || min(start,end) < 1 ||
      any(c(start,end) != floor(c(start,end)))) ggchord_stop(where, "invalid accession or one-based interval")
  selected <- names(panels)
  if (!is.null(selected_panel)) {
    if (!selected_panel %in% names(panels)) ggchord_stop(where, "unknown leaf panel '", selected_panel, "'; available: ", paste(selected, collapse = ", "))
    selected <- selected_panel
  }
  candidates <- available <- list()
  for (id in selected) {
    layout <- panels[[id]]$layout
    if (isTRUE(layout$character)) {
      if (!length(selected_facet))
        candidates <- c(candidates, ggchord_alignment_interval(panels[[id]], id,
          entity, start, end, selected_instance, point = point))
      next
    }
    if (isTRUE(layout$dotplot) || is.null(layout) || is.null(layout$sequence_reference)) {
      if (!is.null(selected_panel)) ggchord_stop(where,
        if (isTRUE(layout$dotplot)) "Dotplot requires a separate dual-axis endpoint contract" else
          "the selected plot has no registered source coordinate provider")
      next
    }
    frames <- layout$panel_references %||% list(layout)
    for (facet in seq_along(frames)) {
      frame <- frames[[facet]]
      facet_at <- match(facet, as.integer(layout$facet_layout$PANEL))
      facet_key <- if (length(layout$facet_keys))
        as.list(layout$facet_layout[facet_at, layout$facet_keys, drop = FALSE]) else list()
      keys <- frame$seqs
      if (!length(keys)) next
      instances <- frame$display_instances
      match_instance <- match(keys, instances$key)
      entities <- instances$accver[match_instance] %||% keys
      instance_ids <- instances$instance_id[match_instance] %||% keys
      source <- ggchord_focus_source(panels[[id]]$seq_data)
      at <- match(entities, panels[[id]]$seq_data$accver)
      originals <- source$accver[at]
      origins <- source$origin[at]
      tab <- frame$sequence_layout
      if (is.data.frame(tab) && all(c("source_accver", "source_origin") %in% names(tab))) {
        tab <- tab[match(keys, tab$accver), , drop = FALSE]
        originals[!is.na(tab$source_accver)] <- tab$source_accver[!is.na(tab$source_accver)]
        origins[!is.na(tab$source_origin)] <- tab$source_origin[!is.na(tab$source_origin)]
      }
      for (j in which(!is.na(originals))) {
        key <- keys[j]
        origin <- origins[j]
        size <- unname(frame$sequence_reference$lens[key])
        if (!is.finite(origin) || !is.finite(size)) next
        endpoint <- list(panel = id, facet = facet, facet_key = facet_key, key = key,
          instance = as.character(instance_ids[j]), entity = originals[j],
          origin = origin, visible_end = origin + size - 1,
          reference_length = source$length[at[j]],
          reference_namespace = panels[[id]]$seq_data[["reference_namespace"]][at[j]],
          start = start, end = end, point = point,
          unit = frame$sequence_unit %||% "bp")
        available[[length(available) + 1L]] <- endpoint
        if (originals[j] == entity && min(start, end) >= origin && max(start, end) <= endpoint$visible_end &&
            (is.null(selected_instance) || identical(endpoint$instance, selected_instance)) &&
            ggchord_fusion_facet_match(facet_key, selected_facet))
          candidates[[length(candidates) + 1L]] <- endpoint
      }
    }
  }
  if (!length(candidates)) {
    matches <- vapply(head(available, 8L), ggchord_fusion_endpoint_label, character(1), panels = panels)
    ggchord_stop(where, "no fully visible matching endpoint for ", entity, " [", start, ", ", end,
      "]; candidates: ", paste(matches, collapse = "; "))
  }

  if (length(candidates) > 1L) {
    matches <- vapply(head(candidates, 8L), ggchord_fusion_endpoint_label, character(1),
      panels = panels)
    action <- "Specify the panel path, facet values or instance ID to identify one match."
    ggchord_stop(where, "ambiguous endpoint for source '", entity, "' [", start, ", ", end,
      "]. Matching candidates: ", paste(matches, collapse = "; "),
      if (length(candidates) > length(matches)) paste0("; +", length(candidates) - length(matches), " more"),
      ". ", action)
  }
  candidates[[1L]]
}

ggchord_fusion_endpoint_label <- function(candidate, panels) {
  layout <- panels[[candidate$panel]]$layout
  path <- paste(encodeString(ggchord_fusion_path_ids(candidate$panel), quote = '"'), collapse = "/")
  label <- paste0(path, " instance ", encodeString(candidate$instance, quote = '"'))
  if (!is.null(candidate$visible_end))
    label <- paste0(label, " ", candidate$entity, " [", candidate$origin, ", ", candidate$visible_end, "]")
  if (is.null(layout$facet_layout)) return(label)

  at <- match(candidate$facet, as.integer(layout$facet_layout$PANEL))
  facet <- layout$facet_layout[at, layout$facet_keys, drop = FALSE]
  values <- vapply(facet, function(x) encodeString(as.character(x), quote = '"'), character(1))
  paste0(label, " (", paste(c(paste0("PANEL=", candidate$facet),
    if (length(values)) paste0(names(values), "=", values)), collapse = ", "), ")")
}

ggchord_fusion_reference_match <- function(q, s, caller) {
  for (field in c("reference_namespace", "reference_length")) {
    if (length(q[[field]]) && length(s[[field]]) && isTRUE(q[[field]] != s[[field]]))
      ggchord_stop(caller, ": conflicting ", field, " for source '", q$entity, "'")
  }
}

ggchord_fusion_front <- function(panel, endpoint) {
  if (isTRUE(panel$layout$character)) return(ggchord_alignment_front(panel, endpoint))
  facet <- endpoint$facet %||% 1L
  layout <- panel$layout$panel_references[[facet]] %||% panel$layout
  boundary <- ggchord_interval_front(endpoint$start, endpoint$end) - endpoint$origin + 1
  pos <- if (isTRUE(endpoint$point)) endpoint$start - endpoint$origin + 1 else
    seq(boundary[1L], boundary[2L], length.out = 40)
  frame <- ggchord_panel_backbone_frame(layout, endpoint$key, pos)
  xy <- frame$xy
  normal <- frame$normal
  project <- function(points) as.data.frame(ggchord_query_transform(panel, points, facet))
  projected <- project(xy)
  if (any(!is.finite(c(projected$x,projected$y))) || any(projected$x < -1e-7 | projected$x > 1+1e-7 | projected$y < -1e-7 | projected$y > 1+1e-7))
    ggchord_stop("Fusion: source front is outside the visible coordinate limits in panel '", endpoint$panel, "'")
  bounds <- panel$context$panel_bounds_mm
  facet_context <- panel$context$facet_contexts[[facet]]
  if (!is.null(facet_context))
    bounds <- facet_context$panel_bounds_mm - (panel$context$device_panel_bounds_mm - bounds)
  normal_tip <- project(xy+normal)
  normal <- cbind((normal_tip$x-projected$x)*diff(bounds[1:2]),
    (normal_tip$y-projected$y)*diff(bounds[3:4]))
  normal <- normal / sqrt(rowSums(normal^2))
  data.frame(position_boundary = if (isTRUE(endpoint$point)) NA_real_ else pos + endpoint$origin - 1,
    source_position = if (isTRUE(endpoint$point)) endpoint$start else NA_real_, local_x=xy[,1],local_y=xy[,2],
    x=bounds[1]+projected$x*diff(bounds[1:2]),y=bounds[3]+projected$y*diff(bounds[3:4]),
    normal_x=normal[,1],normal_y=normal[,2],tangent_x=-normal[,2],tangent_y=normal[,1],
    facet_panel = facet)
}

ggchord_fusion_point_path <- function(q, s, route, curvature) {
  start <- unlist(q[c("x", "y")], use.names = FALSE)
  end <- unlist(s[c("x", "y")], use.names = FALSE)
  toward <- end - start
  control <- function(front, direction) {
    normal <- unlist(front[c("normal_x", "normal_y")], use.names = FALSE)
    if (sum(normal * direction) < 0) -normal else normal
  }
  nq <- control(q, toward)
  ns <- control(s, -toward)
  distance <- sqrt(sum(toward^2)) * curvature
  path <- if (route == "straight") rbind(start, end) else
    bezier_pts(start, end, start + nq * distance, end + ns * distance, n = 60)
  attr(path, "endpoint_normals") <- rbind(q = nq, s = ns)
  path
}

ggchord_fusion_path <- function(q, s, curvature) {
  qc <- colMeans(q[c("x","y")]); sc <- colMeans(s[c("x","y")])
  normal <- function(front, toward) {
    middle <- floor(nrow(front)/2)
    tangent <- unlist(front[middle+1L,c("x","y")])-unlist(front[middle,c("x","y")])
    if (sum(tangent^2) == 0) ggchord_stop("Fusion: degenerate interval front")
    n <- c(-tangent[2],tangent[1])/sqrt(sum(tangent^2))
    if (sum(n*toward)<0) n <- -n
    n
  }
  nq <- normal(q,sc-qc); ns <- normal(s,qc-sc)
  # Near a tangential approach, centre-facing signs are unstable: opposite
  # transverse handles produce an artificial S bend. Keep the well-facing
  # endpoint and choose the other valid backbone normal for a coherent arch.
  axis <- sc - qc
  distance_between <- sqrt(sum(axis^2))
  axis <- if (distance_between > 0) axis / distance_between else c(1, 0)
  transverse <- c(-axis[2],axis[1])
  facing <- c(abs(sum(nq*axis)),abs(sum(ns*axis)))
  if (sum(nq*transverse)*sum(ns*transverse)<0 && min(facing)<.15) {
    if (facing[1]<facing[2]) nq <- -nq else ns <- -ns
  }
  distance <- distance_between * curvature
  edge <- function(a,b) bezier_pts(a,b,a+nq*distance,b+ns*distance,n=60)
  a <- as.matrix(q[c("x","y")]); b <- as.matrix(s[c("x","y")])
  edge1 <- edge(a[1,],b[1,]); edge2 <- edge(a[nrow(a),],b[nrow(b),])
  path <- rbind(a,edge2,b[nrow(b):1,,drop=FALSE],edge1[nrow(edge1):1,,drop=FALSE])
  attr(path,"endpoint_normals") <- rbind(q=nq,s=ns)
  path
}

#' @export
print.ggchord_fusion <- function(x, ...) {
  ggchord_local_console_errors()
  set_last_plot(x)
  grid::grid.newpage()
  grid.draw(x)
  invisible(x)
}

#' @export
#' @importFrom grid grid.draw
grid.draw.ggchord_fusion <- function(x, recording = TRUE) {
  ggchord_local_console_errors()
  # Keep the composition in the display list, rather than fixed-mm primitives.
  # Recorded plots may be replayed on a device with a different size or DPI.
  # Rebuild there so panels, interval fronts and leaders share the final metric.
  grob <- grid::recordGrob(
    grid.draw(ggchord_build_fusion(plot)$grob, recording = FALSE),
    list(plot = x, ggchord_build_fusion = ggchord_build_fusion)
  )
  grid.draw(grob, recording = recording)
  invisible(x)
}

#' @export
ggplot_build.ggchord_fusion <- function(plot, ...) {
  ggchord_local_console_errors()
  ggchord_stop("Fusion is a composition: use grid.draw(), ggsave() or export_ggchord_layout()")
}

# The default cross-panel ruler uses its final orientation and real text metrics.
ggchord_fusion_identity_guide <- function(plot, styles, direction = NULL) {
  theme <- plot$theme
  explicit <- list()
  for (style in styles) if (inherits(style,"theme")) {
    theme <- theme + style
    for (field in c("legend.key.width","legend.key.height","legend.key.size"))
      if (!is.null(style[[field]])) explicit[[field]] <- style[[field]]
  }
  horizontal <- if (!is.null(direction)) direction == "horizontal" else
    (theme$legend.direction %||% if ((theme$legend.position %||% "right") %in% c("top", "bottom")) "horizontal" else "vertical") == "horizontal"
  text <- ggplot2::calc_element("legend.text",ggplot2::theme_get()+theme)
  values <- c(0,50,80,90,95,100)
  dimensions <- if(inherits(text,"element_blank")) rep(0,length(values)) else
    vapply(as.character(values),function(label) {
      grob <- textGrob(label,gp=gpar(fontsize=text$size %||% 9,
        fontfamily=text$family %||% "",fontface=text$face %||% "plain"))
      if(horizontal) convertWidth(grobWidth(grob),"mm",TRUE) else
        convertHeight(grobHeight(grob),"mm",TRUE)
    },numeric(1))
  required <- max(((head(dimensions,-1)+tail(dimensions,-1))/2+1.2)/(diff(values)/100))
  field <- if(horizontal) "legend.key.width" else "legend.key.height"
  thickness <- if(horizontal) "legend.key.height" else "legend.key.width"
  length <- explicit[[field]] %||% explicit$legend.key.size
  if(is.null(length)) length <- unit(max(required,if(horizontal) 100 else 85),"mm") else
    if(is.unit(length) && convertUnit(length,"mm",valueOnly=TRUE)<required)
      warning("Fusion Identity labels exceed the explicitly fixed guide length; increase its length or reduce its explicit text size",call.=FALSE)
  width <- explicit[[thickness]] %||% explicit$legend.key.size %||% ggchord_identity_bar_width()
  guide_theme <- do.call(ggplot2::theme,setNames(list(length,width),c(field,thickness)))
  result <- ggchord_identity_title_gap(guide_coloursteps(even.steps=FALSE,show.limits=TRUE,
    direction=if(horizontal) "horizontal" else "vertical",theme=guide_theme))
  attr(result, "ggchord_fusion_default") <- TRUE
  result
}
