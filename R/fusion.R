# Composition owns device-local panels; source plots remain independent.

#' Compose independent plots
#' @param panels Uniquely named list of ggchord plots, Cartesian ggplots or
#'   nested compositions. Source endpoints require a registered source map.
#' @param layout Panel allocation from [layout_fusion()] or [layout_fusion_grid()].
#' @param guides Keep child guides local, or collect semantically equivalent guides.
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
ggchord_fusion <- function(panels, layout = layout_fusion_grid(), guides = "keep") {
  ids <- names(panels)
  if (!is.list(panels) || !length(panels) || length(ids) != length(panels) ||
      anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    ggchord_stop("ggchord_fusion(): panels must be a uniquely named list")
  if (!all(vapply(panels, inherits, logical(1), "ggplot")))
    ggchord_stop("ggchord_fusion(): each panel must be an independent plot or Fusion")
  if (inherits(layout, "ggchord_fusion_grid"))
    layout <- ggchord_fusion_grid_areas(layout, ids)
  if (!inherits(layout, "ggchord_fusion_layout"))
    ggchord_stop("ggchord_fusion(): layout must come from layout_fusion() or layout_fusion_grid()")
  if (!setequal(ids, layout$areas$panel))
    ggchord_stop("ggchord_fusion(): allocate exactly one area to each named panel")
  layout$areas <- layout$areas[match(ids, layout$areas$panel), , drop = FALSE]
  keys <- vapply(ids, ggchord_fusion_path_key, character(1))
  plot <- ggplot2::ggplot()
  plot$panels <- stats::setNames(lapply(panels, function(child) list(plot = child)), keys)
  plot$fusion_areas <- layout
  plot$ribbons <- list()
  plot$zooms <- list()
  plot$styles <- list()
  plot$scope_scales <- list()
  plot$allow_overlap <- layout$allow_overlap
  plot$fusion_guides <- match.arg(guides, c("keep", "collect"))
  class(plot) <- c("ggchord_fusion", class(plot))
  ggplot2::set_last_plot(plot)
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
  caller <- "layout_fusion()"
  fields <- c("xmin", "xmax", "ymin", "ymax")
  ggchord_require_columns(areas, c("panel", fields), caller)
  units <- match.arg(units, c("relative", "mm", "in"))
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
layout_fusion_grid <- function(design = NULL, widths = 1, heights = 1) {
  if (!is.null(design) && (!is.matrix(design) || !is.character(design) ||
      !length(design) || all(is.na(design)) || any(!is.na(design) & !nzchar(design))))
    ggchord_stop("layout_fusion_grid(): design must be a nonempty character matrix")
  for (weights in list(widths, heights)) {
    if (!is.numeric(weights) || !length(weights) || any(!is.finite(weights) | weights <= 0))
      ggchord_stop("layout_fusion_grid(): widths and heights must be positive finite weights")
  }
  spec <- list(design = design, widths = widths, heights = heights)
  if (!is.null(design)) return(ggchord_fusion_grid_areas(spec, unique(c(t(design))[!is.na(c(t(design)))])))
  structure(spec, class = "ggchord_fusion_grid")
}

ggchord_fusion_grid_areas <- function(spec, ids) {
  design <- spec$design
  if (is.null(design)) {
    columns <- ceiling(sqrt(length(ids)))
    rows <- ceiling(length(ids) / columns)
    design <- matrix(c(ids, rep(NA_character_, rows * columns - length(ids))),
      nrow = rows, ncol = columns, byrow = TRUE)
  }
  edges <- function(weights, n) {
    if (!length(weights) %in% c(1L, n))
      ggchord_stop("layout_fusion_grid(): weights must match the grid dimensions")
    weights <- rep_len(weights, n) / max(weights)
    edges <- c(0, cumsum(weights) / sum(weights))
    edges[n + 1L] <- 1
    edges
  }
  x <- edges(spec$widths, ncol(design))
  y <- 1 - edges(spec$heights, nrow(design))
  rows <- lapply(ids, function(id) {
    cells <- which(design == id, arr.ind = TRUE)
    if (!nrow(cells)) ggchord_stop("layout_fusion_grid(): missing panel '", id, "'")
    left <- min(cells[, 2]); right <- max(cells[, 2])
    top <- min(cells[, 1]); bottom <- max(cells[, 1])
    if (nrow(cells) != (right - left + 1L) * (bottom - top + 1L))
      ggchord_stop("layout_fusion_grid(): panel '", id, "' must occupy a rectangle")
    data.frame(panel = id, xmin = x[left], xmax = x[right + 1L],
      ymin = y[bottom + 1L], ymax = y[top])
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
  if (!is.character(path) || !length(path) || anyNA(path) || any(!nzchar(path)))
    ggchord_stop(label, ": panel paths must contain nonempty string IDs")
  ggchord_fusion_path_key(path)
}

#' Connect real source intervals across panels
#' @param mapping Optional aesthetic mappings, including interval and path roles.
#' @param data Relationship table with q/s endpoints and qpanel_path/spanel_path.
#' @param route Curve or straight ribbon routing.
#' @param curvature Curve strength in the interval `[0, 1]`.
#' @param ... Fixed visual attributes, lineend or key_glyph.
#' @param na.rm Remove incomplete relationships without a warning.
#' @param show.legend Whether the layer contributes to composition guides.
#' @param id Optional root layer ID.
#' @return A component for [ggchord_fusion()].
#' @export
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
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' cross <- links[1, ]
#' cross$qpanel_path <- "overview"
#' cross$spanel_path <- "detail"
#' f <- ggchord_fusion(list(overview = p_basic, detail = p_basic))
#' f + fusion_ribbon(aes(fill = pident), data = cross, id = "cross") +
#'   scale_fill_gradient(limits = c(0, 100))
fusion_ribbon <- function(mapping = NULL, data = NULL, route = "curve",
    curvature = 0.5, ..., na.rm = FALSE, show.legend = NA, id = NULL) {
  caller <- "fusion_ribbon()"
  roles <- c("qpanel_path", "spanel_path", "qaccver", "saccver",
    "qstart", "qend", "sstart", "send", "qinstance", "sinstance", "qstrand", "sstrand")
  if (!is.data.frame(data)) ggchord_stop(caller, ": provide a relationship data frame")
  if (!is.null(mapping) && !inherits(mapping, c("ggplot2::mapping", "uneval")))
    ggchord_stop(caller, ": mapping must be aes()")
  styles <- c("fill", "colour", "alpha", "linewidth", "linetype")
  if (length(setdiff(names(mapping), c(roles, styles, "group"))))
    ggchord_stop(caller, ": unknown mapped role or aesthetic")
  fixed <- list(...)
  if (length(fixed) && (length(names(fixed)) != length(fixed) || anyDuplicated(names(fixed)) ||
      any(!names(fixed) %in% c(styles, "lineend", "key_glyph"))))
    ggchord_stop(caller, ": supply only named visual attributes, lineend or key_glyph")
  if (any(lengths(fixed[intersect(names(fixed), styles)]) != 1L))
    ggchord_stop(caller, ": fixed visual attributes must be scalar")
  for (name in intersect(names(fixed), styles)) mapping[[name]] <- NULL
  route <- match.arg(route, c("curve", "straight"))
  if (!is.numeric(curvature) || length(curvature) != 1L || !is.finite(curvature) ||
      curvature < 0 || curvature > 1 || (route == "straight" && curvature != 0.5))
    ggchord_stop(caller, ": curvature must be in `[0,1]` and only applies to curve routing")
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm))
    ggchord_stop(caller, ": na.rm must be TRUE or FALSE")
  if (!is.logical(show.legend) || !length(show.legend) ||
      (length(show.legend) > 1L && (is.null(names(show.legend)) || anyNA(names(show.legend)) ||
        anyDuplicated(names(show.legend)) || any(!names(show.legend) %in% styles))))
    ggchord_stop(caller, ": show.legend must be a logical scalar or named vector")
  if (!is.null(fixed$lineend) && (length(fixed$lineend) != 1L ||
      !fixed$lineend %in% c("butt", "round", "square")))
    ggchord_stop(caller, ": lineend must be butt, round or square")
  required <- roles[seq_len(8L)]
  ggchord_require_columns(data, setdiff(required, names(mapping)), caller)
  structure(list(data = data, mapping = mapping, roles = roles, required = required,
    route = route, curvature = curvature, fixed = fixed, na.rm = na.rm,
    show.legend = show.legend, id = id), class = "ggchord_fusion_ribbon")
}

ggchord_fusion_input <- function(layer) {
  data <- as.data.frame(layer$data)
  for (role in intersect(names(layer$mapping), layer$roles)) {
    value <- rlang::eval_tidy(layer$mapping[[role]], data = layer$data)
    if (length(value) == 1L) value <- rep(value, nrow(data))
    if (length(value) != nrow(data))
      ggchord_stop("Fusion layer '", layer$id, "': ", role, " must have one value per row")
    data[[role]] <- value
  }
  data
}

#' @export
`+.ggchord_fusion` <- function(e1, e2) {
  if (inherits(e2, c("ggchord_fusion_ribbon", "ggchord_fusion_zoom"))) {
    family <- if (inherits(e2, "ggchord_fusion_ribbon")) "ribbons" else "zooms"
    prefix <- if (family == "ribbons") "ribbon" else "zoom"
    ids <- c(names(e1$ribbons), names(e1$zooms))
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
    e1[[family]][[e2$id]] <- unclass(e2)
  } else if (inherits(e2, "ggchord_scoped_scale")) {
    e1$scope_scales <- c(e1$scope_scales, list(e2))
  } else if (inherits(e2, "ggchord_layout_policy")) {
    e1$fusion_policy <- unclass(e2)
  } else if (inherits(e2, c("Scale", "theme", "Guides", "ggplot2::labels"))) {
    if (inherits(e2, "Scale") && any(e2$aesthetics %in% c("x", "y")))
      ggchord_stop("Fusion: root scales map visual aesthetics")
    e1$styles <- c(e1$styles, list(e2))
  } else ggchord_stop("Fusion accepts ribbons, zooms, render policies and visual plot components")
  ggplot2::set_last_plot(e1)
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
    list(aesthetics = kept, type = class(scale)[1L], limits = scale$get_limits(),
      transform = if (!scale$is_discrete()) scale$get_transformation()$name,
      breaks = scale$get_breaks(), labels = scale$get_labels())
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
      fields <- intersect(all.vars(rlang::get_expr(value)), names(data))
      units <- lapply(data[fields], function(column) list(units = attr(column, "units"),
        unit = attr(column, "unit"), measure = attr(column, "measure")))
      list(expression = rlang::as_label(value), columns = units)
    })
    list(mapping = values, unit = layer$ggchord_params$unit,
      measure = layer$ggchord_params$measure)
  })
  variables <- unique(Filter(Negate(is.null), variables))
  scopes <- Filter(Negate(is.null), details)
  list(title = guide$title, variables = variables, scales = scopes,
    unresolved = if (!length(variables) || !length(scopes)) source)
}

ggchord_fusion_collect_guides <- function(x) {
  boxes <- list()
  walk <- function(plot, path) {
    if (inherits(plot, "ggchord_fusion")) {
      for (id in names(plot$panels)) plot$panels[[id]]$plot <-
        walk(plot$panels[[id]]$plot, paste(path, id, sep = "/"))
      return(plot)
    }
    if (inherits(plot, "ggchord"))
      plot <- ggchord_render_track_source(ggchord_render_inputs(plot))
    built <- ggplot2::ggplot_build(plot)
    table <- ggplot2::ggplot_gtable(built)
    params <- built$plot$guides$params
    positions <- vapply(params, function(guide)
      guide$position %||% built$plot$theme$legend.position %||% "right", character(1))
    at <- which(grepl("^guide-box", table$layout$name) &
      !vapply(table$grobs, inherits, logical(1), "zeroGrob"))
    collected <- character()
    incomplete <- FALSE
    for (i in at) {
      side <- sub("^guide-box-", "", table$layout$name[i])
      if (!side %in% c("top", "bottom", "left", "right")) next
      box <- table$grobs[[i]]
      ids <- which(box$layout$name == "guides")
      candidates <- params[positions == side]
      distinct <- length(candidates) == length(ids)
      for (j in seq_along(ids)) {
        semantics <- if (distinct) ggchord_fusion_guide_semantics(built, candidates[[j]], path) else
          list(unresolved = path)
        boxes[[length(boxes) + 1L]] <<- list(side = side, grob = box$grobs[[ids[j]]],
          template = box, source = path, semantics = semantics)
      }
      incomplete <- incomplete || !distinct
      collected <- union(collected, unlist(lapply(candidates, `[[`, "aesthetic")))
    }
    if (!length(at) || (!length(collected) && !incomplete)) return(plot)
    if (!any(positions %in% "inside"))
      return(plot + ggplot2::theme(legend.position = "none"))
    if (!length(collected)) return(plot)
    plot + do.call(ggplot2::guides, stats::setNames(rep(list("none"), length(collected)), collected))
  }
  for (id in names(x$panels)) x$panels[[id]]$plot <- walk(x$panels[[id]]$plot, id)
  list(plot = x, boxes = boxes)
}

# Native guides keep their physical text/key sizes. Only boxes are arranged.
ggchord_fusion_shared_guides <- function(styles, boxes, x, canvas) {
  if (!length(boxes)) return(styles)
  hidden <- FALSE
  for (style in x$styles) if (inherits(style,"theme") && !is.null(style$legend.position))
    hidden <- identical(style$legend.position,"none")
  if (hidden) return(styles)
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
      grid::grobTree(items[[i]]$grob,vp=grid::viewport(x=pos[1],y=pos[2],
        width=sizes[[i]][1],height=sizes[[i]][2],default.units="mm"))
    })
    box <- items[[1]]$template %||% unique_boxes[[1]]$template
    box$grobs <- grobs
    box$layout <- data.frame(t=1,l=1,b=1,r=1,z=seq_along(grobs),clip="off",name=paste0("shared-",seq_along(grobs)))
    box$widths <- grid::unit(dims[1],"mm");box$heights<-grid::unit(dims[2],"mm")
    styles$guides[[side]] <- box;styles$sizes[side,]<-dims
  }
  styles$side <- if(length(styles$guides)==1L) names(styles$guides) else "multiple"
  styles$report$collected <- list(input=length(boxes),shared=length(unique_boxes),
    sources=vapply(boxes,`[[`,"","source"),mode="collect")
  styles
}

# Train cross-panel styles through ggplot2; no second colour mapping engine.
ggchord_fusion_styles <- function(x) {
  p <- ggplot2::ggplot() + ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(legend.position = "right", plot.margin = ggplot2::margin(0,0,0,0))
  for (layer in x$ribbons) {
    input <- layer$data
    d <- input[rep(seq_len(nrow(input)), each = 3L), , drop = FALSE]
    d$.fusion_x <- rep(c(0, 1, 0), nrow(input))
    d$.fusion_y <- rep(c(0, 0, 1), nrow(input))
    d$.fusion_row <- rep(seq_len(nrow(input)), each = 3L)
    mapping <- layer$mapping %||% ggplot2::aes()
    for (role in intersect(c(layer$roles, "group"), names(mapping))) mapping[[role]] <- NULL
    mapping$x <- ggplot2::aes(x = .fusion_x)$x
    mapping$y <- ggplot2::aes(y = .fusion_y)$y
    mapping$group <- ggplot2::aes(group = .fusion_row)$group
    defaults <- list(fill = "steelblue", colour = NA, alpha = .35, linewidth = .25, linetype = 1)
    fixed <- layer$fixed
    for (name in setdiff(names(defaults), c(names(mapping), names(fixed))))
      fixed[[name]] <- defaults[[name]]
    carrier_layer <- do.call(ggplot2::geom_polygon, c(list(
      data = d, mapping = mapping, inherit.aes = FALSE, show.legend = layer$show.legend), fixed))
    carrier_layer$ggchord_layer_id <- layer$id
    carrier_layer$ggchord_type <- "fusion_ribbon"
    p <- p + carrier_layer
  }
  explicit <- ggplot2::theme()
  for (style in x$styles) {
    p <- suppressMessages(p + style)
    if (inherits(style, "theme")) explicit <- explicit + style
  }
  # Capture user locks separately from the carrier's default right-hand guide.
  p$ggchord <- list(theme_explicit = explicit)
  if (!is.null(x$fusion_policy)) {
    p$ggchord$layout_policy <- x$fusion_policy
    p <- ggchord_apply_adaptive_styles(p)
  }
  scoped <- Filter(function(scope) any(scope$layers %in% names(x$ribbons)), x$scope_scales)
  for (scope in x$scope_scales) {
    if (any(!scope$layers %in% names(x$ribbons)))
      ggchord_stop("Fusion: scale_for() must select root ribbon layers")
  }
  aliases <- lapply(ggchord_scale_scopes(p, scoped), function(scope) list(selected = scope$selected,
    inverse = stats::setNames(names(scope$aliases), scope$aliases)))
  if (length(scoped)) p <- ggchord_apply_scoped_scales(p, scoped)
  p <- ggchord_add_identity_scales(p)
  built <- ggplot2::ggplot_build(p)
  report <- list(status = "fixed", render_policy = p$ggchord$adaptive_report)
  if (identical(x$fusion_policy$mode, "adaptive")) {
    canvas <- c(grid::convertWidth(grid::unit(1, "npc"), "mm", TRUE),
      grid::convertHeight(grid::unit(1, "npc"), "mm", TRUE))
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
    built <- joint$built; table <- joint$table %||% ggplot2::ggplot_gtable(built)
  } else table <- ggplot2::ggplot_gtable(built)
  at <- which(grepl("^guide-box", table$layout$name) & !vapply(table$grobs, inherits, logical(1), "zeroGrob"))
  sides <- sub("guide-box-", "", table$layout$name[at])
  if (any(!sides %in% c("left", "right", "top", "bottom")))
    ggchord_stop("Fusion: cross-panel legends must be outside the panels")
  guides <- stats::setNames(table$grobs[at], sides)
  sizes <- matrix(0, 4, 2, dimnames = list(c("left", "right", "top", "bottom"), c("width", "height")))
  for (side in sides) sizes[side, ] <- ggchord_guide_size(guides[[side]])
  styles <- lapply(seq_along(built$data), function(i) {
    d <- built$data[[i]]
    for (alias in aliases) if (i %in% alias$selected) d <- ggchord_alias_aes(d, alias$inverse)
    d[seq.int(1L, length.out = nrow(d) / 3L, by = 3L), , drop = FALSE]
  })
  list(data = styles, guides = guides, sizes = sizes, report = report,
    guide = if (length(at) == 1L) guides[[1L]] else grid::nullGrob(),
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
  areas <- layout$areas
  bounds <- as.matrix(areas[c("xmin", "xmax", "ymin", "ymax")])
  if (layout$units != "relative") {
    multiplier <- if (layout$units == "mm") 1 else 25.4
    bounds <- sweep(bounds * multiplier, 2L, rep(canvas, each = 2L), "/")
  }
  if (any(bounds > 1 + 1e-10)) ggchord_stop("Fusion: panel allocation lies outside the canvas")
  bounds[bounds > 1] <- 1
  z <- areas$z %||% seq_len(nrow(areas))
  rects <- stats::setNames(lapply(seq_len(nrow(areas)), function(i) bounds[i, ]), names(x$panels))
  rects[order(z)]
}

# Partition backgrounds without deleting explicit fills or repainting foreground.
ggchord_fusion_partition <- function(g, background = FALSE, inherited = FALSE) {
  is_background <- inherited || grepl("^(plot[.]background|panel[.]background|background)([.-]|$)", g$name %||% "")
  if (!is.null(g$grobs)) {
    g$grobs <- lapply(g$grobs, ggchord_fusion_partition, background = background, inherited = is_background)
  } else if (!is.null(g$children) && length(g$children)) {
    children <- lapply(g$children, ggchord_fusion_partition, background = background, inherited = is_background)
    g$children <- do.call(grid::gList, children)
  } else if (background != is_background) return(grid::nullGrob(name = g$name))
  g
}

# Resolve entity + instance + fully visible source window before projection.
ggchord_fusion_endpoint <- function(panels, row, prefix, layer_id, source_row) {
  entity <- as.character(row[[paste0(prefix,"accver")]])
  selected_panel <- ggchord_fusion_check_path(
    row[[paste0(prefix, "panel_path")]][[1L]], "fusion_ribbon()")
  selected_instance <- row[[paste0(prefix, "instance")]]
  if (length(selected_instance) && is.na(selected_instance)) selected_instance <- NULL
  start <- row[[paste0(prefix, "start")]]
  end <- row[[paste0(prefix, "end")]]
  strand <- row[[paste0(prefix, "strand")]]
  if (length(strand) && !is.na(strand)) {
    if (!strand %in% c("+", "-", ".")) ggchord_stop("fusion_ribbon(): invalid strand")
    ordered <- sort(c(start, end), decreasing = strand == "-")
    start <- ordered[1L]
    end <- ordered[2L]
  }
  where <- paste0("Fusion ", layer_id, " row ", source_row, " ", prefix, ": ")
  ggchord_resolve_panel_interval(panels, entity, start, end, selected_panel, selected_instance, where)
}

# Shared locus resolver: decorative focus regions are not HSP relationships.
ggchord_resolve_panel_interval <- function(panels, entity, start, end,
    selected_panel, selected_instance, where) {
  if (length(entity)!=1 || is.na(entity) || !nzchar(entity) ||
      !is.numeric(start) || !is.numeric(end) || any(!is.finite(c(start,end))) || min(start,end) < 1 ||
      any(c(start,end) != floor(c(start,end)))) ggchord_stop(where, "invalid accession or one-based interval")
  selected <- names(panels)
  if (!is.null(selected_panel)) {
    if (!selected_panel %in% names(panels)) ggchord_stop(where, "unknown leaf panel '", selected_panel, "'")
    selected <- selected_panel
  }
  candidates <- list()
  for (id in selected) {
    layout <- panels[[id]]$layout
    if (isTRUE(layout$character)) {
      candidates <- c(candidates, ggchord_alignment_interval(panels[[id]], id, entity, start, end, selected_instance))
      next
    }
    if (is.null(layout) || isTRUE(layout$dotplot) || is.null(layout$sequence_reference)) next
    if (length(layout$panel_references) > 1L) ggchord_stop(where,
      "a faceted panel does not have a unique interval projection")
    keys <- layout$seqs
    instances <- layout$display_instances
    entities <- if (is.null(instances)) keys else instances$accver[match(keys,instances$key)]
    instance_ids <- if (is.null(instances)) keys else instances$instance_id[match(keys,instances$key)]
    seq_data <- panels[[id]]$seq_data
    original_entities <- entities
    source_ids <- seq_data[[".source_accver"]] %||% seq_data[["source_accver"]]
    if (length(source_ids)) original_entities <- as.character(source_ids[match(entities, seq_data$accver)])
    tab <- layout$sequence_layout
    if (is.data.frame(tab) && "source_accver" %in% names(tab)) {
      sources <- as.character(tab$source_accver[match(keys, tab$accver)])
      original_entities[!is.na(sources)] <- sources[!is.na(sources)]
    }
    for (at in which(original_entities == entity)) {
      if (!is.null(selected_instance) && (is.na(selected_instance) || instance_ids[at] != selected_instance)) next
      key <- keys[at]
      origin <- 1
      source_origins <- seq_data[[".source_origin"]] %||% seq_data[["source_origin"]]
      if (length(source_origins)) origin <- source_origins[match(entities[at], seq_data$accver)]
      if (is.data.frame(layout$sequence_layout) && "source_origin" %in% names(layout$sequence_layout))
        origin <- layout$sequence_layout$source_origin[match(key, layout$sequence_layout$accver)]
      length <- unname(layout$sequence_reference$lens[key])
      if (!is.finite(origin) || !is.finite(length) ||
          min(start, end) < origin || max(start, end) > origin + length - 1) next
      candidates[[length(candidates)+1L]] <- list(panel=id, key=key, instance=instance_ids[at],
        entity=entity, origin=origin, start=start, end=end, unit=layout$sequence_unit %||% "bp")
    }
  }
  if (length(candidates) != 1L) ggchord_stop(where,
    if (!length(candidates)) "no fully visible matching interval/instance" else "ambiguous panel or instance; specify panel and instance")
  candidates[[1L]]
}

ggchord_fusion_front <- function(panel, endpoint) {
  if (isTRUE(panel$layout$character)) return(ggchord_alignment_front(panel, endpoint))
  layout <- panel$layout; ref <- layout$sequence_reference; id <- endpoint$key
  boundary <- c(min(endpoint$start,endpoint$end)-1,max(endpoint$start,endpoint$end)) - endpoint$origin + 1
  if (endpoint$end < endpoint$start) boundary <- rev(boundary)
  pos <- seq(boundary[1], boundary[2], length.out = 40)
  frame <- ggchord_panel_backbone_frame(layout, id, pos)
  xy <- frame$xy; normal <- frame$normal
  coord <- panel$built$layout$coord
  project <- function(points) {
    # ggplot transforms position scales before Coord; reflecting a circular
    # panel must also reflect its interval anchors, even with symmetric limits.
    scaled <- data.frame(
      x=panel$built$layout$panel_scales_x[[1L]]$transform(points[,1]),
      y=panel$built$layout$panel_scales_y[[1L]]$transform(points[,2]))
    coord$transform(scaled,panel$built$layout$panel_params[[1L]])
  }
  projected <- project(xy)
  if (any(!is.finite(c(projected$x,projected$y))) || any(projected$x < -1e-7 | projected$x > 1+1e-7 | projected$y < -1e-7 | projected$y > 1+1e-7))
    ggchord_stop("Fusion: source front is outside the visible coordinate limits in panel '", endpoint$panel, "'")
  bounds <- panel$context$panel_bounds_mm
  normal_tip <- project(xy+normal)
  normal <- cbind((normal_tip$x-projected$x)*diff(bounds[1:2]),
    (normal_tip$y-projected$y)*diff(bounds[3:4]))
  normal <- normal / sqrt(rowSums(normal^2))
  data.frame(position_boundary = pos+endpoint$origin-1, local_x=xy[,1],local_y=xy[,2],
    x=bounds[1]+projected$x*diff(bounds[1:2]),y=bounds[3]+projected$y*diff(bounds[3:4]),
    normal_x=normal[,1],normal_y=normal[,2],tangent_x=-normal[,2],tangent_y=normal[,1])
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
  ggplot2::set_last_plot(x)
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}

#' @export
#' @importFrom grid grid.draw
grid.draw.ggchord_fusion <- function(x, recording = TRUE) {
  # Keep the composition in the display list, rather than fixed-mm primitives.
  # Recorded plots may be replayed on a device with a different size or DPI.
  # Rebuild there so panels, interval fronts and leaders share the final metric.
  grob <- grid::recordGrob(
    grid::grid.draw(ggchord_build_fusion(plot)$grob, recording = FALSE),
    list(plot = x, ggchord_build_fusion = ggchord_build_fusion)
  )
  grid::grid.draw(grob, recording = recording)
  invisible(x)
}

#' @export
ggplot_build.ggchord_fusion <- function(plot, ...) {
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
      grob <- grid::textGrob(label,gp=grid::gpar(fontsize=text$size %||% 9,
        fontfamily=text$family %||% "",fontface=text$face %||% "plain"))
      if(horizontal) grid::convertWidth(grid::grobWidth(grob),"mm",TRUE) else
        grid::convertHeight(grid::grobHeight(grob),"mm",TRUE)
    },numeric(1))
  required <- max(((head(dimensions,-1)+tail(dimensions,-1))/2+1.2)/(diff(values)/100))
  field <- if(horizontal) "legend.key.width" else "legend.key.height"
  thickness <- if(horizontal) "legend.key.height" else "legend.key.width"
  length <- explicit[[field]] %||% explicit$legend.key.size
  if(is.null(length)) length <- grid::unit(max(required,if(horizontal) 100 else 85),"mm") else
    if(grid::is.unit(length) && grid::convertUnit(length,"mm",valueOnly=TRUE)<required)
      warning("Fusion Identity labels exceed the explicitly fixed guide length; increase its length or reduce its explicit text size",call.=FALSE)
  width <- explicit[[thickness]] %||% explicit$legend.key.size %||% grid::unit(4,"mm")
  guide_theme <- do.call(ggplot2::theme,stats::setNames(list(length,width),c(field,thickness)))
  result <- ggplot2::guide_coloursteps(even.steps=FALSE,show.limits=TRUE,
    direction=if(horizontal) "horizontal" else "vertical",theme=guide_theme)
  attr(result, "ggchord_fusion_default") <- TRUE
  result
}
