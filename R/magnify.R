#' Mark an expanded source interval in a Fusion
#' @param from,to Character vector identifying the source and detail panel paths.
#' @param region Source intervals with accver, start and end, and optional
#'   from_instance and to_instance selectors.
#' @param highlight Whether to draw a highlighter on the source interval.
#' @param fill,colour,alpha,linewidth Fixed appearance of the decoration.
#' @param id Optional root layer ID.
#' @param connector Join the selected intervals, or lead to the facing
#'   `"border"` of the detail panel.
#' @param linetype Connector line type. By default, border leaders are dashed
#'   and interval outlines are solid.
#' @return A decorative component for [ggchord_fusion()].
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
#' f <- ggchord_fusion(list(overview = p_basic, detail = p_basic))
#' region <- data.frame(accver = "A", start = 100, end = 240)
#' f + fusion_zoom("overview", "detail", region, id = "zoom", connector = "border")
fusion_zoom <- function(from, to, region, highlight = TRUE,
    fill = "#187B77", alpha = 0.12, colour = "#187B77", linewidth = 0.4, id = NULL,
    connector = "interval", linetype = NULL) {
  caller <- "fusion_zoom()"
  connector <- ggchord_match_choice(connector, c("interval", "border"), "connector")
  if (!is.null(linetype)) {
    valid_codes <- c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    if (length(linetype) != 1L || anyNA(linetype) ||
        !(is.numeric(linetype) && linetype %in% 0:6 ||
          is.character(linetype) && (linetype %in% valid_codes ||
            grepl("^([1-9A-Fa-f]{2}){1,4}$", linetype))))
      ggchord_stop(caller, ": use a single grid line type or NULL")
  }
  from_key <- ggchord_fusion_check_path(from, caller)
  to_key <- ggchord_fusion_check_path(to, caller)
  if (identical(from_key, to_key)) ggchord_stop(caller, ": provide distinct panels")
  ggchord_require_columns(region, c("accver", "start", "end"), caller)
  if (!is.character(region$accver) || anyNA(region$accver) || any(!nzchar(region$accver)))
    ggchord_stop(caller, ": accver must contain nonempty source IDs")
  ggchord_project_integer(region$start, "start", caller)
  ggchord_project_integer(region$end, "end", caller)
  if (any(region$start > region$end)) ggchord_stop(caller, ": intervals must have ordered endpoints")
  for (name in intersect(c("from_instance", "to_instance"), names(region))) {
    if (!is.character(region[[name]]) || any(!is.na(region[[name]]) & !nzchar(region[[name]])))
      ggchord_stop(caller, ": instance selectors must be strings")
  }
  if (!is.logical(highlight) || length(highlight) != 1L || is.na(highlight))
    ggchord_stop(caller, ": highlight must be TRUE or FALSE")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1 ||
      !is.numeric(linewidth) || length(linewidth) != 1L || !is.finite(linewidth) || linewidth < 0)
    ggchord_stop(caller, ": alpha must be in [0,1] and linewidth must be nonnegative")
  for (value in list(fill, colour)) {
    if (length(value) != 1L) ggchord_stop(caller, ": colours must be scalar")
    tryCatch(col2rgb(value), error = function(e) ggchord_stop(caller, ": invalid colour"))
  }
  structure(list(from = from, to = to, from_key = from_key, to_key = to_key,
    region = as.data.frame(region), highlight = highlight, fill = fill,
    alpha = alpha, colour = colour, linewidth = linewidth, id = id,
    connector = connector, linetype = linetype),
    class = "ggchord_fusion_zoom")
}

#' Compose an overview and a detail with a source interval marker
#' @param overview,detail Independent plots.
#' @param region Source intervals passed to [fusion_zoom()].
#' @param layout An optional layout for the overview and detail.
#' @param guides Policy passed to [ggchord_fusion()].
#' @param ... Fixed zoom styles or a root decoration ID.
#' @return A Fusion with the two plots and a zoom decoration.
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
#' windows <- data.frame(accver = "A", start = 1, end = 500, window_id = "w1")
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' region <- data.frame(accver = "A", start = 100, end = 240)
#' p_detail <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) +
#'   layout_sequences(windows = data.frame(instance_id = "A_zoom", accver = "A",
#'     start = 80, end = 260, row_id = "r1")) + coord_collinear()
#' magnify_ggchord(p_basic, p_detail, region)
magnify_ggchord <- function(overview, detail, region, layout = NULL,
    guides = "keep", ...) {
  if (is.null(layout))
    layout <- layout_fusion_grid(matrix(c("overview", "detail"), nrow = 1L),
      widths = c(2, 1))
  args <- list(...)
  allowed <- c("highlight", "fill", "alpha", "colour", "linewidth", "id",
    "connector", "linetype")
  if (length(args) && (length(names(args)) != length(args) || anyDuplicated(names(args)) ||
      !all(names(args) %in% allowed)))
    ggchord_stop("magnify_ggchord(): supply only named zoom styles or id in ...")
  ggchord_fusion(list(overview = overview, detail = detail), layout = layout, guides = guides) +
    fusion_zoom("overview", "detail", region, ...)
}

ggchord_render_zoom <- function(layer, panels) {
  connector <- layer$connector %||% "interval"
  line <- layer$linetype %||% if (connector == "border") "dashed" else "solid"
  outline <- gpar(fill = alpha(layer$fill, layer$alpha),
    col = layer$colour, lwd = layer$linewidth * ggplot2::.pt, lty = line)
  resolve <- function(region, panel, name) {
    instance <- region[[name]]
    if (length(instance) && is.na(instance)) instance <- NULL
    ggchord_resolve_panel_interval(panels, region$accver, region$start, region$end,
      panel, instance, paste0("Fusion zoom '", layer$id, "': "))
  }
  source_overlay <- function(front) {
    half <- front$cell_half_height_mm %||% rep(0.7, nrow(front))
    x <- front$x + front$normal_x * half
    y <- front$y + front$normal_y * half
    xx <- front$x - front$normal_x * half
    yy <- front$y - front$normal_y * half
    polygonGrob(c(x, rev(xx)), c(y, rev(yy)), default.units = "mm",
      gp = gpar(fill = alpha(layer$fill, layer$alpha),
        col = layer$colour, lwd = layer$linewidth * ggplot2::.pt))
  }
  parts <- lapply(seq_len(nrow(layer$region)), function(i) {
    region <- layer$region[i, , drop = FALSE]
    q <- resolve(region, layer$from_key, "from_instance")
    s <- resolve(region, layer$to_key, "to_instance")
    ggchord_check_character_origin(q, s, "fusion_zoom()")
    if (q$unit != s$unit) ggchord_stop("fusion_zoom(): source and detail must use the same unit")
    qf <- ggchord_fusion_front(panels[[q$panel]], q)
    sf <- ggchord_fusion_front(panels[[s$panel]], s)
    region$from_instance <- q$instance
    region$to_instance <- s$instance
    region$from_panel_path <- list(layer$from)
    region$to_panel_path <- list(layer$to)
    region$sequence_unit <- q$unit
    region$layer_id <- layer$id
    region$source_row <- i
    region$connector <- connector
    region$component_type <- "decoration"
    region$composition_id <- ""
    path <- grob <- NULL
    if (connector == "interval") {
      xy <- ggchord_fusion_path(qf, sf, curvature = 0)
      path <- data.frame(x = xy[, 1L], y = xy[, 2L], layer_id = layer$id,
        source_row = i, component_type = "decoration", composition_id = "", units = "mm")
      grob <- polygonGrob(path$x, path$y, default.units = "mm", gp = outline)
    }
    list(region = region, path = path, grob = grob, front = qf,
      source_facet = q$facet %||% 1L, target_facet = s$facet %||% 1L,
      overlay = if (layer$highlight) setNames(list(source_overlay(qf)), q$panel) else list())
  })
  if (connector == "interval" || !length(parts)) return(parts)

  # Draw one envelope for all source intervals feeding the same two facets.
  # Preserve each source highlight and interval in the exported region table.
  groups <- split(seq_along(parts), vapply(parts, function(part)
    paste(part$source_facet, part$target_facet, sep = ":"), character(1)))
  bounds <- function(panel, facet) {
    context <- panel$context
    result <- context$panel_bounds_mm
    leaf <- if (length(context$facet_contexts) >= facet) context$facet_contexts[[facet]]
    if (!is.null(leaf))
      result <- leaf$panel_bounds_mm - (context$device_panel_bounds_mm - result)
    matrix(result, nrow = 2L)
  }
  for (where in groups) {
    first <- parts[[where[1L]]]
    from <- bounds(panels[[layer$from_key]], first$source_facet)
    to <- bounds(panels[[layer$to_key]], first$target_facet)
    direction <- colMeans(from) - colMeans(to)
    separating <- pmax(from[1L, ] - to[2L, ], to[1L, ] - from[2L, ], 0)
    normal <- if (any(separating > 0)) which.max(separating) else which.max(abs(direction))
    side <- if (direction[normal] > 0) 2L else 1L
    along <- 3L - normal
    tips <- do.call(rbind, lapply(where, function(i) {
      front <- parts[[i]]$front
      data.frame(x = front$x, y = front$y, source_row = i)
    }))
    pick <- c(which.min(tips[[along]]), which.max(tips[[along]]))
    if (pick[1L] == pick[2L]) pick <- c(which.min(tips[[normal]]), which.max(tips[[normal]]))
    roots <- tips[pick, , drop = FALSE]
    edge <- matrix(rep(to[side, ], each = 2L), nrow = 2L)
    edge[, along] <- to[, along]
    xy <- rbind(as.matrix(roots[c("x", "y")]), edge)
    path <- data.frame(x = xy[, 1L], y = xy[, 2L], layer_id = layer$id,
      source_row = rep(roots$source_row, 2L), path_part = rep(1:2, 2L),
      component_type = "decoration", composition_id = "", units = "mm")
    parts[[where[1L]]]$path <- path
    parts[[where[1L]]]$grob <- segmentsGrob(roots$x, roots$y, edge[, 1L], edge[, 2L],
      default.units = "mm", gp = outline)
    for (i in where[-1L]) parts[[i]]$grob <- nullGrob()
  }
  parts
}
