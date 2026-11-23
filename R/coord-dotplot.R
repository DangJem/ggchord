#' Project HSPs onto two source-reference axes
#'
#' @param ratio Positive ratio of y to x units on the device, or NULL.
#' @param xlim,ylim Optional increasing pairs of view limits.
#' @param expand Whether to expand automatic view limits.
#' @param clip Panel clipping, `"on"` or `"off"`.
#' @param fit Fit geometry, include labels, or require manual limits.
#' @return A coordinate system for a ggchord plot.
#' @export
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' ggchord(data = seqs) +
#'   geom_alignment_segment(data = links, aes(colour = pident)) +
#'   layout_dotplot(x_order = "A", y_order = "B") +
#'   coord_dotplot(ratio = 1)
coord_dotplot <- function(ratio = 1, xlim = NULL, ylim = NULL,
    expand = FALSE, clip = "on", fit = "geometry") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  settings <- ggchord_coord_settings(ratio, xlim, ylim, expand, clip, fit)
  ggproto("CoordDotplot", ggplot2::CoordCartesian, !!!settings, ggchord_dotplot = TRUE)
}

ggchord_dotplot_axis <- function(seqs, ids, orientation, gap, axis) {
  if (!is.character(ids) || !length(ids) || anyNA(ids) || anyDuplicated(ids) ||
      any(!ids %in% seqs$accver))
    ggchord_stop("coord_dotplot(): ", axis, "_seqs must be unique known accessions")
  orientation <- ggchord_layout_values(orientation, ids, field = paste0(axis, "_orientation"), default = 1)
  gap <- ggchord_layout_values(gap, ids, field = paste0(axis, "_gap"), default = 0)
  if (any(!orientation %in% c(-1, 1))) ggchord_stop("Dotplot orientation must be 1 or -1")
  idx <- match(ids, seqs$accver)
  lens <- seqs$length[idx]
  unit <- ggchord_sequence_unit(seqs[idx, , drop = FALSE], paste("Dotplot", axis, "axis"))
  # Position 1 is a point, not the left boundary of a one-base interval.
  # Offsets use length + gap, preserving one source unit between positions.
  offset <- c(0, head(cumsum(lens + gap), -1L))
  data.frame(axis = axis, accver = ids, length = lens,
    offset = offset, orientation = orientation,
    min = offset + 1, max = offset + lens, sequence_unit = unit,
    source_accver = if (".source_accver" %in% names(seqs)) seqs$.source_accver[idx] else ids,
    source_origin = if (".source_origin" %in% names(seqs)) seqs$.source_origin[idx] else 1,
    stringsAsFactors = FALSE)
}

ggchord_dotplot_project <- function(frame, accver, position, inverse = FALSE) {
  idx <- match(accver, frame$accver)
  if (anyNA(idx)) ggchord_stop("Dotplot projection: accession is not on this axis")
  offset <- frame$offset[idx]; len <- frame$length[idx]; dir <- frame$orientation[idx]
  if (inverse) {
    local <- position - offset
    result <- ifelse(dir == 1, local, len + 1 - local)
    result[local < 1 | local > len] <- NA_real_
  } else {
    if (any(!is.finite(position)) || any(position < 1 | position > len))
      ggchord_stop("Dotplot projection: positions must be within sequence length")
    result <- offset + ifelse(dir == 1, position, len + 1 - position)
  }
  unname(result)
}

ggchord_dotplot_geometry <- function(plot, for_drawing = FALSE) {
  ggchord_validate_coord_contract(plot)
  keys <- ggchord_facet_vars(plot)
  seqs <- plot$ggchord$data$seq_data
  if (anyNA(seqs$accver) || anyDuplicated(seqs$accver) || !is.numeric(seqs$length) ||
      any(!is.finite(seqs$length)) || any(seqs$length < 1))
    ggchord_stop("coord_dotplot(): sequence accessions must be unique with positive finite lengths")
  units <- setNames(seqs$sequence_unit, seqs$accver)
  inputs <- list()
  for (layer in plot$layers) {
    if (is.null(layer$ggchord_type)) next
    if (!identical(layer$ggchord_type, "alignment"))
      ggchord_stop("coord_dotplot(): use geom_alignment_segment(); unsupported layer: ", layer$ggchord_type)
    fallback <- ggchord_layer_fallback(layer, plot$ggchord$data)
    source <- ggchord_layer_source(layer, fallback)
    input <- ggchord_resolve_layer_input(layer, fallback, data = source)
    if (is.null(input)) input <- data.frame(qaccver = character(), saccver = character(),
      qstart = numeric(), qend = numeric(), sstart = numeric(), send = numeric())
    ggchord_require_columns(input, c("qaccver", "saccver", "qstart", "qend", "sstart", "send"), "geom_alignment_segment()")
    if (length(keys)) {
      ggchord_check_panel_layout(layer)
      for (key in intersect(keys, names(source))) if (!identical(as.character(input[[key]]), as.character(source[[key]])))
        ggchord_stop(ggchord_layer_input_context(layer), ": prepare canonical ", key, " values before faceting")
    }
    ggchord_check_reference_ids(input, c("qaccver", "saccver"), seqs$accver,
      ggchord_layer_input_context(layer))
    for (side in c("q", "s")) {
      ids <- as.character(input[[paste0(side, "accver")]])
      unit <- input[[paste0(side, "sequence_unit")]] %||% input$sequence_unit
      if (!is.null(unit) && any(!is.na(unit) & unit != units[ids]))
        ggchord_stop("geom_alignment_segment(): endpoint sequence unit disagrees with the reference")
      for (end in c("start", "end")) {
        pos <- input[[paste0(side, end)]]
        if (!is.numeric(pos) || any(!is.finite(pos)) || any(pos < 1 | pos > seqs$length[match(ids, seqs$accver)]))
          ggchord_stop("geom_alignment_segment(): endpoints must be finite and within sequence length")
      }
    }
    inputs[[layer$ggchord_layer_id]] <- list(alignment = input)
    if (for_drawing) inputs[[layer$ggchord_layer_id]]$source <- source
  }
  config <- plot$ggchord$dotplot_layout %||% layout_dotplot()
  facet <- if (length(keys)) {
    frames <- lapply(inputs, `[[`, "alignment")
    native <- Filter(function(l) is.null(l$ggchord_type), plot$layers)
    frames <- c(frames, lapply(native, function(l) l$layer_data(plot$data)))
    ggchord_facet_panels(plot, frames)
  }
  panels <- facet$panels
  free <- plot$facet$params$free
  count <- if (is.null(panels)) 1L else nrow(panels)
  membership <- if (!is.null(panels)) ggchord_facet_backbones(plot, panels, facet$params)
  build_axis <- function(axis) {
    ids <- config[[paste0(axis, "_order")]] %||% seqs$accver
    frame <- ggchord_dotplot_axis(seqs, ids, config[[paste0(axis, "_orientation")]],
      config[[paste0(axis, "_gap")]], axis)
    if (!isTRUE(free[[axis]])) return(rep(list(frame), count))
    key <- paste0("SCALE_", toupper(axis))
    scales <- panels[[key]]
    frames <- lapply(unique(scales), function(group) {
      selected <- unique(unlist(lapply(membership[scales == group], `[[`, "displayed")))
      local <- frame[frame$accver %in% selected, , drop = FALSE]
      if (!nrow(local)) ggchord_stop("No complete reference for free Dotplot ", axis, " axis")
      ggchord_dotplot_axis(seqs, local$accver, setNames(local$orientation, local$accver),
        ggchord_layout_values(config[[paste0(axis, "_gap")]], ids, field = "gap", default = 0)[local$accver], axis)
    })
    frames[match(scales, unique(scales))]
  }
  x <- build_axis("x")
  y <- build_axis("y")
  alignment_layers <- setNames(Filter(function(x) identical(x$ggchord_type, "alignment"),
    plot$layers), names(inputs))
  geometries <- list(); reports <- list()
  for (id in names(inputs)) {
    input <- inputs[[id]]$alignment
    if (any(input$is_summary %in% TRUE) || any(input$aggregated %in% TRUE))
      ggchord_stop("geom_alignment_segment(): aggregated intervals are not HSPs")
    direction <- ggchord_reference_direction(input, "geom_alignment_segment()")
    if (anyNA(direction)) ggchord_stop("geom_alignment_segment(): missing HSP direction",
      ggchord_input_rows(is.na(direction)))
    if (is.null(panels)) {
      members <- data.frame(.ggchord_row = seq_len(nrow(input)), PANEL = rep(1L, nrow(input)))
    } else {
      rows <- input[intersect(keys, names(input))]
      rows$.ggchord_row <- seq_len(nrow(input))
      attr(rows, "layout") <- alignment_layers[[id]]$layout
      members <- plot$facet$map_data(rows, panels, facet$params)
    }
    drawn <- vector("list", count); skipped <- vector("list", count)
    for (i in seq_len(count)) {
      rows <- members$.ggchord_row[as.integer(members$PANEL) == i]
      d <- input[rows, , drop = FALSE]
      keep <- d$qaccver %in% x[[i]]$accver & d$saccver %in% y[[i]]$accver
      skipped[[i]] <- data.frame(layer_id = rep(id, length(rows)), source_row = rows,
        PANEL = rep(i, length(rows)), retained = keep,
        reason = ifelse(keep, "retained", "outside_axis_selection"))
      d <- d[keep, , drop = FALSE]
      slope <- direction[rows[keep]]
      d$x <- ggchord_dotplot_project(x[[i]], d$qaccver, d$qstart)
      d$xend <- ggchord_dotplot_project(x[[i]], d$qaccver, d$qend)
      d$y <- ggchord_dotplot_project(y[[i]], d$saccver, ifelse(slope > 0, d$sstart, d$send))
      d$yend <- ggchord_dotplot_project(y[[i]], d$saccver, ifelse(slope > 0, d$send, d$sstart))
      d$source_row <- rows[keep]; d$layer_id <- rep(id, nrow(d))
      d$component <- rep("alignment", nrow(d))
      d$PANEL <- rep(i, nrow(d))
      drawn[[i]] <- d
    }
    geometries[[id]] <- do.call(rbind, drawn)
    reports[[id]] <- do.call(rbind, skipped)
  }
  report <- ggchord_rbind_fill(reports) %||% data.frame()
  axis_panels <- lapply(seq_len(count), function(i) {
    out <- rbind(x[[i]], y[[i]])
    out$PANEL <- i
    out
  })
  axis <- if (any(unlist(free))) do.call(rbind, axis_panels) else {
    out <- axis_panels[[1L]]
    out$PANEL <- NULL
    out
  }
  list(dotplot = TRUE, axis_layout = axis, alignment = geometries,
    layer_inputs = inputs, connection_report = report, facet_layout = facet$panels,
    sequence_unit = if (length(unique(axis$sequence_unit)) == 1L) axis$sequence_unit[1L] else NA_character_,
    seqs = as.character(seqs$accver),
    axis_panels = axis_panels, free_axes = free,
    xlim = c(1, max(vapply(x, function(d) max(d$max), numeric(1)))),
    ylim = c(1, max(vapply(y, function(d) max(d$max), numeric(1)))),
    rotation = 0, ratio = plot$coordinates$ratio)
}

ggchord_prepare_dotplot <- function(plot) {
  plot <- ggchord_prepare_position_scales(ggchord_register_layers(plot))
  plot$scales <- plot$scales$clone()
  layout <- ggchord_dotplot_geometry(plot, for_drawing = TRUE)
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    if (is.null(layer$ggchord_type)) next
    layer <- ggchord_resolve_style_defaults(layer, plot$mapping)
    layer <- ggchord_bind_visual_mapping(layer, plot)
    inputs <- layout$layer_inputs[[layer$ggchord_layer_id]]
    layout$layer_inputs[[layer$ggchord_layer_id]]$source <- NULL
    layer$ggchord_resolved_input <- inputs$alignment
    geometry <- layout$alignment[[layer$ggchord_layer_id]]
    drawing <- inputs$alignment[geometry$source_row, , drop = FALSE]
    roles <- intersect(layer$ggchord_role_aes, names(inputs$source))
    drawing[roles] <- inputs$source[geometry$source_row, roles, drop = FALSE]
    mapping <- ggchord_effective_mapping(layer)
    # Keep user columns (including x, y or source_row) visible to visual
    # mappings. Only the geom's four position aesthetics use projected data.
    for (axis in c("x", "y", "xend", "yend")) {
      field <- tail(make.unique(c(names(drawing), paste0(".ggchord_", axis))), 1L)
      drawing[[field]] <- geometry[[axis]]
      mapping[[axis]] <- new_quosure(sym(field), emptyenv())
    }
    if (!is.null(layout$facet_layout)) drawing$.ggchord_panel <- geometry$PANEL
    plot$layers[[i]] <- reconstruct_layer(layer, drawing, mapping)
    if (!is.null(layout$facet_layout)) plot$layers[[i]]$layout <- NULL
  }
  plot <- ggchord_add_link_scales(plot)
  for (axis in c("x", "y")) {
    if (plot$scales$has_scale(axis)) next
    axis_data <- layout$axis_panels[[1L]]
    frame <- axis_data[axis_data$axis == axis, , drop = FALSE]
    unit <- frame$sequence_unit[1L]
    if (isTRUE(layout$free_axes[[axis]])) {
      scale <- if (axis == "x") scale_x_continuous else scale_y_continuous
      name <- paste0(if (axis == "x") "Query position" else "Subject position", " (", unit, ")")
      plot$scales$add(scale(name = name))
      next
    }
    breaks <- labels <- vector("list", nrow(frame))
    for (i in seq_len(nrow(frame))) {
      extent <- c(1, frame$length[i])
      pos <- pretty(extent, n = if (nrow(frame) > 1L) 2 else 3)
      pos <- sort(unique(c(extent, pos[pos >= extent[1L] & pos <= extent[2L]])))
      breaks[[i]] <- ggchord_dotplot_project(frame, rep(frame$accver[i], length(pos)), pos)
      labels[[i]] <- format(pos, trim = TRUE, scientific = FALSE,
        big.mark = if (nrow(frame) > 1L) "," else "")
    }
    scale <- if (axis == "x") scale_x_continuous else scale_y_continuous
    # Blank expansion geometry trains complete axis extents without censoring
    # additional ordinary annotation layers or changing the source HSP data.
    segment_labels <- paste0(frame$accver, "\n1-",
      format(frame$length, trim = TRUE, scientific = FALSE, big.mark = ","), " ", unit)
    centers <- (frame$min + frame$max) / 2
    if (nrow(frame) > 1L) {
      # Sequence identity is a separate tier from local numeric positions.
      # The secondary axis is an identity transform, never a new coordinate domain.
      header <- ggplot2::dup_axis(name = NULL, breaks = centers, labels = segment_labels,
        guide = guide_axis(check.overlap = FALSE))
      if (axis == "x") {
        plot$scales$add(scale(name = paste0("Query position (", unit, ")"),
          breaks = unlist(breaks), labels = unlist(labels), sec.axis = header,
          guide = guide_axis(check.overlap = TRUE, n.dodge = 2)))
      } else {
        local <- ggplot2::dup_axis(name = paste0("Subject position (", unit, ")"),
          breaks = unlist(breaks), labels = unlist(labels), guide = guide_axis(check.overlap = TRUE))
        plot$scales$add(scale(name = "Subject sequence", breaks = centers,
          labels = segment_labels, sec.axis = local,
          guide = guide_axis(check.overlap = FALSE)))
      }
    } else {
      plot$scales$add(scale(name = paste0(if (axis == "x") "Query" else "Subject",
        ": ", frame$accver, " (", unit, ")"),
        breaks = unlist(breaks), labels = unlist(labels),
        guide = guide_axis(check.overlap = TRUE)))
    }
  }
  blank <- do.call(rbind, lapply(layout$axis_panels, function(frame) {
    data.frame(x = c(1, max(frame$max[frame$axis == "x"])),
      y = c(1, max(frame$max[frame$axis == "y"])), .ggchord_panel = frame$PANEL[1L])
  }))
  plot$layers[[length(plot$layers) + 1L]] <- ggplot2::geom_blank(
    data = blank, aes(x, y), inherit.aes = FALSE)
  if (!is.null(layout$facet_layout))
    plot$facet <- ggchord_precomputed_facet(plot$facet, layout$facet_layout)
  frame <- layout$axis_panels[[1L]]
  xf <- frame[frame$axis == "x", , drop = FALSE]
  yf <- frame[frame$axis == "y", , drop = FALSE]
  if (!any(unlist(layout$free_axes)) && (nrow(xf) > 1L || nrow(yf) > 1L)) {
    cells <- expand.grid(query = seq_len(nrow(xf)), subject = seq_len(nrow(yf)))
    cells$xmin <- xf$min[cells$query]; cells$xmax <- xf$max[cells$query]
    cells$ymin <- yf$min[cells$subject]; cells$ymax <- yf$max[cells$subject]
    # White gaps represent display spacing only. Each rectangle is the actual
    # query/subject extent, so short sequences retain proportionally smaller cells.
    cells$background <- ifelse((cells$query + cells$subject) %% 2 == 0, "#F0F6F3", "#FFFFFF")
    layer <- ggplot2::geom_rect(data = cells,
      aes(xmin = !!sym("xmin"), xmax = !!sym("xmax"),
          ymin = !!sym("ymin"), ymax = !!sym("ymax"),
          fill = I(!!sym("background"))),
      inherit.aes = FALSE, show.legend = FALSE,
      colour = "#A9BBB4", linewidth = .35)
    plot$layers <- c(list(layer), plot$layers)
  }
  base <- ggplot2::theme_gray()
  axes <- base[grepl("^axis\\.", names(base))]
  plot$theme <- plot$theme + do.call(theme, axes) + theme(
    axis.line = ggchord_theme_element(plot, "ggchord.axis.line"),
    axis.text = ggchord_theme_element(plot, "ggchord.axis.text"),
    axis.ticks = ggchord_theme_element(plot, "ggchord.axis.ticks"))
  if (isTRUE(ggchord_plot_settings(plot)$axis$hidden))
    plot$theme <- plot$theme + theme(axis.line = element_blank(),
      axis.text = element_blank(), axis.ticks = element_blank())
  if (!is.null(plot$ggchord$native_theme)) plot$theme <- plot$theme + plot$ggchord$native_theme
  ggchord_cache_layout(plot, layout)
  plot
}

ggchord_export_dotplot <- function(plot, include, original_data) {
  layout <- get_chord_layout(plot)
  out <- setNames(lapply(include, function(x) data.frame()), include)
  if ("alignment" %in% include && length(layout$alignment)) {
    # Explicit layers may have different metadata columns.
    out$alignment <- ggchord_rbind_fill(layout$alignment)
  }
  if ("axis" %in% include) out$axis <- layout$axis_layout
  out$axis_layout <- layout$axis_layout
  if (!is.null(layout$facet_layout)) out$facet_layout <- layout$facet_layout
  out$connection_report <- layout$connection_report
  out$metadata <- list(coordinate_system = "dotplot", sequence_unit = layout$sequence_unit,
    position_convention = "one-based positions", xlim = plot$coordinates$limits$x %||% layout$xlim,
    ylim = plot$coordinates$limits$y %||% layout$ylim)
  if (original_data) out$original_data <- layout$layer_inputs
  out$source_registry <- ggchord_source_registry(plot$ggchord$sources %||% list())
  if (!is.null(plot$ggchord$layout_policy)) out$adaptive_layout <- list(
    policy=plot$ggchord$layout_policy,measurement="layout_only",guide_layout="target_device_required")
  class(out) <- c("ggchord_layout_export", "list")
  out
}

ggchord_map_dotplot <- function(plot, accver, position, inverse, source, axis) {
  if (is.null(axis)) ggchord_stop("map_ggchord_position(): specify axis = 'x' or 'y' for a Dotplot")
  axis <- ggchord_match_arg(axis, c("x", "y"))
  if (!is.logical(inverse) || length(inverse) != 1L || is.na(inverse) ||
      !is.logical(source) || length(source) != 1L || is.na(source))
    ggchord_stop("map_ggchord_position(): inverse/source must be TRUE or FALSE")
  if (!is.numeric(position) || any(!is.finite(position)) || anyNA(accver) ||
      !length(accver) %in% c(1L, length(position)))
    ggchord_stop("map_ggchord_position(): supply finite positions and matching accessions")
  layout <- get_chord_layout(plot)
  frame <- layout$axis_layout[layout$axis_layout$axis == axis, , drop = FALSE]
  ids <- rep_len(as.character(accver), length(position))
  origin <- frame$source_origin[match(ids, frame$accver)]
  pos <- position - if (source && !inverse) origin - 1 else 0
  value <- ggchord_dotplot_project(frame, ids, pos, inverse)
  local <- if (inverse) value else pos
  global <- if (inverse) position else value
  data.frame(accver = ids, position = local + if (source) origin - 1 else 0,
    axis = axis, x = if (axis == "x") global else NA_real_,
    y = if (axis == "y") global else NA_real_, sequence_unit = layout$sequence_unit)
}
