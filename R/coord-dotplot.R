#' Two-axis HSP Dotplot coordinates
#'
#' Query positions occupy x and subject positions occupy y. Each axis has
#' independently ordered sequence segments, optionally reversed for display.
#' Coordinates are one-based positions in native sequence units (bp or aa).
#' Concatenation gaps are display spacing, not omitted biological sequence.
#' No sequence backbone layer is required. Use [geom_alignment_segment()] to
#' display HSPs; other ggchord geometries and facets are not yet supported.
#' Standard ggplot2 annotations and explicit position scales remain available.
#'
#' @param x_seqs,y_seqs Optional unique accession vectors in axis order.
#'   Defaults to the query/subject accessions occurring in alignment layers,
#'   or all sequences when there are no alignments. Rows outside a selected
#'   axis set are omitted and recorded; self and repeated HSP rows are retained.
#' @param x_orientation,y_orientation Display direction: 1 or -1, either a
#'   scalar or a vector named by every selected accession. Input is unchanged.
#' @param gap Nonnegative spacing in native sequence units between axis
#'   segments. `NULL` uses five percent of the longest sequence.
#' @param ratio Ratio of physical y units to physical x units, default 1.
#' @param xlim,ylim Optional global concatenated position limits. These zoom
#'   the panel without changing exported HSP endpoints.
#' @param expand Logical. Expand panel limits, default FALSE.
#' @param clip Panel clipping, `"on"` (default) or `"off"`.
#' @return A ggplot2 Coord object.
#' @export
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' ggchord(seq_data_example, ribbon_data_example) +
#'   geom_alignment_segment() +
#'   coord_dotplot()
coord_dotplot <- function(x_seqs = NULL, y_seqs = NULL,
    x_orientation = 1, y_orientation = 1, gap = NULL,
    ratio = 1, xlim = NULL, ylim = NULL, expand = FALSE, clip = "on") {
  if (!is.null(gap) && (!is.numeric(gap) || length(gap) != 1L ||
      !is.finite(gap) || gap < 0)) ggchord_stop("coord_dotplot(): gap must be nonnegative")
  if (!is.numeric(ratio) || length(ratio) != 1L || !is.finite(ratio) || ratio <= 0)
    ggchord_stop("coord_dotplot(): ratio must be positive")
  coord <- ggplot2::coord_fixed(ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = match.arg(clip, c("on", "off")))
  coord$ggchord_dotplot <- TRUE
  coord$dotplot <- list(x_seqs = x_seqs, y_seqs = y_seqs,
    x_orientation = x_orientation, y_orientation = y_orientation, gap = gap)
  coord
}

ggchord_dotplot_axis <- function(seqs, ids, orientation, gap, axis) {
  if (!is.character(ids) || !length(ids) || anyNA(ids) || anyDuplicated(ids) ||
      any(!ids %in% seqs$accver))
    ggchord_stop("coord_dotplot(): ", axis, "_seqs must be unique known accessions")
  if (!is.numeric(orientation) || anyNA(orientation) || any(!orientation %in% c(-1, 1)))
    ggchord_stop("coord_dotplot(): axis orientations must be 1 or -1")
  if (length(orientation) == 1L && is.null(names(orientation))) {
    orientation <- rep(orientation, length(ids))
  } else {
    if (is.null(names(orientation)) || anyDuplicated(names(orientation)) ||
        !setequal(names(orientation), ids))
      ggchord_stop("coord_dotplot(): name each axis orientation by selected accession")
    orientation <- unname(orientation[ids])
  }
  idx <- match(ids, seqs$accver)
  lens <- seqs$length[idx]
  # Position 1 is a point, not the left boundary of a one-base interval.
  # Offsets use length + gap, preserving one source unit between positions.
  offset <- c(0, head(cumsum(lens + gap), -1L))
  data.frame(axis = axis, accver = ids, length = lens,
    offset = offset, orientation = orientation,
    min = offset + 1, max = offset + lens,
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

ggchord_dotplot_geometry <- function(plot) {
  plot <- ggchord_register_layers(plot)
  if (!inherits(plot$facet, "FacetNull")) ggchord_stop("coord_dotplot(): facets are not supported")
  seqs <- plot$ggchord$data$seq_data
  if (anyNA(seqs$accver) || anyDuplicated(seqs$accver) || !is.numeric(seqs$length) ||
      any(!is.finite(seqs$length)) || any(seqs$length < 1))
    ggchord_stop("coord_dotplot(): sequence accessions must be unique with positive finite lengths")
  unit <- ggchord_sequence_unit(seqs, "coord_dotplot()")
  inputs <- list(); query <- subject <- character()
  for (layer in plot$layers) {
    if (is.null(layer$ggchord_type)) next
    if (!identical(layer$ggchord_type, "alignment"))
      ggchord_stop("coord_dotplot(): use geom_alignment_segment(); unsupported layer: ", layer$ggchord_type)
    input <- ggchord_resolve_layer_input(layer, plot$ggchord$data$ribbon_data)
    if (is.null(input)) input <- data.frame(qaccver = character(), saccver = character(),
      qstart = numeric(), qend = numeric(), sstart = numeric(), send = numeric())
    ggchord_require_columns(input, c("qaccver", "saccver", "qstart", "qend", "sstart", "send"), "geom_alignment_segment()")
    ggchord_check_interval_unit(input, unit, "geom_alignment_segment()")
    for (side in c("q", "s")) {
      ids <- as.character(input[[paste0(side, "accver")]])
      if (anyNA(ids) || any(!ids %in% seqs$accver)) ggchord_stop("geom_alignment_segment(): unknown accession")
      for (end in c("start", "end")) {
        pos <- input[[paste0(side, end)]]
        if (!is.numeric(pos) || any(!is.finite(pos)) || any(pos < 1 | pos > seqs$length[match(ids, seqs$accver)]))
          ggchord_stop("geom_alignment_segment(): endpoints must be finite and within sequence length")
      }
    }
    inputs[[layer$ggchord_layer_id]] <- list(alignment = input)
    query <- union(query, as.character(input$qaccver)); subject <- union(subject, as.character(input$saccver))
  }
  config <- plot$coordinates$dotplot
  gap <- config$gap %||% (max(seqs$length) * .05)
  x <- ggchord_dotplot_axis(seqs, config$x_seqs %||% if (length(query)) query else as.character(seqs$accver), config$x_orientation, gap, "x")
  y <- ggchord_dotplot_axis(seqs, config$y_seqs %||% if (length(subject)) subject else as.character(seqs$accver), config$y_orientation, gap, "y")
  permissions <- stats::setNames(lapply(Filter(function(x) identical(x$ggchord_type, "alignment"), plot$layers),
    function(x) isTRUE(x$ggchord_params$allow_repeated)), names(inputs))
  geometries <- list(); reports <- list()
  for (id in names(inputs)) {
    input <- inputs[[id]]$alignment
    if ("source_hsp_id" %in% names(input) &&
        (anyNA(input$source_hsp_id) || any(!nzchar(as.character(input$source_hsp_id)))))
      ggchord_stop("source_hsp_id must be nonempty when supplied")
    reason <- rep("retained", nrow(input))
    reason[!input$qaccver %in% x$accver | !input$saccver %in% y$accver] <- "outside_axis_selection"
    unknown <- ggchord_alignment_unknown(input, "geom_alignment_segment()")
    reason[unknown] <- "unknown_direction"
    if (any(unknown)) warning("geom_alignment_segment(): omitted ", sum(unknown), " HSPs with unknown direction", call. = FALSE)
    rows <- which(reason == "retained")
    d <- input[rows, , drop = FALSE]
    d$x <- ggchord_dotplot_project(x, d$qaccver, d$qstart)
    d$xend <- ggchord_dotplot_project(x, d$qaccver, d$qend)
    d$y <- ggchord_dotplot_project(y, d$saccver, d$sstart)
    d$yend <- ggchord_dotplot_project(y, d$saccver, d$send)
    d$source_row <- rows; d$layer_id <- rep(id, nrow(d)); d$component <- rep("alignment", nrow(d))
    geometries[[id]] <- d
    reports[[id]] <- data.frame(layer_id = rep(id, nrow(input)), source_row = seq_len(nrow(input)),
      retained = reason == "retained", reason = reason,
      source_hsp_id = if ("source_hsp_id" %in% names(input)) as.character(input$source_hsp_id) else rep(NA_character_, nrow(input)))
  }
  report <- if (length(reports)) do.call(rbind, reports) else data.frame()
  if (nrow(report)) {
    declared <- !is.na(report$source_hsp_id) & report$retained
    shown <- report$source_hsp_id[declared]
    repeated <- shown %in% shown[duplicated(shown)]
    if (any(declared) && any(repeated & !unlist(permissions[report$layer_id[declared]])))
      ggchord_stop("A declared source_hsp_id is displayed more than once; use disjoint layer data or allow_repeated=TRUE explicitly")
  }
  axis <- rbind(x, y); axis$sequence_unit <- unit
  list(dotplot = TRUE, axis_layout = axis, alignment = geometries,
    layer_inputs = inputs, connection_report = report,
    sequence_unit = unit, seqs = as.character(seqs$accver),
    xlim = c(1, max(x$max)), ylim = c(1, max(y$max)), rotation = 0, ratio = plot$coordinates$ratio)
}

ggchord_prepare_dotplot <- function(plot) {
  plot <- ggchord_register_layers(plot)
  plot$scales <- plot$scales$clone()
  layout <- ggchord_dotplot_geometry(plot)
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    if (is.null(layer$ggchord_type)) next
    layer <- ggchord_resolve_style_defaults(layer, plot$mapping)
    layer$ggchord_resolved_input <- layout$layer_inputs[[layer$ggchord_layer_id]]$alignment
    plot$layers[[i]] <- reconstruct_layer(layer, layout$alignment[[layer$ggchord_layer_id]])
  }
  plot <- ggchord_add_link_scales(plot)
  for (axis in c("x", "y")) {
    if (plot$scales$has_scale(axis)) next
    frame <- layout$axis_layout[layout$axis_layout$axis == axis, , drop = FALSE]
    breaks <- labels <- vector("list", nrow(frame))
    for (i in seq_len(nrow(frame))) {
      pos <- sort(unique(c(1, pretty(c(1, frame$length[i]), n = if (nrow(frame) > 1L) 2 else 3), frame$length[i])))
      pos <- pos[pos >= 1 & pos <= frame$length[i]]
      breaks[[i]] <- ggchord_dotplot_project(frame, rep(frame$accver[i], length(pos)), pos)
      labels[[i]] <- format(pos, trim = TRUE, scientific = FALSE,
        big.mark = if (nrow(frame) > 1L) "," else "")
    }
    scale <- if (axis == "x") ggplot2::scale_x_continuous else ggplot2::scale_y_continuous
    # Blank expansion geometry trains complete axis extents without censoring
    # additional ordinary annotation layers or changing the source HSP data.
    segment_labels <- paste0(frame$accver, "\n1–",
      format(frame$length, trim = TRUE, scientific = FALSE, big.mark = ","), " ", layout$sequence_unit)
    centers <- (frame$min + frame$max) / 2
    if (nrow(frame) > 1L) {
      # Sequence identity is a separate tier from local numeric positions.
      # The secondary axis is an identity transform, never a new coordinate domain.
      header <- ggplot2::dup_axis(name = NULL, breaks = centers, labels = segment_labels,
        guide = ggplot2::guide_axis(check.overlap = FALSE))
      if (axis == "x") {
        plot$scales$add(scale(name = paste0("Query position (", layout$sequence_unit, ")"),
          breaks = unlist(breaks), labels = unlist(labels), sec.axis = header,
          guide = ggplot2::guide_axis(n.dodge = 2)))
      } else {
        local <- ggplot2::dup_axis(name = paste0("Subject position (", layout$sequence_unit, ")"),
          breaks = unlist(breaks), labels = unlist(labels), guide = ggplot2::guide_axis(n.dodge = 2))
        plot$scales$add(scale(name = "Subject sequence", breaks = centers,
          labels = segment_labels, sec.axis = local,
          guide = ggplot2::guide_axis(check.overlap = FALSE)))
      }
    } else {
      plot$scales$add(scale(name = paste0(if (axis == "x") "Query" else "Subject",
        ": ", frame$accver, " (", layout$sequence_unit, ")"),
        breaks = unlist(breaks), labels = unlist(labels),
        guide = ggplot2::guide_axis(check.overlap = TRUE)))
    }
  }
  plot$layers[[length(plot$layers) + 1L]] <- ggplot2::geom_blank(
    data = data.frame(x = layout$xlim, y = layout$ylim), ggplot2::aes(x, y), inherit.aes = FALSE)
  frame <- layout$axis_layout
  xf <- frame[frame$axis == "x", , drop = FALSE]
  yf <- frame[frame$axis == "y", , drop = FALSE]
  if (nrow(xf) > 1L || nrow(yf) > 1L) {
    cells <- expand.grid(query = seq_len(nrow(xf)), subject = seq_len(nrow(yf)))
    cells$xmin <- xf$min[cells$query]; cells$xmax <- xf$max[cells$query]
    cells$ymin <- yf$min[cells$subject]; cells$ymax <- yf$max[cells$subject]
    # White gaps represent display spacing only. Each rectangle is the actual
    # query/subject extent, so short sequences retain proportionally smaller cells.
    layer <- ggplot2::geom_rect(data = cells,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, show.legend = FALSE,
      fill = ifelse((cells$query + cells$subject) %% 2 == 0, "#F0F6F3", "#FFFFFF"),
      colour = "#A9BBB4", linewidth = .35)
    plot$layers <- c(list(layer), plot$layers)
  }
  base <- ggplot2::theme_gray()
  axes <- base[grepl("^axis\\.", names(base))]
  plot$theme <- plot$theme + do.call(ggplot2::theme, axes) + ggplot2::theme(
    axis.line = ggchord_theme_element(plot, "ggchord.axis.line"),
    axis.text = ggchord_theme_element(plot, "ggchord.axis.text"),
    axis.ticks = ggchord_theme_element(plot, "ggchord.axis.ticks"))
  if (isTRUE(ggchord_plot_settings(plot)$axis$hidden))
    plot$theme <- plot$theme + ggplot2::theme(axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  if (!is.null(plot$ggchord$native_theme)) plot$theme <- plot$theme + plot$ggchord$native_theme
  plot$ggchord$ref$layout <- layout
  plot
}

ggchord_export_dotplot <- function(plot, include, original_data) {
  layout <- get_chord_layout(plot)
  out <- stats::setNames(lapply(include, function(x) data.frame()), include)
  if ("alignment" %in% include && length(layout$alignment)) {
    # Explicit layers may have different metadata columns.
    out$alignment <- ggchord_rbind_fill(layout$alignment)
  }
  if ("axis" %in% include) out$axis <- layout$axis_layout
  out$axis_layout <- layout$axis_layout
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
  axis <- match.arg(axis, c("x", "y"))
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
