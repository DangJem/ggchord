#' Query source, local and display positions
#' @param plot A ggchord plot, registered character plot or Fusion.
#' @param data Query rows with source_accver/source_position, accver/pos, or
#'   x/y in normalised panel coordinates, according to from.
#' @param from,to Coordinate spaces: source, local or display.
#' @param axis A required x or y selector in a Dotplot.
#' @param panel The target panel ID path in a Fusion.
#' @param instance Optional exact display instance ID.
#' @param facet Named list of actual facet values. For display queries,
#'   matching panels are returned with their facet keys.
#' @inheritParams export_ggchord_layout
#' @return A correspondence table retaining all matching instances and recording
#'   reasons for unavailable positions.
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
#' map_ggchord_position(p_basic,
#'   data.frame(source_accver = "A", source_position = 150),
#'   width = 8, height = 5, device = "pdf")
map_ggchord_position <- function(plot, data, from = "source", to = "display",
    axis = NULL, panel = NULL, instance = NULL, width = NULL, height = NULL,
    units = "in", device = NULL, dpi = 300, facet = NULL) {
  ggchord_local_console_errors()
  caller <- "map_ggchord_position()"
  ggchord_output_plot(plot, caller)
  spaces <- c("source", "local", "display")
  from <- ggchord_match_arg(from, spaces)
  to <- ggchord_match_arg(to, spaces)
  fields <- list(source = c("source_accver", "source_position"),
    local = c("accver", "pos"), display = c("x", "y"))
  ggchord_require_columns(data, fields[[from]], caller)
  request <- ggchord_output_request(width, height, units, device, dpi, caller)
  needs_display <- "display" %in% c(from, to)
  path <- character()
  facet_selector <- ggchord_fusion_facet(facet, caller)
  facet_number <- NULL
  leaf <- plot
  if (inherits(plot, "ggchord_fusion")) {
    ggchord_fusion_check_path(panel, caller)
    path <- unname(panel)
    leaf <- ggchord_fusion_leaf(plot, path, caller)
  } else if (!is.null(panel)) {
    if (!is.numeric(panel) || length(panel) != 1L || !is.finite(panel) ||
        panel < 1 || panel != trunc(panel)) ggchord_stop(caller, ": invalid panel number")
    facet_number <- panel
  }
  ggchord_output_plot(leaf, caller)
  if (inherits(leaf, "ggchord_fusion")) ggchord_stop(caller, ": select a leaf panel")
  no_instance <- is.null(instance)
  if (!no_instance) ggchord_check_layer_id(instance)
  dotplot <- isTRUE(leaf$coordinates$ggchord_dotplot)
  if (dotplot && is.null(axis)) ggchord_stop(caller, ": supply a Dotplot axis")
  if (!dotplot && !is.null(axis)) ggchord_stop(caller, ": axis applies only to Dotplot")
  if (dotplot) axis <- ggchord_match_arg(axis, c("x", "y"))
  reference <- ggchord_query_reference(leaf, axis)
  if (!no_instance && !instance %in% reference$instance_id)
    ggchord_stop(caller, ": unknown instance '", instance, "'")
  if (!no_instance) reference <- reference[reference$instance_id == instance, , drop = FALSE]
  input_pos <- if (from == "display") if (dotplot) axis else c("x", "y") else fields[[from]][2L]
  if (!all(vapply(data[input_pos], function(x)
      is.numeric(x) && all(is.finite(x) | is.na(x)), logical(1))))
    ggchord_stop(caller, ": query positions must be finite numbers or NA")
  failure <- NULL
  context <- NULL
  if (needs_display && !is.null(request$size)) {
    context <- tryCatch(ggchord_with_output_device(request, function() {
      on_error <- if (length(path)) function(panel_path, error, layer_id = NULL) {
        if (is.null(layer_id) && identical(head(path, length(panel_path)), panel_path))
          failure <<- conditionMessage(error)
      }
      result <- ggchord_output_build_geometry(ggchord_query_clone(plot), on_error)
      if (!length(path)) result$leaf else result$leaves[[ggchord_fusion_path_key(path)]]
    }), error = function(e) {
      failure <<- conditionMessage(e)
      NULL
    })
    if (!is.null(context) && !identical(context$context$status, "measured")) context <- NULL
  }
  if (isTRUE(context$layout$character)) {
    reference <- ggchord_query_character_reference(context$layout$character_data)
    if (!no_instance) reference <- reference[reference$instance_id == instance, , drop = FALSE]
  }
  measurement <- ggchord_output_measurement(request, !is.null(context), failure,
    actual = attr(context, "measurement_device"))
  empty <- data.frame(input_row = integer(), source_accver = character(),
    source_position = numeric(), accver = character(), pos = numeric(),
    x = numeric(), y = numeric(), instance_id = character(), PANEL = integer(),
    sequence_unit = character(), status = character(), reason = character(),
    panel_path = I(list()), facet_key = I(list()))
  if (dotplot) empty$axis <- character()
  result <- lapply(seq_len(nrow(data)), function(i) {
    input <- data[i, , drop = FALSE]
    miss <- function(reason, status = "unavailable") {
      row <- empty[NA_integer_, , drop = FALSE]
      for (name in fields[[from]]) row[[name]] <- input[[name]]
      row$instance_id <- instance %||% NA_character_
      row$status <- status
      row$reason <- reason
      row
    }
    any_display <- !is.null(context)
    if (needs_display && !any_display) {
      value <- miss(if (is.null(request$size)) "missing_dimensions" else
        if (is.null(failure)) "panel_measurement_unavailable" else failure,
        if (is.null(failure)) "unavailable" else "failed")
    } else {
      panel_number <- input[["PANEL"]] %||% facet_number
      candidates <- ggchord_query_match(reference, from, input)
      value <- if (any_display) {
        ggchord_query_project(candidates, context, from, input, axis, panel_number, facet_selector)
      } else candidates
      if (!nrow(value)) value <- miss("not_present") else {
        gap <- if ("is_gap" %in% names(value)) to == "source" & value$is_gap else rep(FALSE, nrow(value))
        value$status <- ifelse(gap, "unavailable", "available")
        value$reason <- ifelse(gap, "gap", NA_character_)
      }
    }
    value$input_row <- rep(i, nrow(value))
    value$panel_path <- rep(list(path), nrow(value))
    if (dotplot) value$axis <- rep(axis, nrow(value))
    for (name in setdiff(names(empty), names(value)))
      value[[name]] <- rep(empty[[name]][NA_integer_], nrow(value))
    value[names(empty)]
  })
  result <- if (length(result)) do.call(rbind, result) else empty
  rownames(result) <- NULL
  attr(result, "measurement") <- measurement
  result
}

ggchord_query_clone <- function(plot) {
  if (inherits(plot, "ggchord_fusion")) {
    for (key in names(plot$panels))
      plot$panels[[key]]$plot <- ggchord_query_clone(plot$panels[[key]]$plot)
    return(plot)
  }
  if (inherits(plot, "ggchord")) return(ggchord_render_plot_copy(plot))
  plot
}

ggchord_query_character_reference <- function(data) data.frame(
  accver = as.character(data$accver), instance_id = data$display_row %||% as.character(data$accver),
  source_accver = data$source_accver, source_position = data$source_position,
  pos = data$alignment_column, sequence_unit = data$sequence_unit, is_gap = data$is_gap)

ggchord_query_reference <- function(plot, axis) {
  data <- plot$ggchord_alignment
  if (!is.null(data)) return(ggchord_query_character_reference(data))
  plot <- ggchord_configure_layout(ggchord_resolve_reference(plot))
  seqs <- plot$ggchord$data$seq_data
  layout <- plot$ggchord$sequence_layout
  windows <- layout$windows
  if (is.null(windows)) {
    ids <- layout$order %||% seqs$accver
    if (!is.null(axis)) {
      axis_layout <- plot$ggchord$dotplot_layout
      ids <- axis_layout[[paste0(axis, "_order")]] %||% ids
    }
    windows <- data.frame(accver = ids, instance_id = ids, start = 1,
      end = seqs$length[match(ids, seqs$accver)])
  }
  if (isTRUE(plot$coordinates$ggchord_sankey)) {
    left <- right <- windows
    left$instance_id <- paste0("left:", left$instance_id)
    right$instance_id <- paste0("right:", right$instance_id)
    windows <- rbind(left, right)
  }
  at <- match(windows$accver, seqs$accver)
  original <- ggchord_focus_source(seqs)
  data.frame(accver = windows$accver, instance_id = windows$instance_id,
    source_accver = original$accver[at], source_origin = original$origin[at],
    local_start = windows$start, local_end = windows$end,
    sequence_unit = seqs$sequence_unit[at])
}

ggchord_query_match <- function(reference, space, input) {
  if (space == "display") return(reference)
  characters <- "is_gap" %in% names(reference)
  from_source <- space == "source"
  accession <- if (from_source) input$source_accver else input$accver
  position <- if (from_source) input$source_position else input$pos
  if (is.na(accession) || is.na(position)) return(reference[FALSE, ])
  id <- if (from_source) reference$source_accver else reference$accver
  if (characters) {
    coordinates <- if (from_source) reference$source_position else reference$pos
    return(reference[which(id == accession & coordinates == position), , drop = FALSE])
  }
  local <- if (from_source) position - reference$source_origin + 1 else rep(position, nrow(reference))
  keep <- which(id == accession & local >= reference$local_start & local <= reference$local_end)
  out <- reference[keep, , drop = FALSE]
  out$pos <- local[keep]
  out$source_position <- out$pos + out$source_origin - 1
  out
}

ggchord_query_transform <- function(leaf, xy, panel, axis = NULL) {
  built <- leaf$built
  index <- built$layout$layout[match(panel, as.integer(built$layout$layout$PANEL)), , drop = FALSE]
  scales <- list(x = built$layout$panel_scales_x[[index$SCALE_X]],
    y = built$layout$panel_scales_y[[index$SCALE_Y]])
  transformed <- lapply(c("x", "y"), function(name) {
    value <- xy[, if (name == "x") 1L else 2L]
    scale <- scales[[name]]
    if (isTRUE(scale$is_discrete())) scale$map(value) else scale$transform(as.numeric(value))
  })
  mapped <- built$layout$coord$transform(data.frame(x = transformed[[1L]], y = transformed[[2L]]),
    built$layout$panel_params[[panel]])
  value <- as.matrix(mapped[c("x", "y")])
  if (!is.null(axis)) value <- value[, axis, drop = FALSE]
  value
}

ggchord_query_project <- function(candidates, leaf, from, input, axis, facet, selector = NULL) {
  built <- leaf$built
  layout <- leaf$layout
  panel_numbers <- as.integer(built$layout$layout$PANEL)
  if (!is.null(facet)) {
    if (!is.numeric(facet) || length(facet) != 1L || is.na(facet) || !facet %in% panel_numbers)
      ggchord_stop("map_ggchord_position(): facet panel does not exist")
    panel_numbers <- facet
  }
  grid <- built$layout$layout
  keys <- setdiff(names(grid), c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y", "COORD"))
  facet_key <- function(panel) as.list(grid[match(panel, as.integer(grid$PANEL)), keys, drop = FALSE])
  if (!is.null(selector)) {
    panel_numbers <- panel_numbers[vapply(panel_numbers, function(panel)
      ggchord_fusion_facet_match(facet_key(panel), selector), logical(1))]
    if (!length(panel_numbers))
      ggchord_stop("map_ggchord_position(): no facet matches; available facets: ",
        paste(vapply(as.list(grid[keys]), function(x) paste(unique(x), collapse = "/"),
          character(1)), collapse = ", "))
  }
  inverse <- from == "display"
  parts <- lapply(panel_numbers, function(panel) {
    projected <- layout$panel_references[[panel]] %||% layout
    frame <- if (!is.null(axis)) layout$axis_panels[[panel]] %||% projected$axis_layout %||% layout$axis_layout
    if (!is.null(axis)) {
      selected <- frame$axis == axis
      if ("PANEL" %in% names(frame)) selected <- selected & as.integer(frame$PANEL) == panel
      frame <- frame[selected, , drop = FALSE]
    }
    rows <- lapply(seq_len(nrow(candidates)), function(i) {
      candidate <- candidates[i, , drop = FALSE]
      id <- candidate$instance_id
      if (isTRUE(layout$character)) {
        at <- which(layout$character_data$accver == candidate$accver &
          layout$character_data$alignment_column == candidate$pos)
        if (length(at) != 1L) return(NULL)
        cell <- layout$character_data[at, , drop = FALSE]
        xy <- ggchord_query_transform(leaf, cbind(cell$display_column, cell$display_row), panel)
        if (inverse && (!all(is.finite(c(input$x, input$y))) ||
            sum((c(input$x, input$y) - xy)^2) > 1e-16)) return(NULL)
      } else {
        key <- if (is.null(layout$display_instances)) id else
          layout$display_instances$key[match(id, layout$display_instances$instance_id)]
        if (is.na(key)) return(NULL)
        if (is.null(axis) && !key %in% projected$seqs) return(NULL)
        if (!is.null(axis) && !key %in% frame$accver) return(NULL)
        project <- function(pos) {
          xy <- if (is.null(axis)) ggchord_panel_backbone_frame(projected, key, pos)$xy else {
            value <- ggchord_dotplot_project(frame, rep(key, length(pos)), pos)
            if (axis == "x") cbind(value, 0) else cbind(0, value)
          }
          ggchord_query_transform(leaf, xy, panel, axis)
        }
        if (inverse) {
          limits <- c(1, candidate$local_end - candidate$local_start + 1)
          point <- if (is.null(axis)) c(input$x, input$y) else input[[axis]]
          samples <- if (!is.null(axis)) 2L else if (isTRUE(projected$collinear)) {
            if (isTRUE(projected$sequence_reference$curvature[[key]] != 0)) 33L else 2L
          } else {
            max(2L, length(projected$sequence_reference$refs[[key]]$angles) + 2L)
          }
          position <- ggchord_query_inverse(point, limits, samples, project)
          if (!is.finite(position)) return(NULL)
          candidate$pos <- candidate$local_start + position - 1
          candidate$source_position <- candidate$source_origin + candidate$pos - 1
        }
        xy <- project(candidate$pos - candidate$local_start + 1)
      }
      if (any(!is.finite(xy) | xy < -1e-10 | xy > 1 + 1e-10)) return(NULL)
      candidate$x <- if (identical(axis, "y")) NA_real_ else xy[1L]
      candidate$y <- if (identical(axis, "x")) NA_real_ else xy[length(xy)]
      candidate$PANEL <- panel
      candidate$facet_key <- list(facet_key(panel))
      candidate
    })
    ggchord_rbind_fill(rows)
  })
  ggchord_rbind_fill(parts) %||% data.frame()
}

ggchord_query_inverse <- function(point, limits, samples, project) {
  if (any(!is.finite(point))) return(NA_real_)
  if (limits[1L] == limits[2L]) return(if
    (sum((project(limits[1L]) - point)^2) <= 1e-16) limits[1L] else NA_real_)
  positions <- seq(limits[1L], limits[2L], length.out = samples)
  path <- project(positions)
  a <- path[-nrow(path), , drop = FALSE]
  b <- path[-1L, , drop = FALSE]
  toward <- sweep(a, 2L, point, function(x, y) y - x)
  direction <- b - a
  squared <- rowSums(direction^2)
  fractions <- pmax(0, pmin(1, rowSums(toward * direction) /
    ifelse(squared > 0, squared, 1)))
  distances <- rowSums((toward - direction * fractions)^2)
  if (!any(is.finite(distances))) return(NA_real_)
  at <- which.min(distances)
  position <- positions[at] + fractions[at] * (positions[at + 1L] - positions[at])
  error <- function(pos) sum((project(pos) - point)^2)
  if (isTRUE(error(position) <= 1e-16)) return(position)
  if (samples > 2L) {
    bounds <- positions[c(max(1L, at - 1L), min(length(positions), at + 2L))]
    fraction <- stats::optimize(function(t) error(bounds[1L] + diff(bounds) * t),
      c(0, 1), tol = 1e-12)
    if (fraction$objective <= 1e-16) return(bounds[1L] + diff(bounds) * fraction$minimum)
  }
  NA_real_
}
