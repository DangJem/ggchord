#' Register an independent character plot for genomic composition
#'
#' Attaches an explicit source/column map to an ordinary Cartesian ggplot.
#' It remains editable and independently printable. Fusion and magnify rebuild
#' the plot at its final physical size, check registered cells against actual
#' geom_residue positions/labels, and project exact source boundaries. No HSP
#' interpolation or implicit alignment is performed.
#'
#' @param plot A single-panel Cartesian ggplot with geom_residue layers mapping
#'   x = alignment_column, y = accver, label = symbol. Identity/reverse x scales
#'   and a discrete y scale are supported. Additional independent annotations
#'   are allowed. Facets, nonlinear scales, flipped coordinates and position
#'   adjustments moving registered cells are rejected.
#'   geom_residue(line_width = ...) is supported. Each selected source interval
#'   must lie within one displayed line; split cross-line selections explicitly.
#' @param data The displayed rows of [as_ggchord_alignment()],
#'   [read_aligned_fasta()] or translate_ggchord_cds()$nucleotides/residues.
#'   One cell per accver/alignment_column; source coordinates must be unique
#'   within each row. Keep original columns when filtering a window. Every
#'   selected source base must be present; split spliced regions explicitly.
#'   Origin-wrapping magnify regions require explicit circular source metadata. Gap columns remain visible but have no source coordinate.
#' @return The ordinary editable plot with a source registration suitable for
#'   [ggchord_fusion()]. Select a display instance if several alignment rows
#'   represent the same source interval.
#'   Registrations assert supplied provenance; source_verified is not an
#'   independent biological validation. Exported source layouts retain data.
#' @export
#' @examples
#' library(ggplot2)
#'
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), alphabet = "DNA")
#' p <- ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(aes(fill = symbol)) + scale_fill_residue("DNA")
#' registered <- register_ggchord_alignment(p, cells)
register_ggchord_alignment <- function(plot, data) {
  ggchord_local_console_errors()
  caller <- "register_ggchord_alignment()"
  if (!inherits(plot, "ggplot") || inherits(plot, "ggchord") ||
      !inherits(plot$coordinates, "CoordCartesian") || inherits(plot$coordinates, "CoordFlip") ||
      !inherits(plot$facet, "FacetNull"))
    ggchord_stop(caller, ": use an ordinary single-panel Cartesian character plot")
  ggchord_require_columns(data, c("accver", "alignment_column", "source_accver", "source_position",
    "symbol", "is_gap", "strand", "sequence_unit", "alphabet"), caller)
  if (!nrow(data)) ggchord_stop(caller, ": data must contain displayed character cells")
  for (name in c("accver", "source_accver", "symbol", "strand", "sequence_unit", "alphabet"))
    if (anyNA(data[[name]]) || any(!nzchar(as.character(data[[name]]))))
      ggchord_stop(caller, ": missing cell metadata: ", name)
  ggchord_project_integer(data$alignment_column, "alignment_column", caller)
  if (!is.logical(data$is_gap) || anyNA(data$is_gap) ||
      any(data$is_gap != is.na(data$source_position)) || any(!data$strand %in% c("+", "-")) ||
      any(!data$sequence_unit %in% c("bp", "aa")) || any(nchar(as.character(data$symbol)) != 1L) ||
      any(data$is_gap != (data$symbol %in% c("-", "."))))
    ggchord_stop(caller, ": inconsistent gap, strand, unit or character metadata")
  ggchord_project_integer(data$source_position[!data$is_gap], "source_position", caller)
  if (any(c(data$alignment_column, data$source_position[!data$is_gap]) > 2^53 - 1))
    ggchord_stop(caller, ": coordinates must be exact positive integers")
  for (id in unique(as.character(data$accver))) {
    d <- data[as.character(data$accver) == id, , drop = FALSE]
    alphabet <- unique(as.character(d$alphabet))
    allowed <- if (length(alphabet) == 1L) switch(alphabet,
      DNA = "ACGTRYSWKMBDHVN-.", RNA = "ACGURYSWKMBDHVN-.", AA = "ACDEFGHIKLMNPQRSTVWYBXZJUO*-.") else NULL
    if (is.null(allowed) || any(!toupper(d$symbol) %in% strsplit(allowed, "", fixed = TRUE)[[1L]]) ||
        any(d$sequence_unit != if (alphabet == "AA") "aa" else "bp"))
      ggchord_stop(caller, ": inconsistent declared alphabet for '", id, "'")
    if (anyDuplicated(d$alignment_column) || anyDuplicated(d$source_position[!d$is_gap]) ||
        length(unique(d$source_accver)) != 1L || length(unique(d$strand)) != 1L ||
        length(unique(d$sequence_unit)) != 1L || all(d$is_gap))
      ggchord_stop(caller, ": ambiguous source/column registration for '", id, "'")
  }
  plot$ggchord_alignment <- as.data.frame(data)
  # Build immediately for useful feedback; also repeat after user edits at render.
  ggchord_alignment_layout(ggplot_build(plot), data)
  plot
}

ggchord_alignment_layout <- function(built, data) {
  caller <- "Registered character panel"
  coord <- built$layout$coord
  if (!inherits(coord, "CoordCartesian") || inherits(coord, "CoordFlip") ||
      !inherits(built$plot$facet, "FacetNull")) ggchord_stop(caller, ": requires single-panel Cartesian coordinates")
  xs <- built$layout$panel_scales_x[[1L]]; ys <- built$layout$panel_scales_y[[1L]]
  if (xs$is_discrete() || !xs$get_transformation()$name %in% c("identity", "reverse") || !ys$is_discrete())
    ggchord_stop(caller, ": use identity/reverse x and discrete accession y scales")
  layers <- which(vapply(built$plot$layers, function(l) inherits(l$geom, "GeomResidue"), logical(1)))
  if (!length(layers)) ggchord_stop(caller, ": geom_residue is required")
  if (any(!vapply(built$plot$layers[layers], function(layer)
      inherits(layer$position, "PositionIdentity"), logical(1))))
    ggchord_stop(caller, ": residue Position must be identity")
  wrapping <- lapply(built$plot$layers[layers], function(l) l$ggchord_residue_wrap)
  wrapped <- !vapply(wrapping, is.null, logical(1))
  data$display_column <- data$alignment_column
  data$display_row <- as.character(data$accver)
  if (any(wrapped)) {
    if (!all(wrapped) || !all(vapply(wrapping, identical, logical(1), wrapping[[1L]])))
      ggchord_stop(caller, ": registered residue layers must share the same line_width and origin")
    specification <- wrapping[[1L]]
    display <- ggchord_residue_wrap(data$alignment_column, data$accver,
      specification$width, specification$origin)
    data$display_column <- display$column
    data$display_row <- as.character(display$row)
  }
  observed <- do.call(rbind, lapply(built$data[layers], function(d) d[c("x", "y", "label")]))
  key <- function(x, y, label) paste(x, y, label, sep = "\r")
  expected <- key(xs$transform(data$display_column), ys$map(data$display_row), data$symbol)
  positions <- coord$transform(observed, built$layout$panel_params[[1L]])
  if (!setequal(expected, key(observed$x, observed$y, observed$label)) ||
      any(!is.finite(positions$x) | !is.finite(positions$y) |
        positions$x < 0 | positions$x > 1 | positions$y < 0 | positions$y > 1))
    ggchord_stop(caller, ": registered cells no longer match rendered residue positions/labels")
  ids <- unique(data$display_row)
  list(character = TRUE, seqs = ids, character_data = data,
    display_instances = data.frame(instance_id = ids,
      accver = as.character(data$accver[match(ids, data$display_row)]), key = ids),
    sequence_unit = unique(data$sequence_unit))
}

ggchord_alignment_interval <- function(panel, id, entity, start, end, instance, point = FALSE) {
  d <- panel$layout$character_data
  candidates <- list()
  for (row in unique(d$display_row[d$source_accver == entity])) {
    if (!is.null(instance) && row != instance) next
    a <- d[d$display_row == row & !d$is_gap, , drop = FALSE]
    lo <- min(start, end); hi <- max(start, end)
    a <- a[a$source_position >= lo & a$source_position <= hi, , drop = FALSE]
    if (!nrow(a) || nrow(a) != hi - lo + 1 || min(a$source_position) != lo || max(a$source_position) != hi) next
    a <- a[order(a$alignment_column), , drop = FALSE]
    sign <- if (a$strand[1L] == "+") 1 else -1
    if (nrow(a) > 1L && any(diff(a$source_position) != sign)) next
    ends <- match(c(start, end), a$source_position)
    direction <- sign * if (end >= start) 1 else -1
    # Source boundaries flank exact endpoint cells; interior gap columns are retained.
    x <- if (point) a$display_column[ends[1L]] else
      a$display_column[ends] + c(-.5, .5) * direction
    candidates[[length(candidates) + 1L]] <- list(panel = id, key = row, instance = row,
      entity = entity, origin = 1, start = start, end = end, unit = a$sequence_unit[1L],
      character_x = x, character_y = row, point = point, facet = 1L, facet_key = list(),
      reference_length = if (length(unique(a$source_length)) == 1L) unique(a$source_length),
      reference_namespace = if (length(unique(a$reference_namespace)) == 1L) unique(a$reference_namespace),
      fragment = identical(unique(a$source_scope), "fragment"))
  }
  candidates
}

ggchord_alignment_front <- function(panel, endpoint) {
  built <- panel$built
  xs <- built$layout$panel_scales_x[[1L]]; ys <- built$layout$panel_scales_y[[1L]]
  xy <- data.frame(x = xs$transform(endpoint$character_x), y = ys$map(endpoint$character_y))
  projected <- built$layout$coord$transform(xy, built$layout$panel_params[[1L]])
  lower <- upper <- xy
  lower$y <- lower$y - .5; upper$y <- upper$y + .5
  lower <- built$layout$coord$transform(lower, built$layout$panel_params[[1L]])
  upper <- built$layout$coord$transform(upper, built$layout$panel_params[[1L]])
  if (any(!is.finite(c(projected$x, projected$y))) ||
      any(projected$x < 0 | projected$x > 1 | projected$y < 0 | projected$y > 1))
    ggchord_stop("Registered character panel: selected source boundary is outside visible limits")
  bounds <- panel$context$panel_bounds_mm
  data.frame(position_boundary = if (isTRUE(endpoint$point)) NA_real_ else
      if (endpoint$end >= endpoint$start) c(endpoint$start - 1, endpoint$end) else
        c(endpoint$start, endpoint$end - 1),
    source_position = if (isTRUE(endpoint$point)) endpoint$start else NA_real_, facet_panel = 1L,
    local_x = endpoint$character_x, local_y = as.numeric(xy$y),
    x = bounds[1] + projected$x * diff(bounds[1:2]), y = bounds[3] + projected$y * diff(bounds[3:4]),
    normal_x = 0, normal_y = 1, tangent_x = -1, tangent_y = 0,
    cell_half_height_mm = abs(upper$y - lower$y) * diff(bounds[3:4]) / 2)
}

ggchord_check_character_origin <- function(q, s, caller) {
  for (endpoint in list(q, s)) {
    if (isTRUE(endpoint$fragment))
      ggchord_stop(caller, ": fragment-local coordinates require an explicit source mapping before linking")
  }
}
