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
#' @param data The displayed rows of [as_ggchord_alignment()],
#'   [read_aligned_fasta()] or translate_ggchord_cds()$nucleotides/residues.
#'   One cell per accver/alignment_column; source coordinates must be unique
#'   within each row. Keep original columns when filtering a window. Every
#'   selected source base must be present; split spliced regions explicitly.
#'   Origin-wrapping magnify regions require explicit circular source metadata. Gap columns remain visible but have no source coordinate.
#' @return The ordinary plot with a source registration. Use it as detail_plot
#'   in [magnify_ggchord()] or as an outer [fusion_panel()]. Region detail_instance
#'   selects an exact alignment row accver when several rows share a source.
#'   Registrations assert supplied provenance; source_verified is not an
#'   independent biological validation. Exported source layouts retain data.
#' @export
#' @examples
#' # Explicit synthetic format example.
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), "DNA")
#' p <- ggplot2::ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(aes(fill = symbol), show.legend = FALSE)
#' register_ggchord_alignment(p, cells)
register_ggchord_alignment <- function(plot, data) {
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
  ggchord_alignment_layout(ggplot2::ggplot_build(plot), data)
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
  observed <- do.call(rbind, lapply(built$data[layers], function(d) d[c("x", "y", "label")]))
  key <- function(x, y, label) paste(x, y, label, sep = "\r")
  expected <- key(xs$transform(data$alignment_column), ys$map(as.character(data$accver)), data$symbol)
  if (any(!expected %in% key(observed$x, observed$y, observed$label)))
    ggchord_stop(caller, ": registered cells no longer match rendered residue positions/labels")
  ids <- unique(as.character(data$accver))
  list(character = TRUE, seqs = ids, character_data = data,
    display_instances = data.frame(instance_id = ids, accver = ids, key = ids),
    sequence_unit = unique(data$sequence_unit))
}

ggchord_alignment_interval <- function(panel, id, entity, start, end, instance) {
  d <- panel$layout$character_data
  candidates <- list()
  for (row in unique(as.character(d$accver[d$source_accver == entity]))) {
    if (!is.null(instance) && row != instance) next
    a <- d[as.character(d$accver) == row & !d$is_gap, , drop = FALSE]
    lo <- min(start, end); hi <- max(start, end)
    a <- a[a$source_position >= lo & a$source_position <= hi, , drop = FALSE]
    if (!nrow(a) || nrow(a) != hi - lo + 1 || min(a$source_position) != lo || max(a$source_position) != hi) next
    a <- a[order(a$alignment_column), , drop = FALSE]
    sign <- if (a$strand[1L] == "+") 1 else -1
    if (nrow(a) > 1L && any(diff(a$source_position) != sign)) next
    ends <- match(c(start, end), a$source_position)
    direction <- sign * if (end >= start) 1 else -1
    # Source boundaries flank exact endpoint cells; interior gap columns are retained.
    x <- a$alignment_column[ends] + c(-.5, .5) * direction
    candidates[[length(candidates) + 1L]] <- list(panel = id, key = row, instance = row,
      entity = entity, origin = 1, start = start, end = end, unit = a$sequence_unit[1L],
      character_x = x, character_y = row)
  }
  candidates
}

ggchord_alignment_front <- function(panel, endpoint) {
  built <- panel$built
  xs <- built$layout$panel_scales_x[[1L]]; ys <- built$layout$panel_scales_y[[1L]]
  xy <- data.frame(x = xs$transform(endpoint$character_x), y = ys$map(endpoint$character_y))
  projected <- built$layout$coord$transform(xy, built$layout$panel_params[[1L]])
  if (any(!is.finite(c(projected$x, projected$y))) ||
      any(projected$x < 0 | projected$x > 1 | projected$y < 0 | projected$y > 1))
    ggchord_stop("Registered character panel: selected source boundary is outside visible limits")
  bounds <- panel$context$panel_bounds_mm
  data.frame(position_boundary = if (endpoint$end >= endpoint$start) c(endpoint$start - 1, endpoint$end) else c(endpoint$start, endpoint$end - 1),
    local_x = endpoint$character_x, local_y = as.numeric(xy$y),
    x = bounds[1] + projected$x * diff(bounds[1:2]), y = bounds[3] + projected$y * diff(bounds[3:4]),
    normal_x = 0, normal_y = 1, tangent_x = -1, tangent_y = 0)
}
