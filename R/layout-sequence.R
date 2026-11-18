#' Declare the arrangement of sequence backbones
#'
#' A plot has one sequence layout, independent of its drawing layers. Add the
#' declaration before or after the layers. It also applies when the backbone
#' is hidden. A named instance inherits its source orientation and gap unless
#' overridden by instance ID. A single circular map uses the opening and
#' direction in [coord_circular()]; [coord_dotplot()] owns both of its axes.
#' @param order Source IDs in display order. Repeats create display instances;
#'   use names such as c(A_high='A', B='B', A_low='A') for stable instance IDs.
#'   Repeated-entity HSPs require a unique between route on each link layer.
#' @param orientation Direction (1 or -1). A named vector selects source or display instance IDs; unspecified sequences use 1.
#' @param gap Optional numeric. Gap proportion between sequences, default
#'   0.03 in chord coordinates. A named vector can set the gap after each
#'   sequence; omitted IDs use the coordinate default. In [coord_collinear()], a positive distance
#'   between consecutive backbones as a proportion of the longest backbone's
#'   width, default 0.3. The final sequence's gap is unused in linear layouts.
#'   Segmented linear displays measure it against the longest display row;
#'   windows in one row must request the same gap, and only gaps between rows
#'   affect vertical placement.
#' @param length_mode `"proportional"` (default) retains source-length ratios;
#'   `"equal"` gives each displayed sequence the same angular allocation or
#'   linear span without modifying source lengths. Shared absolute scale bars
#'   are unavailable in multi-sequence equal mode; local axes retain source units.
#' @return A layout declaration to add to a ggchord plot with `+`.
#' @export
#' @examples
#' reference <- data.frame(accver = c("A", "B"), length = c(100, 200))
#' ggchord(reference) +
#'   geom_seq() +
#'   layout_sequences(order = c("B", "A"), orientation = c(B = -1), gap = 0.2,
#'                    length_mode = "equal") +
#'   coord_collinear()
layout_sequences <- function(order = NULL, orientation = NULL, gap = NULL,
    length_mode = c("proportional", "equal")) {
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!is.null(order) && (!is.character(order) || !length(order) || anyNA(order) || any(!nzchar(order))))
    ggchord_stop("layout_sequences(): order must contain nonempty sequence IDs")
  if (!is.null(names(order)) && (anyNA(names(order)) || any(!nzchar(names(order))) ||
      anyDuplicated(names(order))))
    ggchord_stop("layout_sequences(): display instance names must be nonempty and unique")
  for (value in list(orientation, gap)) {
    if (anyDuplicated(names(value)))
      ggchord_stop("layout_sequences(): a sequence or instance cannot be configured twice")
  }
  structure(list(order = order, orientation = orientation, gap = gap,
    length_mode = match.arg(length_mode)), class = "ggchord_sequence_layout")
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.ggchord_sequence_layout <- function(object, plot, ...) {
  if (!inherits(plot, "ggchord")) ggchord_stop("Add layout_sequences() to a ggchord plot")
  if (!is.null(plot$ggchord$sequence_layout))
    ggchord_stop("Sequence layout already declared; supply one layout_sequences() per plot")
  plot$ggchord$sequence_layout <- object
  plot$ggchord$ref <- new.env(parent = emptyenv())
  plot$ggchord$layout <- NULL
  plot
}

# Sequence allocation, references and track geometry.
ggchord_layout_sequence_step <- function(
    axisGap, axisMajLen, axisMinLen, circular, geometry_cache, labelOffset, lens, orientation,
    rotation, seq_curvature, seq_gap, seqRadius, seqs, length_mode) {
  n <- length(seqs)

  # ====================================================================
  # Step 1: compute angle allocation
  # ====================================================================
  total_circ <- 2 * pi

  total_gap_prop <- sum(seq_gap)

  if (total_gap_prop >= 1) {
    ggchord_stop("The sum of layout_sequences(gap) cannot exceed 1 (no space left for sequences)")
  }

  seq_total_prop <- 1 - total_gap_prop
  weights <- if (identical(length_mode, "equal")) rep(1, length(lens)) else lens
  theta <- (weights / sum(weights)) * total_circ * seq_total_prop
  gap_rads <- total_circ * seq_gap

  # Compute the start and end angles of each sequence
  starts <- numeric(n)
  names(starts) <- seqs
  starts[1] <- 0

  if (n > 1) {
    for (i in 2:n) {
      starts[i] <- starts[i - 1] + theta[i - 1] + gap_rads[i - 1]
    }
  }
  ends <- starts + theta
  names(ends) <- seqs
  if (isTRUE(circular) && n == 1L && orientation[seqs[1L]] == -1) {
    starts[1L] <- -theta[1L]
    ends[1L] <- 0
  }

  # ====================================================================
  # Step 2: convert to radians
  # ====================================================================
  rot_rad <- rotation * pi / 180

  # ====================================================================
  # Step 3: generate reference paths and map_to_curve
  # ====================================================================
  nSeg <- 500
  nRef <- 2000

  # High-resolution reference path for each sequence
  reference_key <- list(
    seqs = seqs, starts = starts, ends = ends,
    radius = seqRadius, curvature = seq_curvature
  )
  if (!is.null(geometry_cache) &&
      identical(geometry_cache$reference_key, reference_key)) {
    seq_refs <- geometry_cache$seq_refs
  } else {
    seq_refs <- lapply(seqs, function(id) {
      path <- generate_curvature_path(
        starts[id], ends[id], seqRadius[id], seq_curvature[id], n_points = nRef
      )
      angles <- seq(starts[id], ends[id], length.out = nRef)
      list(path = path, angles = angles, r0 = seqRadius[id])
    })
    names(seq_refs) <- seqs
    if (!is.null(geometry_cache)) {
      geometry_cache$reference_key <- reference_key
      geometry_cache$seq_refs <- seq_refs
    }
  }

  gene_track_radius <- function(gene, sid, strand, angle) {
    signed_outward <- if (".normal_offset" %in% names(gene)) {
      as.numeric(gene[[".normal_offset"]])
    } else 0
    if (!is.finite(signed_outward) || signed_outward == 0) {
      return(unname(seqRadius[sid]))
    }
    ref <- seq_refs[[sid]]
    base <- map_to_curve(angle, seqRadius[sid], ref)
    plus <- map_to_curve(angle, seqRadius[sid] + 1e-5, ref)
    left_normal <- plus - base
    # Convert the stable signed-outward Position contract to the curve
    # mapper's left-normal radius convention at this exact genomic location.
    # This remains valid for signed and non-unit sequence curvature.
    radius_sign <- if (sum(left_normal * base) >= 0) 1 else -1
    unname(seqRadius[sid]) + radius_sign * signed_outward
  }

  # ====================================================================
  # Step 4: generate sequence arcs (outer layer)
  # ====================================================================
  seq_arcs <- stats::setNames(lapply(seqs, function(id) {
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id], seq_curvature[id], nSeg
    )
    path_data$accver <- id
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    attr(path_data, "ggchord_path_direction") <- unname(orientation[id])
    path_data
  }), seqs)

  axis_reference_x <- unlist(lapply(seq_arcs, `[[`, "x"), use.names = FALSE)
  axis_reference_y <- unlist(lapply(seq_arcs, `[[`, "y"), use.names = FALSE)
  axis_units_per_inch <- ggchord_device_units_per_inch(
    axis_reference_x, axis_reference_y
  )
  axisGap <- axisGap * axis_units_per_inch
  axisMajLen <- axisMajLen * axis_units_per_inch
  axisMinLen <- axisMinLen * axis_units_per_inch
  labelOffset <- labelOffset * axis_units_per_inch

  # ====================================================================
  list(n = n, starts = starts, ends = ends, rot_rad = rot_rad, nSeg = nSeg, seq_refs = seq_refs,
    gene_track_radius = gene_track_radius, seq_arcs = seq_arcs,
    axis_units_per_inch = axis_units_per_inch, axisGap = axisGap, axisMajLen = axisMajLen,
    axisMinLen = axisMinLen, labelOffset = labelOffset)
}
