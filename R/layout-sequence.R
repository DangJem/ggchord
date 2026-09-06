# Sequence allocation, references and track geometry.
ggchord_layout_sequence_step <- function(context) evalq({
  n <- length(seqs)

  # ====================================================================
  # Step 1: compute angle allocation
  # ====================================================================
  total_circ <- 2 * pi

  total_gap_prop <- sum(seq_gap)

  if (total_gap_prop >= 1) {
    ggchord_stop("The sum of seq_gap cannot exceed 1 (no space left for sequences)")
  }

  seq_total_prop <- 1 - total_gap_prop
  sum_lens <- sum(lens)
  theta <- (lens / sum_lens) * total_circ * seq_total_prop
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

  gene_track_radius <- function(gene, sid, strand) {
    side <- if (".feature_stack_side" %in% names(gene)) {
      as.character(gene[[".feature_stack_side"]])
    } else if (strand == "+") {
      "inside"
    } else {
      "outside"
    }
    lane <- if (".feature_stack_lane" %in% names(gene)) {
      as.numeric(gene[[".feature_stack_lane"]])
    } else 0
    spacing <- if (".feature_stack_spacing" %in% names(gene)) {
      as.numeric(gene[[".feature_stack_spacing"]])
    } else 0
    direction <- if (identical(side, "inside")) -1 else 1
    seqRadius[sid] + direction *
      (geneGap[[sid]][strand] + lane * spacing)
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
}, envir = context)
