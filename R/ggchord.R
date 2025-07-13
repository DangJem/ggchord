globalVariables(c(
  "x", "y", "group", "pident", "fill", "strand", "anno", "seq_id",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "label", "label_x", "label_y", "size"
))

#' ggchord: A ggplot2-based tool for multi-sequence alignment chord plots
#'
#' ggchord is used to draw chord plots containing multiple sequences, which can display alignment relationships between sequences and gene annotation information.
#' ggchord supports customizing various parameters such as sequence arrangement, colors, ribbon styles, and gene arrow styles, making it suitable for genome alignment visualization.
#'
#' @param seq_data data.frame/tibble, required. A data frame containing basic sequence information, must include columns:
#'   - seq_id: Unique sequence identifier (character)
#'   - length: Sequence length (numeric, > 0)
#' @param ribbon_data data.frame/tibble, optional. BLAST alignment result data frame, must include columns:
#'   - qaccver: Query sequence ID (matching seq_id)
#'   - saccver: Subject sequence ID (matching seq_id)
#'   - length: Alignment length
#'   - pident: Percentage of sequence identity (0-100)
#'   - qstart/qend: Start/end positions of the query sequence in the alignment
#'   - sstart/send: Start/end positions of the subject sequence in the alignment
#' @param gene_data data.frame/tibble, optional. Gene annotation data frame, must include columns:
#'   - seq_id: ID of the associated sequence (matching seq_id)
#'   - start/end: Gene start/end positions (numeric)
#'   - strand: Strand direction (only "+" or "-")
#'   - anno: Gene annotation name (character)
#' @param title Character. Main title of the plot, default NULL (no title displayed)
#' @param seq_order Character vector, optional. Specifies the drawing order of sequences (must be a subset of seq_id), default follows the order in seq_data
#' @param seq_labels Character vector/named vector, optional. Sequence labels (length matching the number of sequences or named to match seq_id), default uses seq_id
#' @param seq_orientation Numeric (1 or -1), optional. Sequence direction (1 = forward, -1 = reverse), supports single value/vector/named vector, default 1
#' @param seq_gap Numeric (0 <= x < 0.5), optional. Proportion of gap between sequences, supports single value/vector/named vector, default 0.03
#' @param seq_radius Numeric (> 0), optional. Radius of sequence arcs, supports single value/vector/named vector, default 1.0
#' @param seq_curvature Numeric, optional. Curvature of sequence arcs (0 = straight line, 1 = standard arc, > 1 = more curved), default 1.0
#' @param seq_colors Color vector/named vector, optional. Colors of sequence arcs, default auto-generated based on RColorBrewer Set1
#' @param gene_offset Numeric/vector/list, optional. Radial offset of gene arrows from sequences (positive values outward, negative values inward), supports:
#'   - single value: shared by all sequences/strands
#'   - vector: length matching the number of sequences (assigned by sequence)
#'   - list: named list (elements are single values or vectors with "+"/"-" to distinguish strands), default 0.1
#' @param gene_width Numeric/vector/list, optional. Width of gene arrows, format same as gene_offset, default 0.05
#' @param gene_label_show Logical. Whether to display gene labels, default FALSE
#' @param gene_label_rotation Numeric/vector/list, optional. Rotation angle (degrees) of gene labels, format same as gene_offset, default 0
#' @param gene_label_size Numeric. Font size of gene labels, default 2.5
#' @param gene_label_radial_offset Numeric/vector/list, optional. Radial offset of gene labels relative to arrows, format same as gene_offset, default 0
#' @param gene_label_circum_offset Numeric/vector/list, optional. Circumferential offset proportion of gene labels along sequences, format same as gene_offset, default 0
#' @param gene_label_circum_limit Logical/vector/list, optional. Whether to limit circumferential offset to half the gene length, format same as gene_offset, default TRUE
#' @param gene_color_scheme Character. Color scheme for genes, optional "strand" (by strand direction) or "manual" (by annotation), default "strand"
#' @param gene_colors Color vector, optional. Fill colors for gene arrows, format depends on gene_color_scheme:
#'   - "strand": named vector ("+"/"-"), unnamed vector of length 1/2
#'   - "manual": named vector (matching anno), unnamed vector (recycled), default auto-generated
#' @param gene_order Character vector, optional. Display order of genes in the legend (matching anno), default follows the order in data
#' @param ribbon_color_scheme Character. Color scheme for ribbons, optional "pident" (gradient by identity), "query" (by query sequence), "single" (uniform color), default "pident"
#' @param ribbon_colors Color vector, optional. Color parameters for ribbons:
#'   - "single": uniform color
#'   - "query": color vector matching seq_id
#'   - "pident": color gradient (at least 2 colors), default blue-to-yellow gradient
#' @param ribbon_alpha Numeric (0-1). Transparency of ribbons, default 0.35
#' @param ribbon_ctrl_point Vector/list, optional. Control points for Bézier curves (adjust ribbon shape):
#'   - vector: length 2 (single control point) or 4 (c1x,c1y,c2x,c2y for two control points)
#'   - list: each element is a sublist with 1-2 control points, default c(0,0)
#' @param ribbon_gap Numeric/vector/named vector, optional. Radial distance between sequences and ribbons, default 0.15
#' @param axis_gap Numeric/vector/named vector, optional. Radial distance between sequences and axes, default 0.04
#' @param axis_tick_major_number Integer/vector/named vector, optional. Number of major ticks, default 5
#' @param axis_tick_major_length Numeric/vector/named vector, optional. Length proportion of major ticks, default 0.02
#' @param axis_tick_minor_number Integer/vector/named vector, optional. Number of minor ticks per major tick, default 4
#' @param axis_tick_minor_length Numeric/vector/named vector, optional. Length proportion of minor ticks, default 0.01
#' @param axis_label_size Numeric/vector/named vector, optional. Font size of axis labels, default 3
#' @param axis_label_offset Numeric/vector/named vector, optional. Offset proportion of labels relative to ticks, default 1.5
#' @param axis_label_orientation Character/numeric/vector, optional. Orientation of axis labels:
#'   - "horizontal": horizontal
#'   - numeric: rotation angle (degrees)
#'   - vector: length matching the number of sequences or named vector (matching seq_id), default "horizontal"
#' @param rotation Numeric. Overall rotation angle of the plot (degrees), default 45
#' @param panel_margin Numeric/list, optional. Margin around the plot panel (t=top, r=right, b=bottom, l=left):
#'   - single value: same margin for all sides
#'   - list: named list (e.g., list(t=1,r=1)), default 0
#' @param show_legend Logical. Whether to display legends, default TRUE
#' @param show_axis Logical. Whether to display axes and ticks, default TRUE
#' @param debug Logical. Whether to output debugging information (e.g., number of valid ribbons), default FALSE
#'
#' @return A ggplot2 graph object, which can be further adjusted with ggplot2 functions
#' @export
#' @examples
#' # Example code
#' p <- ggchord(
#'   seq_data = seq_data_example,
#'   ribbon_data = ribbon_data_example,
#'   gene_data = gene_data_example
#' )
#' print(p)
#'
#' @import ggplot2
#' @import ggnewscale
#' @import RColorBrewer
#' @import grDevices
#' @import grid


ggchord <- function(
    seq_data,
    ribbon_data = NULL,
    gene_data = NULL,
    title = NULL,
    seq_order = NULL,
    seq_labels = NULL,
    seq_orientation = NULL,
    seq_gap = 0.03,
    seq_radius = 1.0,
    seq_colors = NULL,
    seq_curvature = 1.0,
    gene_offset = 0.1,
    gene_width = 0.05,
    gene_label_show = FALSE,
    gene_label_rotation = 0,
    gene_label_size = 2.5,
    gene_label_radial_offset = 0,
    gene_label_circum_offset = 0,
    gene_label_circum_limit = TRUE,
    gene_color_scheme = c("strand", "manual"),
    gene_colors = NULL,
    gene_order = NULL,
    ribbon_color_scheme = c("pident", "query", "single"),
    ribbon_colors = NULL,
    ribbon_alpha = 0.35,
    ribbon_ctrl_point = c(0,0),
    ribbon_gap = 0.15,
    axis_gap = 0.04,
    axis_tick_major_number = 5,
    axis_tick_major_length = 0.02,
    axis_tick_minor_number = 4,
    axis_tick_minor_length = 0.01,
    axis_label_size = 3,
    axis_label_offset = 1.5,
    axis_label_orientation = "horizontal",
    rotation = 45,
    panel_margin = 0,
    show_legend = TRUE,
    show_axis = TRUE,
    debug = FALSE
) {
  ribbon_color_scheme <- match.arg(ribbon_color_scheme)
  gene_color_scheme <- match.arg(gene_color_scheme)

  required_seq_cols <- c("seq_id", "length")
  if (!all(required_seq_cols %in% colnames(seq_data))) {
    stop("seq_data must contain the following columns: ", paste(required_seq_cols, collapse = ", "))
  }
  if (any(seq_data$length <= 0)) {
    stop("The 'length' values in seq_data must be positive numbers")
  }

  fb_all <- NULL
  if (!is.null(ribbon_data)) {
    required_ribbon_cols <- c("qaccver", "saccver", "length", "pident",
                              "qstart", "qend", "sstart", "send")
    if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
      stop("ribbon_data must contain the following columns: ", paste(required_ribbon_cols, collapse = ", "))
    }
    fb_all <- ribbon_data # Use the preprocessed alignment data directly
    if (nrow(fb_all) == 0) warning("No valid alignment data found in ribbon_data")
    if (debug) cat("Number of alignment data rows used:", nrow(fb_all), "\n")
  }

  # Validate gene annotation data (optional)
  if (!is.null(gene_data)) {
    required_gene_cols <- c("seq_id", "start", "end", "strand", "anno")
    if (!all(required_gene_cols %in% colnames(gene_data))) {
      stop("gene_data must contain the following columns: ", paste(required_gene_cols, collapse = ", "))
    }
    if (nrow(gene_data) == 0) warning("No valid gene annotation data found in gene_data")
    if (debug) cat("Number of gene annotation data rows used:", nrow(gene_data), "\n")

    # Validate gene_order parameter
    if (!is.null(gene_order)) {
      unknown_genes <- setdiff(gene_order, unique(gene_data$anno))
      if (length(unknown_genes) > 0) {
        stop("gene_order contains unknown gene annotations: ", paste(unknown_genes, collapse = ", "))
      }
    }
  }

  # 2. Process sequence information
  seqs <- seq_data$seq_id
  lens <- setNames(seq_data$length, seqs)
  if (length(seqs) < 1) stop("seq_data must contain at least one sequence")

  # Process sequence order
  if (!is.null(seq_order)) {
    if (!all(seq_order %in% seqs)) {
      stop("seq_order contains unknown sequence IDs: ", paste(setdiff(seq_order, seqs), collapse = ", "))
    }
    seqs <- seq_order
    lens <- lens[seqs] # Reorder lengths according to the new sequence order
  }
  n <- length(seqs) # Number of sequences

  # 3. Process sequence-related parameters
  seq_labels <- process_sequence_param(seq_labels, seqs, "seq_labels", default_value = seqs)
  seqRadius <- process_sequence_param(seq_radius, seqs, "seq_radius", default_value = 1.0)
  ribbonGap <- process_sequence_param(ribbon_gap, seqs, "ribbon_gap", default_value = 0.1)
  axisGap <- process_sequence_param(axis_gap, seqs, "axis_gap", default_value = 0.05)
  axisMaj <- process_sequence_param(axis_tick_major_number, seqs, "axis_tick_major_number", default_value = 5)
  axisMajLen <- process_sequence_param(axis_tick_major_length, seqs, "axis_tick_major_length", default_value = 0.02)
  axisMin <- process_sequence_param(axis_tick_minor_number, seqs, "axis_tick_minor_number", default_value = 4)
  axisMinLen <- process_sequence_param(axis_tick_minor_length, seqs, "axis_tick_minor_length", default_value = 0.01)
  axisLabelOrientation <- process_axis_orientation(axis_label_orientation, seqs)
  labelSize <- process_sequence_param(axis_label_size, seqs, "axis_label_size", default_value = 3)
  labelOffset <- process_sequence_param(axis_label_offset, seqs, "axis_label_offset", default_value = 0)
  orientation <- process_sequence_param(seq_orientation, seqs, "seq_orientation", default_value = 1)
  rot_rad <- rotation * pi / 180  # Convert to radians
  # Process gene offset parameters
  geneGap <- process_gene_param(
    param = gene_offset,
    seqs = seqs,
    param_name = "gene_offset",
    default_value = 0.03,
    is_logical = FALSE
  )
  geneWidth <- process_gene_param(
    param = gene_width,
    seqs = seqs,
    param_name = "gene_width",
    default_value = 0.1,
    is_logical = FALSE
  )
  # Process gene label offset parameters
  geneLabelRadialOffset <- process_gene_param(
    param = gene_label_radial_offset,
    seqs = seqs,
    param_name = "gene_label_radial_offset",
    default_value = 0,
    is_logical = FALSE
  )
  geneLabelCircumOffset <- process_gene_param(
    param = gene_label_circum_offset,
    seqs = seqs,
    param_name = "gene_label_circum_offset",
    default_value = 0,
    is_logical = FALSE
  )
  geneLabelCircumLimit <- process_gene_param(
    param = gene_label_circum_limit,
    seqs = seqs,
    param_name = "gene_label_circum_limit",
    default_value = TRUE,
    is_logical = TRUE
  )
  # Process gene label rotation angle parameters
  geneLabelRotation <- process_gene_param(
    param = gene_label_rotation,
    seqs = seqs,
    param_name = "gene_label_rotation",
    default_value = 0,
    is_logical = FALSE
  )

  # Process sequence gap parameters
  seq_gap <- process_sequence_param(seq_gap, seqs, "seq_gap")
  if (any(seq_gap < 0 | seq_gap >= 0.5)) {
    stop("seq_gap must be in the range [0, 0.5)")
  }

  # Process sequence curvature parameters
  seq_curvature <- process_sequence_param(seq_curvature, seqs, "seq_curvature", default_value = 1.0)

  # Process sequence colors
  if (is.null(seq_colors)) {
    pal <- if (length(seqs) <= 9) {
      brewer.pal(length(seqs), "Set1")
    } else {
      colorRampPalette(brewer.pal(9, "Set1"))(length(seqs))
    }
    seq_colors <- setNames(pal, seqs)
  } else {
    seq_colors <- process_sequence_param(seq_colors, seqs, "seq_colors")
  }

  # 4. Process ribbon color parameters (only if ribbon_data is provided)
  if (!is.null(ribbon_data)) {
    if (is.null(ribbon_colors)) {
      ribbon_colors <- switch(ribbon_color_scheme,
                              single = "steelblue",
                              query = {
                                # Lightened version of seq_colors (mixed with white at 'mix' ratio)
                                mix <- 0.5
                                sapply(seq_colors, function(col) {
                                  cols <- col2rgb(col)
                                  light_cols <- cols + (255 - cols) * mix
                                  rgb(light_cols[1,], light_cols[2,], light_cols[3,], maxColorValue = 255)
                                })
                              },
                              pident = c(
                                "#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF",
                                "#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF"
                              )
      )
    }

    # Validate ribbon color parameters
    if (ribbon_color_scheme == "single") {
      singleCol <- if (length(ribbon_colors) > 1) ribbon_colors[[1]] else ribbon_colors
    } else if (ribbon_color_scheme == "query") {
      queryCols <- process_sequence_param(ribbon_colors, seqs, "ribbon_colors")
    } else if (ribbon_color_scheme == "pident") {
      if (length(ribbon_colors) < 2) {
        stop("At least two colors are required for the 'pident' color scheme")
      }
      rampFunc <- colorRampPalette(ribbon_colors)
    }
  }

  # 5. Calculate gene colors (based on gene_color_scheme)
  gene_pal <- NULL
  final_gene_order <- NULL
  if (!is.null(gene_data) && nrow(gene_data) > 0) {
    valid_genes <- gene_data[gene_data$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      # Get unique gene annotations
      unique_anno <- unique(valid_genes$anno)

      # Determine final gene order
      if (!is.null(gene_order)) {
        final_gene_order <- c(gene_order, setdiff(unique_anno, gene_order))
      } else {
        final_gene_order <- unique_anno
      }

      # Process colors according to the scheme
      if (gene_color_scheme == "strand") {
        gene_pal <- process_strand_colors(gene_colors)
      } else if (gene_color_scheme == "manual") {
        gene_pal <- process_manual_colors(gene_colors, unique_anno, gene_order)
      }
    }
  }

  # 6. Calculate sequence radians and gap radians
  total_circ <- 2 * pi  # Total circumference in radians
  total_gap_prop <- sum(seq_gap)  # Total gap proportion

  if (total_gap_prop >= 1) {
    stop("Sum of seq_gap cannot exceed 1 (insufficient space for all sequences)")
  }

  seq_total_prop <- 1 - total_gap_prop  # Total proportion occupied by sequences
  sum_lens <- sum(lens)
  theta <- (lens / sum_lens) * total_circ * seq_total_prop  # Radians for each sequence
  gap_rads <- total_circ * seq_gap  # Radians for each gap (actual angle)

  # 7. Calculate sequence start and end angles
  starts <- numeric(n)
  starts[1] <- 0  # First sequence starts at 0

  if (n > 1) {
    for (i in 2:n) {
      starts[i] <- starts[i-1] + theta[i-1] + gap_rads[i-1]
    }
  }

  ends <- starts + theta
  names(starts) <- names(ends) <- seqs

  # 8. Prepare sequence outer and inner coordinates
  nSeg <- 500  # Number of segments per sequence (controls smoothness)
  seqArcs <- lapply(seqs, function(id) {
    # Outer arc of the sequence (unchanged)
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id], seq_curvature[id], nSeg
    )
    path_data$seq_id <- id
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    path_data
  })

  innerArcs <- lapply(seqs, function(id) {
    # Inner arc for ribbons (correction: seqRadius + ribbonGap → positive outward, negative inward)
    # Principle: When ribbonGap > 0, radius = sequence radius + gap (outward); when ribbonGap < 0, radius = sequence radius + negative value = smaller (inward)
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id] + ribbonGap[id], seq_curvature[id], nSeg
    )
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    path_data
  })
  names(innerArcs) <- seqs


  # Generate high-resolution reference paths for each sequence
  seq_refs <- lapply(seqs, function(id) {
    ref_n <- 2000
    # Reference radius uses the outermost layer: sequence radius
    r0 <- seqRadius[id]
    path <- generate_curvature_path(starts[id], ends[id], r0, seq_curvature[id], n_points = ref_n)
    angles <- seq(starts[id], ends[id], length.out = ref_n)
    list(path = path, angles = angles, r0 = r0)
  })
  names(seq_refs) <- seqs

  # Mapping function
  map_to_curve <- function(angle, radius, ref) {
    # Find the reference point closest to the angle
    idx <- which.min(abs(ref$angles - angle))
    base <- ref$path[idx, ]
    # Calculate tangent and normal
    if (idx < nrow(ref$path)) {
      dx <- ref$path$x[idx+1] - base$x
      dy <- ref$path$y[idx+1] - base$y
    } else {
      dx <- base$x - ref$path$x[idx-1]
      dy <- base$y - ref$path$y[idx-1]
    }
    norm <- c(-dy, dx)
    norm <- norm / sqrt(sum(norm^2))
    # Offset along normal
    offset <- radius - ref$r0
    c(x = base$x + norm[1] * offset,
      y = base$y + norm[2] * offset)
  }


  # Generate tick and label coordinates using map_to_curve
  axisTicks <- do.call(rbind, lapply(seqs, function(id) {
    ref <- seq_refs[[id]]
    # Reference radius: sequence radius + axisGap
    r0 <- ref$r0 - axisGap[id]

    # Major and minor tick positions
    majors <- breakPointsFunc(lens[id], axisMaj[id])
    minors <- unlist(lapply(seq_len(length(majors)-1), function(i) {
      seq(majors[i], majors[i+1], length.out = axisMin[id] + 2)[-c(1, axisMin[id] + 2)]
    }))
    pts <- data.frame(pos = c(majors, minors),
                      is_major = c(rep(TRUE, length(majors)), rep(FALSE, length(minors))))

    # Map each tick point
    do.call(rbind, lapply(seq_len(nrow(pts)), function(j) {
      p <- pts[j,]
      frac <- if (orientation[id] == 1) p$pos / lens[id] else 1 - p$pos / lens[id]
      angle <- starts[id] + frac * (ends[id] - starts[id])

      # Base point, tick tip, label position
      base <- map_to_curve(angle, r0, ref)
      dir <- if (axisGap[id] >= 0) -1 else 1
      len <- if (p$is_major) axisMajLen[id] else axisMinLen[id]
      tip <- map_to_curve(angle, r0 + len * dir, ref)
      lbl <- map_to_curve(angle, r0 + len * (1.5 + labelOffset[id]) * dir, ref)

      data.frame(
        x0 = base[1], y0 = base[2],
        x1 = tip[1], y1 = tip[2],
        label = if (p$is_major) as.character(p$pos) else NA,
        label_x = lbl[1], label_y = lbl[2],
        size = labelSize[id],
        seq_id = id
      )
    }))
  }))


  # Generate full axis lines using map_to_curve
  axisLines <- do.call(rbind, lapply(seqs, function(id) {
    ref <- seq_refs[[id]]
    # Reference radius: sequence radius + axisGap
    r0 <- ref$r0 - axisGap[id]
    # Divide the entire angle range into nSeg segments
    angles <- seq(starts[id], ends[id], length.out = nSeg)
    # Map each angle
    pts <- t(sapply(angles, function(angle) {
      map_to_curve(angle, r0, ref)
    }))
    data.frame(x = pts[,1], y = pts[,2], seq_id = id)
  }))

  # 10. Generate ribbon data (only if ribbon_data is provided)
  allRibbon <- NULL
  if (!is.null(ribbon_data) && !is.null(fb_all) && nrow(fb_all) > 0) {
    ribbons <- list()
    cntValid <- cntInvalid <- 0
    for (i in seq_len(nrow(fb_all))) {
      row <- fb_all[i,]
      q <- row$qaccver
      s <- row$saccver
      if (q == s || !q %in% seqs || !s %in% seqs) {
        cntInvalid <- cntInvalid + 1
        next
      }

      # Correct query sequence ribbon start/end coordinates (using seqRadius + ribbonGap to ensure correct direction)
      q_ref <- seq_refs[[q]]
      q_frac_start <- if (orientation[q] == 1) (row$qstart - 1)/lens[q] else 1 - (row$qstart - 1)/lens[q]
      q_angle_start <- starts[q] + q_frac_start * (ends[q] - starts[q])
      # Ribbon radius = sequence radius + ribbonGap (positive outward, negative inward)
      q_start_coord <- map_to_curve(q_angle_start, seqRadius[q] + ribbonGap[q], q_ref)

      q_frac_end <- if (orientation[q] == 1) (row$qend - 1)/lens[q] else 1 - (row$qend - 1)/lens[q]
      q_angle_end <- starts[q] + q_frac_end * (ends[q] - starts[q])
      q_end_coord <- map_to_curve(q_angle_end, seqRadius[q] + ribbonGap[q], q_ref)

      # Generate query sequence alignment segment (curved path)
      q_angles <- seq(q_angle_start, q_angle_end, length.out = 50)
      q_coords <- do.call(rbind, lapply(q_angles, function(angle) {
        map_to_curve(angle, seqRadius[q] + ribbonGap[q], q_ref)
      }))
      segQ <- data.frame(x = q_coords[,1], y = q_coords[,2])

      # Correct subject sequence ribbon start/end coordinates
      s_ref <- seq_refs[[s]]
      s_frac_start <- if (orientation[s] == 1) (row$sstart - 1)/lens[s] else 1 - (row$sstart - 1)/lens[s]
      s_angle_start <- starts[s] + s_frac_start * (ends[s] - starts[s])
      s_start_coord <- map_to_curve(s_angle_start, seqRadius[s] + ribbonGap[s], s_ref)

      s_frac_end <- if (orientation[s] == 1) (row$send - 1)/lens[s] else 1 - (row$send - 1)/lens[s]
      s_angle_end <- starts[s] + s_frac_end * (ends[s] - starts[s])
      s_end_coord <- map_to_curve(s_angle_end, seqRadius[s] + ribbonGap[s], s_ref)

      # Generate subject sequence alignment segment (curved path)
      s_angles <- seq(s_angle_start, s_angle_end, length.out = 50)
      s_coords <- do.call(rbind, lapply(s_angles, function(angle) {
        map_to_curve(angle, seqRadius[s] + ribbonGap[s], s_ref)
      }))
      segS <- data.frame(x = s_coords[,1], y = s_coords[,2])

      # Determine Bézier curve control points (fix ribbon_ctrl_point parameter handling)
      if (!is.null(ribbon_ctrl_point)) {
        # Handle list format: each element corresponds to two control points (c1, c2) for one ribbon
        if (is.list(ribbon_ctrl_point)) {
          # Ensure list index does not exceed bounds
          cp_idx <- ifelse(i > length(ribbon_ctrl_point), length(ribbon_ctrl_point), i)
          cp <- ribbon_ctrl_point[[cp_idx]]
          # Ensure each control point element contains two points
          if (length(cp) >= 2) {
            c1 <- cp[[1]]  # Start curve control point
            c2 <- cp[[2]]  # End curve control point
          } else {
            # If insufficient elements, use the first point as default
            c1 <- c2 <- if (length(cp) == 1) cp[[1]] else c(0, 0)
          }
        } else {
          # Handle vector format: length 2 for single point, length 4 for two points (c1x,c1y,c2x,c2y)
          if (length(ribbon_ctrl_point) == 2) {
            # Single control point applied to all curves
            c1 <- c2 <- ribbon_ctrl_point
          } else if (length(ribbon_ctrl_point) == 4) {
            # Specify start and end control points separately
            c1 <- ribbon_ctrl_point[1:2]
            c2 <- ribbon_ctrl_point[3:4]
          } else {
            # Use default (center) for invalid lengths
            warning("ribbon_ctrl_point vector must have length 2 or 4; default value used")
            c1 <- c2 <- c(0, 0)
          }
        }
      } else {
        # Automatically calculate suitable control points based on curvature (original logic retained)
        mid_angle_q <- (q_angle_start + q_angle_end) / 2
        mid_angle_s <- (s_angle_start + s_angle_end) / 2
        mid_point_q <- map_to_curve(mid_angle_q, seqRadius[q] + ribbonGap[q] * 0.5, q_ref)
        mid_point_s <- map_to_curve(mid_angle_s, seqRadius[s] + ribbonGap[s] * 0.5, s_ref)
        c1 <- colMeans(rbind(mid_point_q, mid_point_s))
        c2 <- c1
      }

      # Generate Bézier curves (using fixed c1 and c2)
      b1 <- bezier_pts(
        as.numeric(segQ[1, ]),  # Start of start curve (query sequence start)
        as.numeric(segS[1, ]),  # End of start curve (subject sequence start)
        c1, c1,  # Use start control point c1
        n = 50
      )
      b2 <- bezier_pts(
        as.numeric(segQ[nrow(segQ), ]),  # Start of end curve (query sequence end)
        as.numeric(segS[nrow(segS), ]),  # End of end curve (subject sequence end)
        c2, c2,  # Use end control point c2
        n = 50
      )

      # Construct closed polygon
      poly <- rbind(
        segQ,  # 1. Query sequence segment: query start → query end (forward)
        b2,  # 2. End Bézier curve: query end → subject end (forward)
        segS[nrow(segS):1, ],  # 3. Subject sequence segment: subject end → subject start (reverse, connects to end of b2)
        b1[nrow(b1):1, ]  # 4. Start Bézier curve: subject start → query start (reverse, closes back to start)
      )

      # Set fill color
      if (ribbon_color_scheme == "pident") {
        poly$pident <- row$pident
      } else {
        poly$fill <- switch(ribbon_color_scheme,
                            single = singleCol,
                            query = queryCols[q])
      }
      poly$group <- i  # Group identifier
      ribbons[[length(ribbons) + 1]] <- poly
      cntValid <- cntValid + 1
    }

    if (debug) {
      cat("Valid ribbons:", cntValid, "Invalid ribbons:", cntInvalid, "\n")
    }
    if (cntValid > 0) {
      allRibbon <- do.call(rbind, ribbons)
    } else {
      allRibbon <- NULL
      warning("No valid ribbons to plot")
    }
  }


  # 11. Process gene annotation arrows — key modification for gene offset direction
  gene_polys <- data.frame()
  if (!is.null(gene_data) && nrow(gene_data) > 0) {
    valid_genes <- gene_data[gene_data$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      for (i in seq_len(nrow(valid_genes))) {
        gene <- valid_genes[i, ]
        sid <- gene$seq_id
        strand <- gene$strand
        anno <- gene$anno

        # Validate width
        width <- geneWidth[[sid]][ strand ]
        if (!is.numeric(width) || width <= 0) width <- 0.1

        # Calculate start/end angles (unchanged)
        seq_len <- lens[sid]
        sp <- min(gene$start, gene$end); ep <- max(gene$start, gene$end)
        if (ep <= sp) next
        frac_sp <- if (orientation[sid]==1) sp/seq_len else 1-sp/seq_len
        frac_ep <- if (orientation[sid]==1) ep/seq_len else 1-ep/seq_len
        a_start <- starts[sid] + frac_sp*(ends[sid]-starts[sid])
        a_end <- starts[sid] + frac_ep*(ends[sid]-starts[sid])
        if (strand == "-") { tmp <- a_start; a_start <- a_end; a_end <- tmp }

        # Generate angle vector and width factor (unchanged)
        n_body <- 30; n_head <- 15
        body_ang <- seq(a_start, a_start + 0.6*(a_end - a_start), length.out = n_body)
        head_ang <- seq(tail(body_ang,1), a_end, length.out = n_head)
        angs <- c(body_ang, head_ang)
        total_pt <- length(angs)
        widths <- c(rep(1, n_body), seq(1, 0, length.out = n_head))

        # Reference radius calculation — core modification: reverse gene_offset direction based on strand
        # Positive strand (+): gene_offset > 0 means outward offset (seqRadius - offset)
        # Negative strand (-): gene_offset > 0 means inward offset (seqRadius + offset)
        if (strand == "+") {
          r0 <- seqRadius[sid] - geneGap[[sid]][strand]  # Positive strand: minus sign keeps outward
        } else {
          r0 <- seqRadius[sid] + geneGap[[sid]][strand]  # Negative strand: plus sign achieves inward
        }

        # Calculate inner and outer radii (unchanged)
        outer_r <- r0 + (width/2)*widths
        inner_r <- r0 - (width/2)*widths

        # Original polar coordinate table (unchanged)
        orig_ang <- c(angs, rev(angs))
        orig_rad <- c(outer_r, rev(inner_r))

        # Map each point (unchanged)
        ref <- seq_refs[[sid]]
        pts_list <- mapply(function(ang, rad) {
          map_to_curve(ang, rad, ref)
        }, orig_ang, orig_rad, SIMPLIFY = FALSE)

        mapped <- do.call(rbind, pts_list)

        # Construct polygon (unchanged)
        gene_poly <- data.frame(
          x = mapped[,1],
          y = mapped[,2],
          group = i,
          anno = anno,
          strand = strand,
          ord = seq_len(2 * total_pt)
        )

        gene_polys <- rbind(gene_polys, gene_poly)
      }
    }
  }


  # 12. Process gene labels to ensure synchronization with arrow offsets
  gene_arrows <- data.frame()
  if (!is.null(gene_data) && nrow(gene_data) > 0 && gene_label_show) {
    valid_genes <- gene_data[gene_data$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      gene_arrows <- do.call(rbind, lapply(seq_len(nrow(valid_genes)), function(i) {
        gene <- valid_genes[i, ]
        sid <- gene$seq_id
        strand <- gene$strand
        seq_len <- lens[sid]
        ref <- seq_refs[[sid]]
        orient <- orientation[sid]  # 1=forward, -1=reverse

        # 1. Calculate relative position of gene midpoint (unchanged)
        sp <- min(gene$start, gene$end); ep <- max(gene$start, gene$end)
        frac_mid <- (sp + ep) / (2 * seq_len)

        # Apply circumferential offset (unchanged)
        circum_ratio <- geneLabelCircumOffset[[sid]][strand]
        if (geneLabelCircumLimit[[sid]][strand]) {
          gene_length_ratio <- (ep - sp) / seq_len
          max_offset_ratio <- gene_length_ratio * 0.5
          circum_ratio <- pmin(max_offset_ratio, pmax(-max_offset_ratio, circum_ratio))
        }
        frac_mid <- frac_mid + circum_ratio
        frac_mid <- pmin(1, pmax(0, frac_mid))

        # Adjust midpoint for reverse sequences (unchanged)
        if (orient != 1) frac_mid <- 1 - frac_mid

        # 2. Find reference path index and tangent vector (dx,dy) (unchanged)
        ref_n <- length(ref$angles)
        idx <- round(frac_mid * (ref_n - 1)) + 1
        idx <- pmin(ref_n, pmax(1, idx))
        if (idx < ref_n) {
          dx <- ref$path$x[idx + 1] - ref$path$x[idx]
          dy <- ref$path$y[idx + 1] - ref$path$y[idx]
        } else {
          dx <- ref$path$x[idx] - ref$path$x[idx - 1]
          dy <- ref$path$y[idx] - ref$path$y[idx - 1]
        }
        dx <- dx * orient  # Correct tangent direction based on sequence orientation
        dy <- dy * orient

        # 3. Calculate arrow center point (core modification: use same reference radius logic as arrows)
        width <- geneWidth[[sid]][ strand ]

        # Consistent reference radius calculation with arrows
        if (strand == "+") {
          r0 <- seqRadius[sid] - geneGap[[sid]][strand]  # Positive strand: outward offset
        } else {
          r0 <- seqRadius[sid] + geneGap[[sid]][strand]  # Negative strand: inward offset
        }

        # Calculate arrow center radius (midpoint between inner and outer radii of arrow)
        center_r <- r0  # Arrow center is at the reference radius

        # Calculate center coordinates
        center_pt <- map_to_curve(angle = ref$angles[idx], radius = center_r, ref = ref)

        # 4. Calculate normal vector (ensure alignment with arrow offset direction)
        normal_x <- -dy  # Basic normal direction (perpendicular to tangent)
        normal_y <- dx
        nl <- sqrt(normal_x^2 + normal_y^2)
        if (nl > 0) {
          normal_x <- normal_x / nl  # Normalize
          normal_y <- normal_y / nl
        }

        # Adjust normal direction based on strand and sequence orientation
        direction_factor <- ifelse(strand == "+", 1, -1) * orient
        normal_x <- normal_x * direction_factor
        normal_y <- normal_y * direction_factor

        # 5. Calculate label position (synchronized with arrow offset)
        text_x <- center_pt[1] - normal_x * geneLabelRadialOffset[[sid]][strand]
        text_y <- center_pt[2] - normal_y * geneLabelRadialOffset[[sid]][strand]

        # 6. Calculate text angle and alignment (unchanged)
        base_angle <- atan2(dy, dx) * 180 / pi
        text_angle <- base_angle + 90 + geneLabelRotation[[sid]][strand]

        # Automatic alignment logic (adjust left/right alignment based on direction)
        if (strand == "+" && orient == 1) {
          hjust <- 1
        } else if (strand == "+" && orient != 1) {
          hjust <- 0
        } else if (strand == "-" && orient == 1) {
          hjust <- 0
        } else {
          hjust <- 1
        }

        # Ensure text direction is correct (not upside down)
        text_angle <- (text_angle + 360) %% 360
        if (text_angle > 90 && text_angle < 270) {
          text_angle <- text_angle + 180
          hjust <- 1 - hjust
        }
        text_angle <- text_angle %% 360

        data.frame(
          text = gene$anno,
          text_x = text_x,
          text_y = text_y,
          text_angle = text_angle,
          hjust = hjust,
          vjust = 0.5,
          seq_id = sid,
          group = i,
          stringsAsFactors = FALSE
        )
      }))
    }
  }


  # Define rotation function
  rotate_df <- function(df) {
    if (all(c("x","y") %in% names(df))) {
      x0 <- df$x; y0 <- df$y
      df$x <- x0 * cos(rot_rad) - y0 * sin(rot_rad)
      df$y <- x0 * sin(rot_rad) + y0 * cos(rot_rad)
    }
    if (all(c("x0","y0","x1","y1") %in% names(df))) {
      X <- df$x0; Y <- df$y0
      df$x0 <- X * cos(rot_rad) - Y * sin(rot_rad)
      df$y0 <- X * sin(rot_rad) + Y * cos(rot_rad)
      X1 <- df$x1; Y1 <- df$y1
      df$x1 <- X1 * cos(rot_rad) - Y1 * sin(rot_rad)
      df$y1 <- X1 * sin(rot_rad) + Y1 * cos(rot_rad)
    }
    if (all(c("label_x","label_y") %in% names(df))) {
      LX <- df$label_x; LY <- df$label_y
      df$label_x <- LX * cos(rot_rad) - LY * sin(rot_rad)
      df$label_y <- LX * sin(rot_rad) + LY * cos(rot_rad)
    }
    if (all(c("x_start","y_start","x_end","y_end") %in% names(df))) {
      Xs <- df$x_start; Ys <- df$y_start
      df$x_start <- Xs * cos(rot_rad) - Ys * sin(rot_rad)
      df$y_start <- Xs * sin(rot_rad) + Ys * cos(rot_rad)
      Xe <- df$x_end; Ye <- df$y_end
      df$x_end <- Xe * cos(rot_rad) - Ye * sin(rot_rad)
      df$y_end <- Xe * sin(rot_rad) + Ye * cos(rot_rad)
    }
    if (all(c("text_x","text_y") %in% names(df))) {
      TX <- df$text_x; TY <- df$text_y
      df$text_x <- TX * cos(rot_rad) - TY * sin(rot_rad)
      df$text_y <- TX * sin(rot_rad) + TY * cos(rot_rad)
      df$text_angle <- df$text_angle + rotation
    }
    df
  }

  # apply rotation
  seqArcs <- lapply(seqArcs, rotate_df)
  innerArcs <- lapply(innerArcs, rotate_df)
  axisLines <- rotate_df(axisLines)
  axisTicks <- rotate_df(axisTicks)
  if (!is.null(allRibbon)) {
    allRibbon <- rotate_df(allRibbon)
  }
  if (nrow(gene_arrows) > 0) {
    gene_arrows <- rotate_df(gene_arrows)
  }
  if (nrow(gene_polys) > 0) {
    gene_polys <- rotate_df(gene_polys)
    gene_polys <- gene_polys[with(gene_polys, order(group, ord)), ]
  }

  extremes <- get_plot_extremes(
    allRibbon = allRibbon,
    seqArcs = seqArcs,
    axisLines = axisLines,
    axisTicks = axisTicks,
    gene_polys = gene_polys,
    gene_arrows = NULL,
    show_axis = show_axis
  )

  chordPlot <- chordPlotFunc(allRibbon,ribbon_alpha,ribbon_color_scheme,ribbon_colors,show_legend,gene_polys,gene_pal,gene_color_scheme,final_gene_order,seqArcs,gene_arrows,gene_label_show,gene_label_size,show_axis,axisLines,axisTicks,axisLabelOrientation,seq_colors,seq_labels,seqs,extremes,panel_margin,title)

  return(chordPlot)
}
