# ggchord.R - constructor
# Entry point of the layered ggplot2 API for chord diagrams
# v0.4.0: ggchord() keeps only data arguments and global style parameters,
# layout parameters moved to the geom_* layers, and coordinate computation is deferred to print time

# Global variable declarations (to avoid R CMD check NOTEs)
globalVariables(c(
  "x", "y", "group", "pident", "fill", "strand", "anno", "seq_id",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "label", "label_x", "label_y", "size"
))

#' ggchord: layered multi-sequence alignment chord diagrams for ggplot2
#'
#' ggchord visualizes multi-sequence BLAST alignment results using ggplot2's layered grammar.
#' The \code{ggchord()} constructor handles data validation and global settings;
#' the \code{geom_*} layers are stacked as needed, each responsible for its own layout parameters and visual rendering.
#' The layout is computed lazily at \code{print()} time.
#'
#' @param seq_data data.frame/tibble, required. Basic sequence information
#' @param ribbon_data data.frame/tibble, optional. BLAST alignment results
#' @param gene_data data.frame/tibble, optional. Gene annotation data
#' @param title Character. Main title of the plot, default NULL
#' @param rotation Numeric. Global rotation angle (degrees), default 45
#' @param panel_margin Optional numeric/list. Panel margin, default 0
#' @param show_legend Logical. Whether to show legends, default TRUE
#' @param debug Logical. Whether to output debug information, default FALSE
#'
#' @return A ggchord object (inherits from ggplot) to which geom_* layers can be added with +
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' data(gene_data_example)
#'
#' p <- ggchord(
#'   seq_data = seq_data_example,
#'   ribbon_data = ribbon_data_example,
#'   gene_data = gene_data_example
#' ) +
#'   geom_seq() +
#'   geom_ribbon() +
#'   geom_gene() +
#'   geom_axis()
#' print(p)
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import grDevices
#' @import grid
ggchord <- function(
    seq_data,
    ribbon_data = NULL,
    gene_data = NULL,
    title = NULL,
    rotation = 45,
    panel_margin = 0,
    show_legend = TRUE,
    debug = FALSE
) {
  # Clear parameters from the previous run
  clear_chord_env()

  # ====================================================================
  # 1. Validate data
  # ====================================================================
  required_seq_cols <- c("seq_id", "length")
  if (!all(required_seq_cols %in% colnames(seq_data))) {
    stop("seq_data must contain the following columns: ", paste(required_seq_cols, collapse = ", "))
  }
  if (any(seq_data$length <= 0)) {
    stop("The 'length' values in seq_data must be positive")
  }
  if (anyDuplicated(seq_data$seq_id)) {
    stop("The 'seq_id' values in seq_data must be unique")
  }

  if (!is.null(ribbon_data)) {
    required_ribbon_cols <- c("qaccver", "saccver", "length", "pident",
                              "qstart", "qend", "sstart", "send")
    if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
      stop("ribbon_data must contain the following columns: ",
           paste(required_ribbon_cols, collapse = ", "))
    }
    if (nrow(ribbon_data) == 0) warning("No valid alignment data in ribbon_data")
    if (debug) cat("Number of alignment data rows: ", nrow(ribbon_data), "\n")
  }

  if (!is.null(gene_data)) {
    required_gene_cols <- c("seq_id", "start", "end", "strand", "anno")
    if (!all(required_gene_cols %in% colnames(gene_data))) {
      stop("gene_data must contain the following columns: ",
           paste(required_gene_cols, collapse = ", "))
    }
    if (nrow(gene_data) == 0) warning("No valid gene annotation data in gene_data")
    if (any(!gene_data$strand %in% c("+", "-"))) {
      stop("The 'strand' values in gene_data can only be '+' or '-'")
    }
    if (debug) cat("Number of gene annotation rows: ", nrow(gene_data), "\n")
  }

  # ====================================================================
  # 2. Store raw data and global parameters
  # ====================================================================
  set_chord_data(seq_data, ribbon_data, gene_data)

  set_global_params(list(
    rotation     = rotation,
    panel_margin = panel_margin,
    show_legend  = show_legend,
    debug        = debug
  ))

  # ====================================================================
  # 3. Build the base ggplot object (the coord is replaced with the correct range at print time)
  # ====================================================================
  margin_vals <- process_panel_margin(panel_margin)

  p <- ggplot() +
    coord_chord() +
    labs(title = title) +
    theme(
      plot.title       = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.margin      = margin(t = margin_vals$t, r = margin_vals$r,
                                b = margin_vals$b, l = margin_vals$l),
      legend.background = element_blank(),
      legend.box.spacing  = unit(10, "mm"),
      legend.spacing      = unit(5, "mm"),
      legend.text         = element_text(size = 8),
      legend.title        = element_text(size = 10, face = "bold"),
      axis.title          = element_blank(),
      axis.line           = element_blank(),
      axis.ticks          = element_blank(),
      axis.text           = element_blank(),
      panel.background    = element_blank(),
      legend.position     = if (isTRUE(show_legend)) "right" else "none"
    )

  class(p) <- c("ggchord", class(p))
  p
}

# ====================================================================
# +.ggchord method
# ====================================================================

#' Combine a ggchord plot with ggplot2 objects
#'
#' Supports stacking ggplot2 layers, lists of layers, scales, and themes
#' onto a ggchord plot using the \code{+} operator.
#'
#' @param e1 A ggchord object
#' @param e2 A ggplot2 layer, a list of layers, a scale, or a theme
#' @return A ggchord object
#' @export
`+.ggchord` <- function(e1, e2) {
  if (is.list(e2) && !inherits(e2, "ggplot") && !inherits(e2, "LayerInstance")) {
    p <- e1
    for (elem in e2) {
      p <- p + elem
    }
  } else {
    p <- NextMethod()
  }
  class(p) <- unique(c("ggchord", class(p)))
  p
}

# ====================================================================
# print.ggchord: lazily compute the layout, inject data into layers, and render
# ====================================================================

#' @export
print.ggchord <- function(x, ...) {
  # Step 1: collect all parameters
  global <- get_global_params()
  data_list <- get_chord_data()

  seq_params    <- get_seq_params()
  ribbon_params <- get_ribbon_params()
  gene_params   <- get_gene_params()
  axis_params   <- get_axis_params()

  seq_layer_requested <- !is.null(seq_params)
  if (is.null(seq_params)) seq_params <- list()

  # --- Process sequences ---
  seq_data <- data_list$seq_data
  seqs     <- seq_data$seq_id
  lens     <- setNames(seq_data$length, seqs)

  if (!is.null(seq_params$seq_order)) {
    if (!all(seq_params$seq_order %in% seqs)) {
      stop("seq_order contains unknown sequence IDs")
    }
    seqs <- seq_params$seq_order
    lens <- lens[seqs]
  }
  n <- length(seqs)

  seq_labels    <- process_sequence_param(seq_params$seq_labels, seqs,
                                          "seq_labels", default_value = seqs)
  seqRadius     <- process_sequence_param(seq_params$seq_radius, seqs,
                                          "seq_radius", 1.0)
  orientation   <- process_sequence_param(seq_params$seq_orientation, seqs,
                                          "seq_orientation", 1)
  seq_gap       <- process_sequence_param(seq_params$seq_gap, seqs,
                                          "seq_gap", 0.03)
  seq_curvature <- process_sequence_param(seq_params$seq_curvature, seqs,
                                          "seq_curvature", 1.0)

  if (any(seqRadius <= 0)) stop("seq_radius must be positive")
  if (any(!orientation %in% c(-1, 1))) {
    stop("seq_orientation can only be 1 or -1")
  }
  if (any(seq_gap < 0 | seq_gap >= 0.5)) {
    stop("seq_gap must be in the [0, 0.5) range")
  }

  if (!is.null(seq_params$seq_colors)) {
    seq_colors <- process_sequence_param(seq_params$seq_colors, seqs, "seq_colors")
  } else {
    pal <- if (n <= 9) {
             # brewer.pal(n<3) errors/warns; take the first n colors
             brewer.pal(max(n, 3), "Set1")[seq_len(n)]
           } else {
             colorRampPalette(brewer.pal(9, "Set1"))(n)
           }
    seq_colors <- setNames(pal, seqs)
  }

  # --- Process ribbons ---
  ribbonGap  <- process_sequence_param(ribbon_params$ribbon_gap %||% 0.15,
                                       seqs, "ribbon_gap", 0.15)
  ribbon_color_scheme <- ribbon_params$ribbon_color_scheme %||% "pident"
  ribbon_alpha    <- ribbon_params$ribbon_alpha %||% 0.35
  ribbon_ctrl_pt  <- ribbon_params$ribbon_ctrl_point %||% c(0, 0)

  ribbon_colors <- ribbon_params$ribbon_colors
  if (!ribbon_color_scheme %in% c("pident", "query", "single")) {
    stop("ribbon_color_scheme must be 'pident', 'query', or 'single'")
  }
  if (!is.numeric(ribbon_alpha) || length(ribbon_alpha) != 1 ||
      ribbon_alpha < 0 || ribbon_alpha > 1) {
    stop("ribbon_alpha must be in the [0, 1] range")
  }

  # ribbon_colors validation only runs when ribbon_data is actually present
  has_ribbon_data <- !is.null(data_list$ribbon_data) &&
                     nrow(data_list$ribbon_data) > 0

  if (has_ribbon_data) {
    if (is.null(ribbon_colors)) {
      ribbon_colors <- switch(ribbon_color_scheme,
        single = "steelblue",
        query  = {
          mix <- 0.5
          sapply(seq_colors, function(col) {
            cols <- col2rgb(col)
            light_cols <- cols + (255 - cols) * mix
            rgb(light_cols[1,], light_cols[2,], light_cols[3,],
                maxColorValue = 255)
          })
        },
        pident = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF",
                    "#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF",
                    "#B4DE2CFF","#FDE725FF"))
    }
    if (ribbon_color_scheme == "query") {
      ribbon_colors <- process_sequence_param(ribbon_colors, seqs,
                                              "ribbon_colors")
    } else if (ribbon_color_scheme == "pident" && length(ribbon_colors) < 2) {
      stop("The 'pident' scheme requires at least two ribbon_colors")
    } else if (ribbon_color_scheme == "single" && length(ribbon_colors) < 1) {
      stop("The 'single' scheme requires at least one ribbon_colors")
    }
  }

  # --- Process genes ---
  gene_off  <- gene_params$gene_offset %||% 0.1
  gene_w    <- gene_params$gene_width %||% 0.05
  gene_cs   <- gene_params$gene_color_scheme %||% "strand"
  gene_cols <- gene_params$gene_colors
  gene_ord  <- gene_params$gene_order
  gene_ls   <- gene_params$show_label_override %||%
    gene_params$gene_label_show %||% FALSE
  gene_lsz  <- gene_params$label_size_override %||%
    gene_params$gene_label_size %||% 2.5
  gene_lr   <- gene_params$gene_label_rotation %||% 0
  gene_lro  <- gene_params$gene_label_radial_offset %||% 0
  gene_lco  <- gene_params$gene_label_circum_offset %||% 0
  gene_lcl  <- if (is.null(gene_params$gene_label_circum_limit)) TRUE
               else gene_params$gene_label_circum_limit

  if (!gene_cs %in% c("strand", "manual")) {
    stop("gene_color_scheme must be 'strand' or 'manual'")
  }

  geneGap    <- process_gene_param(gene_off, seqs, "gene_offset", 0.1, FALSE)
  geneWidth  <- process_gene_param(gene_w, seqs, "gene_width", 0.05, FALSE)
  geneLabelRadialOffset <- process_gene_param(gene_lro, seqs,
                                              "gene_label_radial_offset", 0, FALSE)
  geneLabelCircumOffset <- process_gene_param(gene_lco, seqs,
                                              "gene_label_circum_offset", 0, FALSE)
  geneLabelCircumLimit  <- process_gene_param(gene_lcl, seqs,
                                              "gene_label_circum_limit", TRUE, TRUE)
  geneLabelRotation     <- process_gene_param(gene_lr, seqs,
                                              "gene_label_rotation", 0, FALSE)

  # --- Process axes ---
  show_axis  <- axis_params$show_axis %||% TRUE
  axisGap    <- process_sequence_param(axis_params$axis_gap %||% 0.04,
                                       seqs, "axis_gap", 0.04)
  axisMaj    <- process_sequence_param(axis_params$axis_tick_major_number %||% 5,
                                       seqs, "axis_tick_major_number", 5)
  axisMajLen <- process_sequence_param(axis_params$axis_tick_major_length %||% 0.02,
                                       seqs, "axis_tick_major_length", 0.02)
  axisMin    <- process_sequence_param(axis_params$axis_tick_minor_number %||% 4,
                                       seqs, "axis_tick_minor_number", 4)
  axisMinLen <- process_sequence_param(axis_params$axis_tick_minor_length %||% 0.01,
                                       seqs, "axis_tick_minor_length", 0.01)
  labelSize  <- process_sequence_param(axis_params$axis_label_size %||% 3,
                                       seqs, "axis_label_size", 3)
  labelOffset <- process_sequence_param(axis_params$axis_label_offset %||% 1.5,
                                        seqs, "axis_label_offset", 1.5)
  axisLabelOri <- process_axis_orientation(
    axis_params$axis_label_orientation %||% "horizontal", seqs
  )

  # ====================================================================
  # Step 2: compute the layout
  # ====================================================================
  layout <- compute_chord_layout(
    seqs = seqs, lens = lens, seq_labels = seq_labels,
    seq_colors = seq_colors, seqRadius = seqRadius,
    seq_curvature = seq_curvature, orientation = orientation,
    seq_gap = seq_gap,
    ribbon_data = data_list$ribbon_data, ribbonGap = ribbonGap,
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors = ribbon_colors, ribbon_alpha = ribbon_alpha,
    ribbon_ctrl_point = ribbon_ctrl_pt,
    gene_data = data_list$gene_data,
    geneGap = geneGap, geneWidth = geneWidth,
    geneLabelRadialOffset = geneLabelRadialOffset,
    geneLabelCircumOffset = geneLabelCircumOffset,
    geneLabelCircumLimit = geneLabelCircumLimit,
    geneLabelRotation = geneLabelRotation,
    gene_label_show = gene_ls, gene_label_size = gene_lsz,
    gene_color_scheme = gene_cs, gene_colors = gene_cols,
    gene_order = gene_ord,
    axisGap = axisGap, axisMaj = axisMaj, axisMajLen = axisMajLen,
    axisMin = axisMin, axisMinLen = axisMinLen,
    labelSize = labelSize, labelOffset = labelOffset,
    axisLabelOrientation = axisLabelOri,
    show_axis = show_axis,
    rotation = global$rotation, debug = global$debug
  )

  set_chord_layout(layout)

  # ====================================================================
  # Step 3: inject data into layers
  # ====================================================================
  # Classify layers by geom type
  seq_indices    <- integer(0)
  ribbon_indices <- integer(0)
  gene_poly_indices <- integer(0)
  gene_text_indices <- integer(0)
  axis_line_indices <- integer(0)
  axis_seg_indices  <- integer(0)
  axis_text_indices <- integer(0)

  seq_path_assigned <- FALSE
  for (i in seq_along(x$layers)) {
    lyr <- x$layers[[i]]
    gname <- class(lyr$geom)[1]
    data_names <- names(lyr$data)

    if (gname == "GeomPath") {
      # Sequence and axis paths share a geom; their placeholder columns do not.
      if ("seq_id" %in% data_names) {
        # geom_seq(), when present, is the first such path.  A plot containing
        # only geom_axis() must retain its axis path as an axis path.
        if (!seq_path_assigned && seq_layer_requested) {
          seq_indices <- c(seq_indices, i)
          seq_path_assigned <- TRUE
        } else {
          axis_line_indices <- c(axis_line_indices, i)
        }
      }
    } else if (gname %in% c("GeomPolygon", "NewGeomPolygon", "GeomChordPolygon")) {
      # The gene placeholder is distinguished from the ribbon placeholder
      # by its strand column (both use GeomPolygon).
      if ("strand" %in% data_names) {
        gene_poly_indices <- c(gene_poly_indices, i)
      } else {
        ribbon_indices <- c(ribbon_indices, i)
      }
    } else if (gname == "GeomSegment") {
      axis_seg_indices <- c(axis_seg_indices, i)
    } else if (gname == "GeomText") {
      # Gene labels use text_x/text_y; axis labels use label_x/label_y.
      if ("text_x" %in% data_names) {
        gene_text_indices <- c(gene_text_indices, i)
      } else {
        axis_text_indices <- c(axis_text_indices, i)
      }
    }
  }

  # Inject sequence arc data
  if (length(seq_indices) > 0 && length(layout$seq_arcs) > 0) {
    arc_df <- do.call(rbind, layout$seq_arcs)
    for (idx in seq_indices) x$layers[[idx]]$data <- arc_df
  }

  # Inject ribbon data
  if (length(ribbon_indices) > 0 && !is.null(layout$ribbon_polys)) {
    for (idx in ribbon_indices) x$layers[[idx]]$data <- layout$ribbon_polys
  }

  # Inject gene data
  if (length(gene_poly_indices) > 0 && nrow(layout$gene_polys) > 0) {
    for (idx in gene_poly_indices) x$layers[[idx]]$data <- layout$gene_polys
  }

  # Inject gene labels
  if (length(gene_text_indices) > 0 && nrow(layout$gene_labels) > 0) {
    for (idx in gene_text_indices) x$layers[[idx]]$data <- layout$gene_labels
  }

  # Inject axis line data
  if (length(axis_line_indices) > 0 && nrow(layout$axis_lines) > 0) {
    for (idx in axis_line_indices) x$layers[[idx]]$data <- layout$axis_lines
  }

  # Inject tick data
  if (length(axis_seg_indices) > 0 && nrow(layout$axis_ticks) > 0) {
    for (idx in axis_seg_indices) x$layers[[idx]]$data <- layout$axis_ticks
  }

  # Inject tick label data
  label_data <- if (nrow(layout$axis_ticks) > 0) {
    subset(layout$axis_ticks, !is.na(label))
  } else {
    layout$axis_ticks
  }
  if (length(axis_text_indices) > 0 && nrow(label_data) > 0) {
    for (idx in axis_text_indices) x$layers[[idx]]$data <- label_data
  }

  # ====================================================================
  # Step 4: dynamically set up scales
  # Note: temporarily remove the ggchord class, otherwise +.ggchord conflicts with the +.gg method
  # ====================================================================
  cls <- class(x)
  class(x) <- setdiff(cls, "ggchord")

  # Sequence color scale (only added when a sequence arc layer exists, to avoid warnings with no colour data)
  if (length(seq_indices) > 0) {
    x <- x + scale_color_manual(
      name   = "Seq ID",
      values = layout$seq_colors,
      labels = layout$seq_labels,
      breaks = layout$seqs
    )
  }

  # Ribbon fill scale
  # Note: when the gene layer and the ribbon share the fill aesthetic, the
  # ribbon layer's fill mapping is renamed to the internal aesthetic
  # "fill_ggnewscale_1" below so the gene scale does not overwrite the ribbon
  # scale.
  ribbon_fill_scale <- NULL
  if (!is.null(layout$ribbon_polys)) {
    if (layout$ribbon_color_scheme == "pident") {
      ribbon_fill_scale <- scale_fill_stepsn(
        name    = "Identity(%)",
        colours = layout$ribbon_colors,
        limits  = c(0, 100),
        breaks  = c(0, 50, 80, 90, 95, 100),
        guide   = guide_colorbar(
          theme = theme(legend.title.position = "top",
                        legend.key.height = unit(1, "null")),
          order = 1, position = "left"
        )
      )
    } else {
      ribbon_fill_scale <- scale_fill_identity()
    }
  }

  # Gene fill scale
  gene_fill_scale <- NULL
  if (nrow(layout$gene_polys) > 0) {
    if (layout$gene_color_scheme == "strand") {
      gene_fill_scale <- scale_fill_manual(
        name   = "Strand",
        breaks = c("+", "-"),
        values = layout$gene_pal
      )
    } else {
      gene_fill_scale <- scale_fill_manual(
        name   = "Gene Annotation",
        breaks = layout$final_gene_order,
        values = layout$gene_pal
      )
    }
  }

  # Merge the two fill scales:
  # - ribbon only: add the ribbon scale directly (aesthetic = fill)
  # - gene only: add the gene scale directly (aesthetic = fill)
  # - both: ggplot2 allows only one scale per aesthetic, so the ribbon layers'
  #   fill mapping is renamed to the internal aesthetic "fill_ggnewscale_1" and
  #   the ribbon scale is attached to it; the gene scale keeps the plain "fill"
  #   aesthetic. The two scales therefore do not overwrite each other.
  if (!is.null(ribbon_fill_scale) && !is.null(gene_fill_scale)) {
    gns_aes <- "fill_ggnewscale_1"
    for (r_idx in ribbon_indices) {
      lyr <- x$layers[[r_idx]]
      if (!is.null(lyr$mapping)) {
        mp_names <- names(lyr$mapping)
        mp_names[mp_names == "fill"] <- gns_aes
        names(lyr$mapping) <- mp_names
      }
    }
    s <- ribbon_fill_scale
    s$aesthetics <- gns_aes
    # guide_colorbar's available_aes only recognizes fill by default; after the
    # aesthetic rename it must be synchronized, otherwise the colorbar is dropped.
    if (inherits(s$guide, "Guide")) {
      s$guide$available_aes <- gsub("^fill$", gns_aes, s$guide$available_aes)
      if (!is.null(s$guide$params$override.aes)) {
        names(s$guide$params$override.aes) <-
          gsub("^fill$", gns_aes, names(s$guide$params$override.aes))
      }
    }
    x <- x + s
    # The gene layer uses the fill aesthetic, so add the gene scale directly (no conflict)
    x <- x + gene_fill_scale
  } else if (!is.null(ribbon_fill_scale)) {
    x <- x + ribbon_fill_scale
  } else if (!is.null(gene_fill_scale)) {
    x <- x + gene_fill_scale
  }

  # Ribbon alpha is a preset value; use an identity scale so it renders as specified
  if (!is.null(layout$ribbon_polys)) {
    x <- x + scale_alpha_identity()
  }

  # Axis text size scale
  x <- x + scale_size_identity()

  # Restore the ggchord class
  class(x) <- unique(c("ggchord", class(x)))

  # ====================================================================
  # Step 5: update the coord range
  # ====================================================================
  ext <- layout$extremes
  pad <- 0.05 * max(ext$x_max - ext$x_min, ext$y_max - ext$y_min, 1)
  x$coordinates <- coord_fixed(
    ratio = 1,
    xlim  = c(ext$x_min - pad, ext$x_max + pad),
    ylim  = c(ext$y_min - pad, ext$y_max + pad),
    clip  = "off"
  )

  # ====================================================================
  # Step 6: render
  # ====================================================================
  cls <- class(x)
  class(x) <- setdiff(cls, c("ggchord"))
  print(x)
  invisible()
}

# ====================================================================
# ggplot_build.ggchord
# ====================================================================

#' @export
ggplot_build.ggchord <- function(plot, ...) {
  # The layout is computed lazily by print.ggchord(). A direct ggplot_build()
  # call on a fresh plot therefore has no layout yet; warn in that case.
  try_layout <- tryCatch(get_chord_layout(), error = function(e) NULL)
  if (is.null(try_layout)) {
    calls <- sys.calls()
    internal_build <- any(vapply(calls, function(cl) {
      grepl("ggplot_add|new_aes", deparse(cl)[1])
    }, logical(1)))
    if (!internal_build) {
      warning("Layout has not been computed; please render via print().")
    }
  }
  NextMethod()
}
