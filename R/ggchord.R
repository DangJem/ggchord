# ggchord.R - constructor
# Entry point of the layered ggplot2 API for chord diagrams
# v0.6.0: plot objects are self-contained. Data and parameters are stored on
# the plot object itself, and the layout is computed at build time so that
# print(), ggsave(), ggplot_build() and other ggplot2 workflows all work.

# Global variable declarations (to avoid R CMD check NOTEs)
globalVariables(c(
  "x", "y", "group", "pident", "fill", "strand", "anno", "seq_id",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "label", "label_x", "label_y", "size"
))

#' ggchord: layered multi-sequence alignment chord diagrams for ggplot2
#'
#' ggchord visualizes multi-sequence alignment results using ggplot2's layered grammar.
#' The \code{ggchord()} constructor handles data validation and global settings;
#' the \code{geom_*} layers are stacked as needed, each responsible for its own layout parameters and visual rendering.
#' The layout is computed lazily when the plot is built (e.g. via \code{print()},
#' \code{ggsave()}, or \code{ggplot_build()}).
#'
#' @param seq_data data.frame/tibble, required. Basic sequence information
#' @param ribbon_data data.frame/tibble, optional. Alignment results
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

  seq_lens <- setNames(seq_data$length, seq_data$seq_id)

  if (!is.null(ribbon_data)) {
    required_ribbon_cols <- c("qaccver", "saccver", "length", "pident",
                              "qstart", "qend", "sstart", "send")
    if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
      stop("ribbon_data must contain the following columns: ",
           paste(required_ribbon_cols, collapse = ", "))
    }
    if (nrow(ribbon_data) == 0) warning("No valid alignment data in ribbon_data")
    if (nrow(ribbon_data) > 0) {
      unknown <- setdiff(unique(c(ribbon_data$qaccver, ribbon_data$saccver)),
                         seq_data$seq_id)
      if (length(unknown) > 0) {
        warning("ribbon_data contains sequence IDs not present in seq_data: ",
                paste(unknown, collapse = ", "))
      }
      if (any(ribbon_data$qstart > ribbon_data$qend |
              ribbon_data$sstart > ribbon_data$send, na.rm = TRUE)) {
        warning("ribbon_data contains rows where start > end; these may render abnormally")
      }
      both_known <- ribbon_data$qaccver %in% seq_data$seq_id &
                    ribbon_data$saccver %in% seq_data$seq_id
      if (any(both_known)) {
        out_of_range <- (ribbon_data$qstart[both_known] < 1 |
                         ribbon_data$qend[both_known] > seq_lens[ribbon_data$qaccver[both_known]] |
                         ribbon_data$sstart[both_known] < 1 |
                         ribbon_data$send[both_known] > seq_lens[ribbon_data$saccver[both_known]])
        if (any(out_of_range, na.rm = TRUE)) {
          warning("ribbon_data contains alignment positions outside the sequence length")
        }
      }
      if (any(ribbon_data$pident < 0 | ribbon_data$pident > 100, na.rm = TRUE)) {
        warning("ribbon_data contains pident values outside [0, 100]")
      }
    }
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
    if (nrow(gene_data) > 0) {
      unknown <- setdiff(unique(gene_data$seq_id), seq_data$seq_id)
      if (length(unknown) > 0) {
        warning("gene_data contains sequence IDs not present in seq_data: ",
                paste(unknown, collapse = ", "))
      }
      known <- gene_data$seq_id %in% seq_data$seq_id
      if (any(known)) {
        if (any(gene_data$start[known] > gene_data$end[known])) {
          warning("gene_data contains rows where start > end")
        }
        out_of_range <- gene_data$start[known] < 1 |
                        gene_data$end[known] > seq_lens[gene_data$seq_id[known]]
        if (any(out_of_range, na.rm = TRUE)) {
          warning("gene_data contains positions outside the sequence length")
        }
      }
    }
    if (debug) cat("Number of gene annotation rows: ", nrow(gene_data), "\n")
  }

  # ====================================================================
  # 2. Build the base ggplot object and store data + global parameters
  #    on the plot itself so the object is fully self-contained.
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

  p$ggchord <- list(
    data   = list(seq_data = seq_data, ribbon_data = ribbon_data,
                  gene_data = gene_data),
    global = list(rotation = rotation, panel_margin = panel_margin,
                  show_legend = show_legend, debug = debug)
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
# ggplot_build.ggchord: compute the layout, inject data into cloned
# layers, add scales, and set the coordinate system.  Everything is
# driven by the plot object itself, so the plot can be printed, saved
# with ggsave(), or built with ggplot_build() any number of times and
# in any order, without cross-talk between plots.
# ====================================================================

#' @export
ggplot_build.ggchord <- function(plot, ...) {
  # Step 1: collect data and parameters from the plot object
  chord <- plot$ggchord
  if (is.null(chord)) {
    stop("Not a valid ggchord object: no data stored on the plot. ",
         "Please build the plot with ggchord().")
  }
  data_list <- chord$data
  global    <- chord$global

  # Work on a copy so the user's plot object is never modified: the scales are
  # cloned because scales are added in place below.
  plot$scales <- plot$scales$clone()

  seq_params    <- list()
  ribbon_params <- list()
  gene_params   <- list()
  axis_params   <- list()
  seq_layer_requested <- FALSE

  for (i in seq_along(plot$layers)) {
    pp <- plot$layers[[i]]$ggchord_params
    if (is.null(pp)) next
    switch(pp$type,
      seq    = { seq_params <- pp; seq_layer_requested <- TRUE },
      ribbon = ribbon_params <- pp,
      gene   = gene_params <- pp,
      axis   = axis_params <- pp
    )
  }

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
    pal <- chord_default_palette(n)
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

  # Cache the layout so the get_chord_layout() accessor can inspect it.
  set_chord_layout(layout)

  # ====================================================================
  # Step 3: classify layers and inject data into CLONED layers
  # (cloning keeps the user's plot object untouched)
  # ====================================================================
  seq_indices    <- integer(0)
  ribbon_indices <- integer(0)
  gene_poly_indices <- integer(0)
  gene_text_indices <- integer(0)
  axis_line_indices <- integer(0)
  axis_seg_indices  <- integer(0)
  axis_text_indices <- integer(0)

  seq_path_assigned <- FALSE
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
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
      # by its strand column (both use polygon geoms).
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

  # Reconstruct a layer with the given data (and optional remapped mapping).
  # LayerInstance objects cannot be cloned with ggproto(NULL, .), so the layer
  # is rebuilt through layer() with the same geom/stat/mapping/params.
  reconstruct_layer <- function(lyr, data, mapping = NULL) {
    params <- c(lyr$geom_params, lyr$stat_params, lyr$aes_params)
    params <- params[!duplicated(names(params))]
    ggplot2::layer(
      geom = lyr$geom, stat = lyr$stat, data = data,
      mapping = mapping %||% lyr$mapping, position = lyr$position,
      params = params,
      inherit.aes = lyr$inherit.aes,
      show.legend = lyr$show.legend,
      check.aes = FALSE
    )
  }

  new_layers <- plot$layers

  if (length(seq_indices) > 0 && length(layout$seq_arcs) > 0) {
    arc_df <- do.call(rbind, layout$seq_arcs)
    for (idx in seq_indices) new_layers[[idx]] <- reconstruct_layer(plot$layers[[idx]], arc_df)
  }

  if (length(ribbon_indices) > 0 && !is.null(layout$ribbon_polys)) {
    for (idx in ribbon_indices) {
      new_layers[[idx]] <- reconstruct_layer(plot$layers[[idx]], layout$ribbon_polys)
    }
  }

  if (length(gene_poly_indices) > 0 && nrow(layout$gene_polys) > 0) {
    for (idx in gene_poly_indices) new_layers[[idx]] <- reconstruct_layer(plot$layers[[idx]], layout$gene_polys)
  }

  if (length(gene_text_indices) > 0 && nrow(layout$gene_labels) > 0) {
    for (idx in gene_text_indices) new_layers[[idx]] <- reconstruct_layer(plot$layers[[idx]], layout$gene_labels)
  }

  if (length(axis_line_indices) > 0 && nrow(layout$axis_lines) > 0) {
    for (idx in axis_line_indices) new_layers[[idx]] <- reconstruct_layer(plot$layers[[idx]], layout$axis_lines)
  }

  if (length(axis_seg_indices) > 0 && nrow(layout$axis_ticks) > 0) {
    for (idx in axis_seg_indices) new_layers[[idx]] <- reconstruct_layer(plot$layers[[idx]], layout$axis_ticks)
  }

  label_data <- if (nrow(layout$axis_ticks) > 0) {
    subset(layout$axis_ticks, !is.na(label))
  } else {
    layout$axis_ticks
  }
  if (length(axis_text_indices) > 0 && nrow(label_data) > 0) {
    for (idx in axis_text_indices) new_layers[[idx]] <- reconstruct_layer(plot$layers[[idx]], label_data)
  }

  plot$layers <- new_layers

  # ====================================================================
  # Step 4: build and attach scales
  # ====================================================================
  if (length(seq_indices) > 0) {
    plot$scales$add(scale_color_manual(
      name   = "Seq ID",
      values = layout$seq_colors,
      labels = layout$seq_labels,
      breaks = layout$seqs
    ))
  }

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

  # ggplot2 allows only one scale per aesthetic. When both the ribbon and the
  # gene layers are present, the ribbon layers' fill mapping is renamed to the
  # internal aesthetic "fill_ribbon" and the ribbon scale is attached to it;
  # the gene scale keeps the plain "fill" aesthetic, so the two scales do not
  # overwrite each other.
  if (!is.null(ribbon_fill_scale) && !is.null(gene_fill_scale)) {
    ribbon_aes <- "fill_ribbon"
    for (idx in ribbon_indices) {
      lyr <- plot$layers[[idx]]
      mp <- lyr$mapping
      mp_names <- names(mp)
      mp_names[mp_names == "fill"] <- ribbon_aes
      names(mp) <- mp_names
      plot$layers[[idx]] <- reconstruct_layer(lyr, layout$ribbon_polys, mapping = mp)
    }
    s <- ribbon_fill_scale
    s$aesthetics <- ribbon_aes
    # guide_colorbar's available_aes only recognizes fill by default; after the
    # aesthetic rename it must be synchronized, otherwise the colorbar is dropped.
    if (inherits(s$guide, "Guide")) {
      s$guide$available_aes <- gsub("^fill$", ribbon_aes, s$guide$available_aes)
      if (!is.null(s$guide$params$override.aes)) {
        names(s$guide$params$override.aes) <-
          gsub("^fill$", ribbon_aes, names(s$guide$params$override.aes))
      }
    }
    plot$scales$add(s)
    plot$scales$add(gene_fill_scale)
  } else if (!is.null(ribbon_fill_scale)) {
    plot$scales$add(ribbon_fill_scale)
  } else if (!is.null(gene_fill_scale)) {
    plot$scales$add(gene_fill_scale)
  }

  # Ribbon alpha is a preset value; use an identity scale so it renders as specified
  if (!is.null(layout$ribbon_polys)) {
    plot$scales$add(scale_alpha_identity())
  }

  # Axis text size scale
  plot$scales$add(scale_size_identity())

  # ====================================================================
  # Step 5: update the coord range
  # ====================================================================
  ext <- layout$extremes
  pad <- 0.05 * max(ext$x_max - ext$x_min, ext$y_max - ext$y_min, 1)
  plot$coordinates <- coord_fixed(
    ratio = 1,
    xlim  = c(ext$x_min - pad, ext$x_max + pad),
    ylim  = c(ext$y_min - pad, ext$y_max + pad),
    clip  = "off"
  )

  # Run the standard ggplot2 build on the prepared plot.  The ggchord class is
  # removed first so that dispatch proceeds to the base ggplot2 method instead
  # of recursing into this method.
  class(plot) <- setdiff(class(plot), "ggchord")
  ggplot_build(plot)
}
