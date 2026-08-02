# ggchord.R - constructor
# Entry point of the layered ggplot2 API for chord diagrams
# v0.6.0: plot objects are self-contained. Data and parameters are stored on
# the plot object itself, and the layout is computed at build time so that
# print(), ggsave(), ggplot_build() and other ggplot2 workflows all work.

# Global variable declarations (to avoid R CMD check NOTEs)
globalVariables(c(
  "x", "y", "group", "pident", "fill", "strand", "anno", "seq_id",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "label", "label_x", "label_y", "size",
  "fill_col", "alpha"
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
      # Keep the legend key background fixed (white) instead of inheriting the
      # panel background (ggplot2 4.x lets unset legend keys follow
      # panel.background, so a grey/coloured panel would bleed into the legend).
      legend.key        = element_rect(fill = "white", colour = NA),
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
                  show_legend = show_legend, debug = debug),
    # Shared reference environment: layers use it to reach the (latest) plot
    # and to lazily fetch their computed geometry (e.g. for plotly::ggplotly).
    ref    = new.env(),
    layout = NULL
  )
  p$ggchord$ref$plot <- p

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
  # Our geom_* functions return plain lists of layers; flatten them.  Other
  # list-like objects (e.g. themes) must be handled by the standard ggplot2
  # `+` method instead.
  is_layer_list <- is.list(e2) && !inherits(e2, "theme") &&
    length(e2) > 0 && all(vapply(e2, function(el) inherits(el, "LayerInstance"), logical(1)))
  if (is_layer_list) {
    p <- e1
    for (elem in e2) {
      wire_ggchord_layer(elem, e1)
      p <- p + elem
    }
    # Eagerly prepare the plot so that it is self-contained (layout, scales and
    # coordinates attached) even before it is built.  Tools such as
    # plotly::ggplotly() clone the plot before building, so the plot must
    # already carry its scales.
    if (!is.null(p$ggchord$ref)) {
      p <- tryCatch(prepare_ggchord_plot(p), error = function(e) p)
    }
  } else if (inherits(e2, "Scale")) {
    # A user-supplied scale intentionally replaces the ggchord-managed default
    # scale of the same aesthetic; muffle ggplot2's "already present" message.
    p <- withCallingHandlers(
      NextMethod(),
      message = function(m) {
        if (grepl("already present", conditionMessage(m)) &&
            !is.null(findRestart("muffleMessage"))) {
          invokeRestart("muffleMessage")
        }
      }
    )
  } else {
    p <- NextMethod()
  }
  # Keep the shared reference pointing at the latest plot so that lazy layer
  # data (used by plotly::ggplotly and friends) sees the complete plot.
  if (!is.null(e1$ggchord$ref)) e1$ggchord$ref$plot <- p
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

compute_chord_geometry <- function(plot) {
  # Step 1: collect data and parameters from the plot object
  chord <- plot$ggchord
  if (is.null(chord)) {
    stop("Not a valid ggchord object: no data stored on the plot. ",
         "Please build the plot with ggchord().")
  }
  data_list <- chord$data
  global    <- chord$global

  seq_params    <- list()
  ribbon_params <- list()
  gene_params   <- list()
  axis_params   <- list()
  seq_label_params <- list()
  seq_layer_requested <- FALSE

  for (i in seq_along(plot$layers)) {
    pp <- plot$layers[[i]]$ggchord_params
    if (is.null(pp)) next
    switch(pp$type,
      seq       = { seq_params <- pp; seq_layer_requested <- TRUE },
      ribbon    = ribbon_params <- pp,
      gene      = gene_params <- pp,
      axis      = axis_params <- pp,
      seq_label = seq_label_params <- pp
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
  if (!ribbon_color_scheme %in% c("pident", "query", "subject", "single")) {
    stop("ribbon_color_scheme must be 'pident', 'query', 'subject', or 'single'")
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
        subject = {
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
    if (ribbon_color_scheme %in% c("query", "subject")) {
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
  axisGap    <- process_sequence_param(axis_params$axis_gap %||% 0.05,
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

  # --- Process sequence labels ---
  seq_label_text <- NULL
  seq_label_radius <- NULL
  seq_label_rotation <- NULL
  seq_label_size <- NULL
  if (length(seq_label_params) > 0) {
    seq_label_text <- seq_label_params$seq_labels %||% seq_labels
    seq_label_radius <- process_sequence_param(
      seq_label_params$seq_label_radius, seqs, "seq_label_radius", 1.15)
    seq_label_rotation <- process_sequence_param(
      seq_label_params$seq_label_rotation, seqs, "seq_label_rotation", 0)
    seq_label_size <- process_sequence_param(
      seq_label_params$seq_label_size, seqs, "seq_label_size", 3)
  }

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
    seq_label_text = seq_label_text,
    seq_label_radius = seq_label_radius,
    seq_label_rotation = seq_label_rotation,
    seq_label_size = seq_label_size,
    axisGap = axisGap, axisMaj = axisMaj, axisMajLen = axisMajLen,
    axisMin = axisMin, axisMinLen = axisMinLen,
    labelSize = labelSize, labelOffset = labelOffset,
    axisLabelOrientation = axisLabelOri,
    show_axis = show_axis,
    rotation = global$rotation, debug = global$debug
  )

  # Cache the layout so the get_chord_layout() accessor can inspect it and so
  # layers can lazily fetch their geometry (e.g. for plotly::ggplotly).  The
  # shared reference environment is used because it is shared by the plot and
  # all of its layers.
  set_chord_layout(layout)
  if (!is.null(plot$ggchord$ref)) plot$ggchord$ref$layout <- layout
  plot$ggchord$layout <- layout
  layout
}


# ====================================================================
# Shared helpers used by both ggplot_build.ggchord() and the lazy layer
# data path (so that plotly::ggplotly() sees the same scales and geometry).
# ====================================================================

#' Reconstruct a layer with the given data (and optional remapped mapping).
#'
#' LayerInstance objects cannot be cloned with \code{ggproto(NULL, .)}, so the
#' layer is rebuilt through \code{layer()} with the same geom/stat/mapping/params.
#' @keywords internal
reconstruct_layer <- function(lyr, data, mapping = NULL) {
  params <- c(lyr$geom_params, lyr$stat_params, lyr$aes_params)
  params <- params[!duplicated(names(params))]
  new <- ggplot2::layer(
    geom = lyr$geom, stat = lyr$stat, data = data,
    mapping = mapping %||% lyr$mapping, position = lyr$position,
    params = params,
    inherit.aes = lyr$inherit.aes,
    show.legend = lyr$show.legend,
    check.aes = FALSE
  )
  # Preserve the ggchord custom fields on the reconstructed layer
  for (fld in c("ggchord_type", "ggchord_params", "ggchord_placeholder", "ggchord_ref")) {
    if (!is.null(lyr[[fld]])) new[[fld]] <- lyr[[fld]]
  }
  new
}

#' Classify the ggchord layers of a plot by their ggchord_type marker
#' @keywords internal
classify_ggchord_layers <- function(plot) {
  idx <- list(seq = integer(0), ribbon = integer(0), gene_poly = integer(0),
              gene_text = integer(0), axis_line = integer(0),
              axis_seg = integer(0), axis_text = integer(0),
              seq_label = integer(0))
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    type <- lyr$ggchord_type %||% ""
    if (type %in% names(idx)) idx[[type]] <- c(idx[[type]], i)
  }
  idx
}

#' Resolve the per-legend position overrides for a plot
#'
#' Each legend can be moved independently with the `legend_position` argument
#' of `geom_seq()`, `geom_ribbon()` and `geom_gene()`. A NULL entry means that
#' legend follows the theme's `legend.position` together with the others.
#' @keywords internal
ggchord_legend_positions <- function(plot) {
  pos <- list(seq = NULL, ribbon = NULL, gene = NULL)
  for (lyr in plot$layers) {
    pp <- lyr$ggchord_params
    if (is.null(pp)) next
    lp <- pp$legend_position
    if (is.null(lp)) next
    if (!lp %in% c("left", "right", "top", "bottom", "inside")) {
      stop("legend_position must be one of 'left', 'right', 'top', 'bottom' or 'inside'")
    }
    if (!is.null(pp$type) && pp$type %in% names(pos)) pos[[pp$type]] <- lp
  }
  pos
}

#' Build the list of scales for a computed layout
#'
#' @param legend_position The plot theme's `legend.position` (character).
#' @param legend_box The plot theme's `legend.box` setting. When the legend is
#'   at the top/bottom or the legend box is laid out horizontally
#'   (`"horizontal"`), a `unit(1, "null")` colorbar key height collapses to zero
#'   height in ggplot2 (the Identity(%) bar becomes invisible). A fixed size is
#'   used in those cases so the colorbar stays visible; otherwise the colorbar
#'   fills the available height.
#' @param positions Named list with per-legend position overrides
#'   (`seq`, `ribbon`, `gene`), each `NULL` or one of "left", "right", "top",
#'   "bottom", "inside". Overrides make that legend sit in its own legend box at
#'   the given position instead of following the theme's `legend.position`.
#' @keywords internal
make_ggchord_scales <- function(layout, has_seq = FALSE, has_gene = FALSE,
                                legend_position = NULL, legend_box = NULL,
                                positions = list()) {
  scales <- list()

  if (has_seq) {
    scales[[length(scales) + 1]] <- scale_color_manual(
      name   = "Seq ID",
      values = layout$seq_colors,
      labels = layout$seq_labels,
      breaks = layout$seqs,
      guide  = guide_legend(position = positions$seq %||% NULL, order = 1)
    )
  }

  ribbon_fill_scale <- NULL
  if (!is.null(layout$ribbon_polys)) {
    if (layout$ribbon_color_scheme == "pident") {
      # The colorbar follows its effective legend position: vertical and
      # filling the available height at the left/right, horizontal (with a
      # fixed size) at the top/bottom. A "null" key height collapses to zero
      # inside horizontal legend boxes, so a fixed size is used whenever the
      # legend sits at the top/bottom or the box is horizontal.
      ribbon_pos <- positions$ribbon %||% legend_position %||% "right"
      horizontal_legend <- ribbon_pos %in% c("top", "bottom") ||
        identical(legend_box, "horizontal")
      ribbon_fill_scale <- scale_fill_stepsn(
        name    = "Identity(%)",
        colours = layout$ribbon_colors,
        limits  = c(0, 100),
        breaks  = c(0, 50, 80, 90, 95, 100),
        guide   = guide_colorbar(
          position = positions$ribbon %||% NULL,
          theme = theme(
            legend.title.position = "top",
            legend.key.height = if (horizontal_legend) {
              unit(1.5, "cm")
            } else {
              unit(1, "null")
            },
            # A horizontal colorbar needs a longer key; the vertical bar keeps
            # the default key width.
            legend.key.width = if (horizontal_legend) unit(4, "cm") else NULL
          ),
          order = 2
        )
      )
    } else {
      ribbon_fill_scale <- scale_fill_identity()
    }
  }

  gene_fill_scale <- NULL
  if (has_gene) {
    if (layout$gene_color_scheme == "strand") {
      gene_fill_scale <- scale_fill_manual(
        name   = "Strand",
        breaks = c("+", "-"),
        values = layout$gene_pal,
        guide  = guide_legend(position = positions$gene %||% NULL, order = 3)
      )
    } else {
      gene_fill_scale <- scale_fill_manual(
        name   = "Gene Annotation",
        breaks = layout$final_gene_order,
        values = layout$gene_pal,
        guide  = guide_legend(position = positions$gene %||% NULL, order = 3)
      )
    }
  }

  # ggplot2 allows only one scale per aesthetic. When both the ribbon and the
  # gene layers are present, the ribbon scale is attached to the internal
  # aesthetic "fill_ribbon" and the gene scale keeps "fill", so they do not
  # overwrite each other.
  ribbon_aes <- "fill"
  if (!is.null(ribbon_fill_scale) && !is.null(gene_fill_scale)) {
    ribbon_aes <- "zfill"
    s <- ribbon_fill_scale
    s$aesthetics <- ribbon_aes
    if (inherits(s$guide, "Guide")) {
      s$guide$available_aes <- gsub("^fill$", ribbon_aes, s$guide$available_aes)
      if (!is.null(s$guide$params$override.aes)) {
        names(s$guide$params$override.aes) <-
          gsub("^fill$", ribbon_aes, names(s$guide$params$override.aes))
      }
    }
    scales[[length(scales) + 1]] <- s
    scales[[length(scales) + 1]] <- gene_fill_scale
  } else if (!is.null(ribbon_fill_scale)) {
    scales[[length(scales) + 1]] <- ribbon_fill_scale
  } else if (!is.null(gene_fill_scale)) {
    scales[[length(scales) + 1]] <- gene_fill_scale
  }

  # Ribbon alpha is a preset value; use an identity scale so it renders as specified
  if (!is.null(layout$ribbon_polys)) {
    scales[[length(scales) + 1]] <- scale_alpha_identity()
  }
  # Axis text size scale
  scales[[length(scales) + 1]] <- scale_size_identity()

  list(scales = scales, ribbon_aes = ribbon_aes)
}

#' Add scales to a plot, respecting user-supplied scales
#' @keywords internal
attach_ggchord_scales <- function(plot, scales) {
  for (s in scales) {
    aes <- s$aesthetics[1]
    if (!is.null(aes) && !plot$scales$has_scale(aes)) {
      s$ggchord_managed <- TRUE
      plot$scales$add(s)
    }
  }
  plot
}

#' Rename the ribbon layers' fill mapping to the internal ribbon aesthetic
#' @keywords internal
rename_ribbon_layers <- function(plot, ribbon_indices, ribbon_aes, layout) {
  if (ribbon_aes != "fill" && length(ribbon_indices) > 0 &&
      !is.null(layout$ribbon_polys)) {
    for (idx in ribbon_indices) {
      lyr <- plot$layers[[idx]]
      mp <- lyr$mapping
      mp_names <- names(mp)
      mp_names[mp_names == "fill"] <- ribbon_aes
      names(mp) <- mp_names
      plot$layers[[idx]] <- reconstruct_layer(lyr, layout$ribbon_polys, mapping = mp)
    }
  }
  plot
}

#' Set the fixed coordinate system from the layout extremes
#' @keywords internal
set_ggchord_coord <- function(plot, layout) {
  ext <- layout$extremes
  pad <- 0.05 * max(ext$x_max - ext$x_min, ext$y_max - ext$y_min, 1)
  plot$coordinates <- coord_fixed(
    ratio = 1,
    xlim  = c(ext$x_min - pad, ext$x_max + pad),
    ylim  = c(ext$y_min - pad, ext$y_max + pad),
    clip  = "off"
  )
  plot
}

#' Fully prepare a ggchord plot and return it (compute layout, rename ribbon
#' mappings, attach scales, set coordinates). The layout is cached on the plot
#' (and on the shared reference environment) during preparation. Used by the
#' lazy layer data path so that plotly::ggplotly() sees the same state as a
#' normal build.
#' @keywords internal
prepare_ggchord_plot <- function(plot) {
  plot$scales$scales <- Filter(function(s) is.null(s$ggchord_managed),
                               plot$scales$scales)
  layout <- compute_chord_geometry(plot)
  cls <- classify_ggchord_layers(plot)
  sc <- make_ggchord_scales(layout,
                            has_seq = length(cls$seq) > 0,
                            has_gene = nrow(layout$gene_polys) > 0,
                            legend_position = plot$theme$legend.position,
                            legend_box = plot$theme$legend.box,
                            positions = ggchord_legend_positions(plot))
  plot <- rename_ribbon_layers(plot, cls$ribbon, sc$ribbon_aes, layout)
  plot <- attach_ggchord_scales(plot, sc$scales)
  plot <- set_ggchord_coord(plot, layout)
  plot
}

#' @export
ggplot_build.ggchord <- function(plot, ...) {
  chord <- plot$ggchord
  if (is.null(chord)) {
    stop("Not a valid ggchord object: no data stored on the plot. ",
         "Please build the plot with ggchord().")
  }
  # The plot object is self-contained: it carries its own scales (tagged with
  # ggchord_managed) so that tools such as plotly::ggplotly() that clone the
  # plot before building see the correct scales.  ggchord-managed scales are
  # refreshed on every build (user-supplied scales are kept).
  plot$scales$scales <- Filter(function(s) is.null(s$ggchord_managed),
                               plot$scales$scales)

  layout <- compute_chord_geometry(plot)

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
  seq_label_indices <- integer(0)

  # Layers are tagged with a ggchord_type marker at creation, so they can be
  # classified even before their (lazily computed) data exists.
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    switch(lyr$ggchord_type %||% "",
      seq       = seq_indices <- c(seq_indices, i),
      ribbon    = ribbon_indices <- c(ribbon_indices, i),
      gene_poly = gene_poly_indices <- c(gene_poly_indices, i),
      gene_text = gene_text_indices <- c(gene_text_indices, i),
      axis_line = axis_line_indices <- c(axis_line_indices, i),
      axis_seg  = axis_seg_indices <- c(axis_seg_indices, i),
      axis_text = axis_text_indices <- c(axis_text_indices, i),
      seq_label = seq_label_indices <- c(seq_label_indices, i)
    )
  }

  # Reconstruct every ggchord layer with its computed geometry (or its
  # placeholder data when the geometry is empty).  This replaces the lazy data
  # functions so that a normal build leaves the plot fully concrete.
  new_layers <- plot$layers
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (is.null(lyr$ggchord_type)) next
    new_layers[[i]] <- reconstruct_layer(lyr, extract_ggchord_layer_data(lyr, layout))
  }
  plot$layers <- new_layers


  # ====================================================================
  # Step 4: build and attach scales
  # ====================================================================
  sc <- make_ggchord_scales(layout,
                            has_seq = length(seq_indices) > 0,
                            has_gene = nrow(layout$gene_polys) > 0,
                            legend_position = plot$theme$legend.position,
                            legend_box = plot$theme$legend.box,
                            positions = ggchord_legend_positions(plot))
  plot <- rename_ribbon_layers(plot, ribbon_indices, sc$ribbon_aes, layout)
  plot <- attach_ggchord_scales(plot, sc$scales)

  # ====================================================================
  # Step 5: update the coord range
  # ====================================================================
  plot <- set_ggchord_coord(plot, layout)

  # Run the standard ggplot2 build on the prepared plot.  The ggchord class is
  # removed first so that dispatch proceeds to the base ggplot2 method instead
  # of recursing into this method.
  class(plot) <- setdiff(class(plot), "ggchord")
  ggplot_build(plot)
}
