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
  "fill_col", "alpha", "label_hjust", "label_vjust", "label_angle",
  "linetype", "zcolour", "zregionfill", "zoutline", "zlinetype",
  "outline_col", "linetype_val", "value", "source_row", "direction"
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
#' @param validate Character, default \code{"warn"}. How to run the structured
#'   input-data validation (see \code{\link{validate_ggchord_data}}):
#'   \code{"warn"} emits a single summary warning when the data has problems
#'   and caches the full report on the plot object
#'   (\code{p$ggchord$validation}); \code{"error"} stops on severe problems;
#'   \code{"none"} skips the diagnostic validation (the cheap structural
#'   checks that prevent crashes are still performed).
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
    debug = FALSE,
    validate = c("warn", "error", "none")
) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  validate <- match.arg(validate)
  # ====================================================================
  # 1. Validate data
  # ====================================================================
  if (!is.numeric(rotation) || length(rotation) != 1 || !is.finite(rotation)) {
    ggchord_stop("rotation must be a finite numeric value")
  }
  if (!is.logical(show_legend) || length(show_legend) != 1 || is.na(show_legend)) {
    ggchord_stop("show_legend must be TRUE or FALSE")
  }
  if (!is.logical(debug) || length(debug) != 1 || is.na(debug)) {
    ggchord_stop("debug must be TRUE or FALSE")
  }
  if (!is.data.frame(seq_data) || nrow(seq_data) == 0) {
    ggchord_stop("seq_data must be a non-empty data.frame")
  }
  required_seq_cols <- c("seq_id", "length")
  if (!all(required_seq_cols %in% colnames(seq_data))) {
    ggchord_stop("seq_data must contain the following columns: ", paste(required_seq_cols, collapse = ", "))
  }
  if (anyNA(seq_data$seq_id) || any(!nzchar(as.character(seq_data$seq_id)))) {
    ggchord_stop("The 'seq_id' values in seq_data must be non-missing and non-empty")
  }
  if (!is.numeric(seq_data$length) || any(!is.finite(seq_data$length)) ||
      any(seq_data$length <= 0)) {
    ggchord_stop("The 'length' values in seq_data must be finite positive numbers")
  }
  if (anyDuplicated(seq_data$seq_id)) {
    ggchord_stop("The 'seq_id' values in seq_data must be unique")
  }

  seq_lens <- setNames(seq_data$length, seq_data$seq_id)

  if (!is.null(ribbon_data)) {
    if (!is.data.frame(ribbon_data)) {
      ggchord_stop("ribbon_data must be a data.frame")
    }
    required_ribbon_cols <- c("qaccver", "saccver", "length", "pident",
                              "qstart", "qend", "sstart", "send")
    if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
      ggchord_stop("ribbon_data must contain the following columns: ",
           paste(required_ribbon_cols, collapse = ", "))
    }
    numeric_ribbon_cols <- c("length", "pident", "qstart", "qend", "sstart", "send")
    if (any(vapply(ribbon_data[numeric_ribbon_cols],
                   function(x) !is.numeric(x) || any(!is.finite(x)), logical(1)))) {
      ggchord_stop("ribbon_data numeric columns (length, pident, qstart, qend, sstart, send) must contain finite numbers")
    }
    if (anyNA(ribbon_data$qaccver) || anyNA(ribbon_data$saccver) ||
        any(!nzchar(as.character(ribbon_data$qaccver))) ||
        any(!nzchar(as.character(ribbon_data$saccver)))) {
      ggchord_stop("ribbon_data sequence IDs must be non-missing and non-empty")
    }
    if (any(ribbon_data$length <= 0)) {
      ggchord_stop("The 'length' values in ribbon_data must be positive")
    }
    if (nrow(ribbon_data) == 0) warning("No valid alignment data in ribbon_data")
    if (debug) cat("Number of alignment data rows: ", nrow(ribbon_data), "\n")
  }

  if (!is.null(gene_data)) {
    if (!is.data.frame(gene_data)) {
      ggchord_stop("gene_data must be a data.frame")
    }
    required_gene_cols <- c("seq_id", "start", "end", "strand", "anno")
    if (!all(required_gene_cols %in% colnames(gene_data))) {
      ggchord_stop("gene_data must contain the following columns: ",
           paste(required_gene_cols, collapse = ", "))
    }
    if (!is.numeric(gene_data$start) || !is.numeric(gene_data$end) ||
        any(!is.finite(gene_data$start)) || any(!is.finite(gene_data$end))) {
      ggchord_stop("The 'start' and 'end' values in gene_data must be finite numbers")
    }
    if (anyNA(gene_data$seq_id) || any(!nzchar(as.character(gene_data$seq_id)))) {
      ggchord_stop("gene_data sequence IDs must be non-missing and non-empty")
    }
    if (anyNA(gene_data$strand) || any(!gene_data$strand %in% c("+", "-"))) {
      ggchord_stop("The 'strand' values in gene_data can only be '+' or '-'")
    }
    if (nrow(gene_data) == 0) warning("No valid gene annotation data in gene_data")
    if (debug) cat("Number of gene annotation rows: ", nrow(gene_data), "\n")
  }

  # ====================================================================
  # 1b. Structured validation (v0.7.0).  A single summary warning is emitted
  #     for "warn" (never one warning per row); "error" stops on severe
  #     problems; "none" keeps the fast path (structural checks above still
  #     prevent internal crashes).  The full report is cached on the plot.
  # ====================================================================
  validation <- NULL
  if (validate != "none") {
    validation <- validate_ggchord_data(seq_data, ribbon_data, gene_data,
                                        strict = FALSE)
    if (validate == "error" && !validation$valid) {
      ggchord_stop(sprintf(
        "ggchord(): input data failed validation (%d severe error(s); first: %s). Run validate_ggchord_data(..., strict = FALSE) for the full report.",
        nrow(validation$errors), validation$errors$message[1]),
        call. = FALSE)
    }
    if (validate == "warn") {
      n_err <- nrow(validation$errors)
      n_warn <- nrow(validation$warnings)
      if (n_err > 0) {
        warning(sprintf(
          "ggchord(): input data has %d severe validation error(s) (e.g. \"%s\"). The plot may be misleading; run validate_ggchord_data(..., strict = FALSE) for details.",
          n_err, validation$errors$message[1]), call. = FALSE)
      } else if (n_warn > 0) {
        warning(sprintf(
          "ggchord(): input data has %d validation warning(s) (e.g. \"%s\"). Run validate_ggchord_data(...) for details.",
          n_warn, validation$warnings$message[1]), call. = FALSE)
      }
    }
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
      # Transparent legend keys: the background (if any) is the plot/page
      # background, not the panel, so legends blend into any theme.
      legend.key        = element_rect(fill = NA, colour = NA),
      legend.box.spacing  = unit(10, "mm"),
      legend.spacing      = unit(5, "mm"),
      legend.text         = element_text(size = 8),
      legend.title        = element_text(size = 10, face = "bold"),
      axis.title          = element_blank(),
      axis.line           = element_blank(),
      axis.ticks          = element_blank(),
      axis.text           = element_blank(),
      panel.background    = element_blank(),
      panel.grid          = element_blank(),
      legend.position     = if (isTRUE(show_legend)) "right" else "none"
    )

  p$ggchord <- list(
    data   = list(seq_data = seq_data, ribbon_data = ribbon_data,
                  gene_data = gene_data),
    global = list(rotation = rotation, panel_margin = panel_margin,
                  show_legend = show_legend, debug = debug,
                  validate = validate),
    validation = validation,
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
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

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
    # Layout/scales/coordinates are computed lazily on the first build or
    # plotly conversion (ggchord_layer_data() and ggchord_plotly_ggplot() both
    # call prepare_ggchord_plot() on demand). Eagerly preparing here would
    # recompute the full chord geometry after every added layer, which is very
    # expensive for large alignment tables.
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
    ggchord_stop("Not a valid ggchord object: no data stored on the plot. ",
         "Please build the plot with ggchord().")
  }
  data_list <- chord$data
  global    <- chord$global

  seq_params    <- list()
  ribbon_params <- list()
  gene_params   <- list()
  gene_data_override <- NULL
  gene_label_params <- list()
  gene_repel_params <- list()
  axis_params   <- list()
  seq_label_params <- list()
  seq_region_params <- list()
  ribbon_highlight_params <- list()
  seq_layer_requested <- FALSE
  gene_label_layer <- FALSE
  gene_repel_layer <- FALSE

  for (i in seq_along(plot$layers)) {
    pp <- plot$layers[[i]]$ggchord_params
    if (is.null(pp)) next
    switch(pp$type,
      seq               = { seq_params <- pp; seq_layer_requested <- TRUE },
      ribbon            = ribbon_params <- pp,
      gene              = {
        gene_params <- pp
        if (!is.null(pp$gene_data_override)) gene_data_override <- pp$gene_data_override
      },
      gene_label        = { gene_label_params <- pp; gene_label_layer <- TRUE },
      gene_label_repel  = { gene_repel_params <- pp; gene_repel_layer <- TRUE },
      axis              = axis_params <- pp,
      seq_label         = seq_label_params <- pp,
      seq_region        = seq_region_params <- pp,
      ribbon_highlight  = ribbon_highlight_params <- pp
    )
  }

  # --- Process sequences ---
  seq_data <- data_list$seq_data
  seqs     <- seq_data$seq_id
  lens     <- setNames(seq_data$length, seqs)

  if (!is.null(seq_params$seq_order)) {
    if (!all(seq_params$seq_order %in% seqs)) {
      ggchord_stop("seq_order contains unknown sequence IDs")
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

  if (!is.numeric(seqRadius) || any(!is.finite(seqRadius)) || any(seqRadius <= 0)) {
    ggchord_stop("seq_radius must contain finite positive numbers")
  }
  if (!is.numeric(orientation) || any(!is.finite(orientation)) ||
      any(!orientation %in% c(-1, 1))) {
    ggchord_stop("seq_orientation can only be 1 or -1")
  }
  if (!is.numeric(seq_gap) || any(!is.finite(seq_gap)) ||
      any(seq_gap < 0 | seq_gap >= 0.5)) {
    ggchord_stop("seq_gap must be in the [0, 0.5) range")
  }
  if (!is.numeric(seq_curvature) || any(!is.finite(seq_curvature))) {
    ggchord_stop("seq_curvature must contain finite numbers")
  }

  if (!is.null(seq_params$seq_colors)) {
    seq_colors <- process_sequence_param(seq_params$seq_colors, seqs, "seq_colors")
  } else {
    pal <- chord_default_palette(n)
    seq_colors <- setNames(pal, seqs)
  }

  # --- Process sequence grouping (v0.8.0) ---
  seq_group <- resolve_ggchord_seq_group(seq_data, seqs, seq_params$seq_group)
  seq_group_gap <- seq_params$seq_group_gap %||% 0.08
  seq_group_labels <- seq_params$seq_group_labels %||% TRUE
  seq_group_label_radius <- seq_params$seq_group_label_radius %||% 1.35
  seq_group_colors <- seq_params$seq_group_colors

  if (!is.numeric(seq_group_gap) || length(seq_group_gap) != 1 ||
      !is.finite(seq_group_gap) || seq_group_gap < 0) {
    ggchord_stop("seq_group_gap must be a finite non-negative number")
  }
  if (!is.numeric(seq_group_label_radius) || length(seq_group_label_radius) != 1 ||
      !is.finite(seq_group_label_radius)) {
    ggchord_stop("seq_group_label_radius must be a finite number")
  }
  if (!is.null(seq_group)) {
    if (length(seq_group) != length(seqs) ||
        anyNA(seq_group) || any(!nzchar(as.character(seq_group)))) {
      ggchord_stop("seq_group must provide a non-missing, non-empty group for every sequence")
    }
    seq_group <- as.character(seq_group)
  } else if (!is.null(seq_group_colors)) {
    ggchord_stop("seq_group_colors requires seq_group (add a seq_group column or pass geom_seq(seq_group = ...))")
  }

  # --- Process ribbons ---
  ribbonGap  <- process_sequence_param(ribbon_params$ribbon_gap %||% 0.15,
                                       seqs, "ribbon_gap", 0.15)
  ribbon_color_scheme <- ribbon_params$ribbon_color_scheme %||% "pident"
  ribbon_color_by     <- ribbon_params$ribbon_color_by
  ribbon_color_limits <- ribbon_params$ribbon_color_limits
  ribbon_color_breaks <- ribbon_params$ribbon_color_breaks
  ribbon_color_name   <- ribbon_params$ribbon_color_name
  ribbon_alpha    <- ribbon_params$ribbon_alpha %||% 0.35
  ribbon_alpha_by <- ribbon_params$ribbon_alpha_by
  ribbon_alpha_range <- ribbon_params$ribbon_alpha_range %||% c(0.15, 0.9)
  ribbon_ctrl_pt  <- ribbon_params$ribbon_ctrl_point %||% c(0, 0)
  ribbon_outline_by   <- ribbon_params$ribbon_outline_by
  ribbon_outline_colors <- ribbon_params$ribbon_outline_colors
  ribbon_linetype_by  <- ribbon_params$ribbon_linetype_by
  ribbon_linetypes    <- ribbon_params$ribbon_linetypes
  ribbon_direction    <- ribbon_params$ribbon_direction %||% "none"
  ribbon_direction_colors <- ribbon_params$ribbon_direction_colors %||% c(same = "black", reverse = "grey50")
  ribbon_direction_linetypes <- ribbon_params$ribbon_direction_linetypes %||% c(same = "solid", reverse = "dashed")
  ribbon_direction_alpha <- ribbon_params$ribbon_direction_alpha %||% c(same = 1, reverse = 0.45)

  ribbon_colors <- ribbon_params$ribbon_colors
  if (!ribbon_color_scheme %in% c("pident", "query", "subject", "single")) {
    ggchord_stop("ribbon_color_scheme must be 'pident', 'query', 'subject', or 'single'")
  }
  if (!is.null(ribbon_color_by)) {
    ribbon_color_scheme <- "value"
    ribbon_color_name <- ribbon_color_name %||% ribbon_color_by
  }
  if (!is.numeric(ribbon_alpha) || length(ribbon_alpha) != 1 ||
      !is.finite(ribbon_alpha) || ribbon_alpha < 0 || ribbon_alpha > 1) {
    ggchord_stop("ribbon_alpha must be in the [0, 1] range")
  }
  if (!is.numeric(ribbonGap) || any(!is.finite(ribbonGap))) {
    ggchord_stop("ribbon_gap must contain finite numbers")
  }
  if (!is.null(ribbon_color_limits) &&
      (!is.numeric(ribbon_color_limits) || length(ribbon_color_limits) != 2 ||
       !is.finite(ribbon_color_limits[1]) || !is.finite(ribbon_color_limits[2]) ||
       ribbon_color_limits[1] >= ribbon_color_limits[2])) {
    ggchord_stop("ribbon_color_limits must be a length-2 increasing numeric vector")
  }
  if (!is.null(ribbon_color_breaks) &&
      (!is.numeric(ribbon_color_breaks) || any(!is.finite(ribbon_color_breaks)))) {
    ggchord_stop("ribbon_color_breaks must be a finite numeric vector")
  }
  if (!is.numeric(ribbon_alpha_range) || length(ribbon_alpha_range) != 2 ||
      any(!is.finite(ribbon_alpha_range)) || ribbon_alpha_range[1] < 0 ||
      ribbon_alpha_range[2] > 1 || ribbon_alpha_range[1] > ribbon_alpha_range[2]) {
    ggchord_stop("ribbon_alpha_range must be two increasing values within [0, 1]")
  }
  if (!ribbon_direction %in% c("none", "alpha", "outline", "linetype")) {
    ggchord_stop("ribbon_direction must be 'none', 'alpha', 'outline', or 'linetype'")
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
                    "#B4DE2CFF","#FDE725FF"),
        value = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF",
                    "#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF",
                    "#B4DE2CFF","#FDE725FF"))
    }
    if (ribbon_color_scheme %in% c("query", "subject")) {
      ribbon_colors <- process_sequence_param(ribbon_colors, seqs,
                                              "ribbon_colors")
    } else if (ribbon_color_scheme %in% c("pident", "value") && length(ribbon_colors) < 2) {
      ggchord_stop("The 'pident'/'value' scheme requires at least two ribbon_colors")
    } else if (ribbon_color_scheme == "single" && length(ribbon_colors) < 1) {
      ggchord_stop("The 'single' scheme requires at least one ribbon_colors")
    }

    if (!is.null(ribbon_color_by)) {
      if (!ribbon_color_by %in% colnames(data_list$ribbon_data)) {
        ggchord_stop("ribbon_color_by column '", ribbon_color_by, "' not found in ribbon_data")
      }
      if (!is.numeric(data_list$ribbon_data[[ribbon_color_by]]) ||
          any(!is.finite(data_list$ribbon_data[[ribbon_color_by]]))) {
        ggchord_stop("ribbon_color_by column '", ribbon_color_by, "' must be numeric and finite")
      }
    }
    if (!is.null(ribbon_alpha_by)) {
      if (!ribbon_alpha_by %in% colnames(data_list$ribbon_data)) {
        ggchord_stop("ribbon_alpha_by column '", ribbon_alpha_by, "' not found in ribbon_data")
      }
      if (!is.numeric(data_list$ribbon_data[[ribbon_alpha_by]]) ||
          any(!is.finite(data_list$ribbon_data[[ribbon_alpha_by]]))) {
        ggchord_stop("ribbon_alpha_by column '", ribbon_alpha_by, "' must be numeric and finite")
      }
    }
    if (!is.null(ribbon_outline_by) &&
        !ribbon_outline_by %in% colnames(data_list$ribbon_data)) {
      ggchord_stop("ribbon_outline_by column '", ribbon_outline_by, "' not found in ribbon_data")
    }
    if (!is.null(ribbon_linetype_by) &&
        !ribbon_linetype_by %in% colnames(data_list$ribbon_data)) {
      ggchord_stop("ribbon_linetype_by column '", ribbon_linetype_by, "' not found in ribbon_data")
    }
  }

  # --- Process genes ---
  gene_off  <- gene_params$gene_offset %||% 0.1
  gene_w    <- gene_params$gene_width %||% 0.05
  gene_cs   <- gene_params$gene_color_scheme %||% "strand"
  gene_cols <- gene_params$gene_colors
  gene_ord  <- gene_params$gene_order
  # Gene label settings come from the dedicated geom_gene_label() layer, with
  # the legacy geom_gene() arguments as fallback.
  # The repel layer takes priority over the fixed label layer
  lbl <- if (gene_repel_layer) gene_repel_params else gene_label_params
  gene_ls   <- gene_label_layer || gene_repel_layer ||
    isTRUE(gene_params$show_label_override) ||
    isTRUE(gene_params$gene_label_show)
  gene_lsz  <- lbl$gene_label_size %||%
    gene_params$label_size_override %||%
    gene_params$gene_label_size %||% 2.5
  gene_lr   <- lbl$gene_label_rotation %||%
    gene_params$gene_label_rotation %||% 0
  gene_lro  <- lbl$gene_label_radial_offset %||%
    gene_params$gene_label_radial_offset %||% 0
  gene_lco  <- lbl$gene_label_circum_offset %||%
    gene_params$gene_label_circum_offset %||% 0
  gene_lcl  <- if (is.null(lbl$gene_label_circum_limit)) {
    if (is.null(gene_params$gene_label_circum_limit)) TRUE
    else gene_params$gene_label_circum_limit
  } else lbl$gene_label_circum_limit
  gene_lwrap  <- lbl$gene_label_wrap %||% gene_params$gene_label_wrap
  gene_lrepel_layer <- gene_repel_layer
  gene_lrepel_maxov <- gene_repel_params$max_overlaps %||% Inf
  gene_lrepel_box   <- gene_repel_params$box_padding %||% 0.25
  gene_lrepel_pt    <- gene_repel_params$point_padding %||% 0.1
  gene_lrepel_minseg <- gene_repel_params$min_segment_length %||% 0.5
  gene_lrepel_force <- gene_repel_params$force %||% 1
  gene_lrepel_seed  <- gene_repel_params$seed %||% 123
  gene_lrepel_orient <- gene_repel_params$gene_label_orientation %||% "horizontal"
  gene_lrepel_seg    <- gene_repel_params$gene_label_segment %||% "elbow"
  gene_lrepel_side   <- gene_repel_params$gene_label_side %||% "outside"
  gene_lrepel_ltype  <- gene_repel_params$gene_label_segment_linetype %||% "auto"

  if (!gene_cs %in% c("strand", "manual")) {
    ggchord_stop("gene_color_scheme must be 'strand' or 'manual'")
  }
  if (!is.numeric(gene_lsz) || length(gene_lsz) != 1 || !is.finite(gene_lsz) ||
      gene_lsz <= 0) {
    ggchord_stop("gene_label_size must be a finite positive number")
  }
  if (!is.null(gene_lwrap) && (!is.numeric(gene_lwrap) ||
      length(gene_lwrap) != 1 || !is.finite(gene_lwrap) || gene_lwrap < 0)) {
    ggchord_stop("gene_label_wrap must be NULL or a finite non-negative number")
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
  axisMaj    <- process_sequence_param(axis_params$axis_tick_major_number %||% 3,
                                       seqs, "axis_tick_major_number", 3)
  axisMajLen <- process_sequence_param(axis_params$axis_tick_major_length %||% 0.02,
                                       seqs, "axis_tick_major_length", 0.02)
  axisMin    <- process_sequence_param(axis_params$axis_tick_minor_number %||% 4,
                                       seqs, "axis_tick_minor_number", 4)
  axisMinLen <- process_sequence_param(axis_params$axis_tick_minor_length %||% 0.01,
                                       seqs, "axis_tick_minor_length", 0.01)
  labelSize  <- process_sequence_param(axis_params$axis_label_size %||% 3,
                                       seqs, "axis_label_size", 3)
  labelOffset <- process_sequence_param(axis_params$axis_label_offset %||% 2,
                                        seqs, "axis_label_offset", 2)
  axisLabelHide <- isTRUE(axis_params$axis_label_hide_overlaps)
  axisLabelOri <- process_axis_orientation(
    axis_params$axis_label_orientation %||% "parallel", seqs
  )
  if (!is.logical(show_axis) || length(show_axis) != 1 || is.na(show_axis)) {
    ggchord_stop("show_axis must be TRUE or FALSE")
  }
  axis_numeric <- list(
    axis_gap = axisGap,
    axis_tick_major_length = axisMajLen,
    axis_tick_minor_length = axisMinLen,
    axis_label_size = labelSize,
    axis_label_offset = labelOffset
  )
  for (nm in names(axis_numeric)) {
    value <- axis_numeric[[nm]]
    if (!is.numeric(value) || any(!is.finite(value))) {
      ggchord_stop(nm, " must contain finite numbers")
    }
  }
  if (!is.numeric(axisMaj) || any(!is.finite(axisMaj)) ||
      any(axisMaj < 1 | axisMaj != as.integer(axisMaj))) {
    ggchord_stop("axis_tick_major_number must contain positive integers")
  }
  if (!is.numeric(axisMin) || any(!is.finite(axisMin)) ||
      any(axisMin < 0 | axisMin != as.integer(axisMin))) {
    ggchord_stop("axis_tick_minor_number must contain non-negative integers")
  }

  # --- Process sequence labels ---
  seq_label_text <- NULL
  seq_label_radius <- NULL
  seq_label_rotation <- NULL
  seq_label_size <- NULL
  seq_label_orientation <- "arc"
  seq_label_hjust <- NULL
  seq_label_vjust <- NULL
  if (length(seq_label_params) > 0) {
    seq_label_text <- if (is.null(seq_label_params$seq_labels)) {
      seq_labels
    } else {
      # Process through the standard parameter helper so that unnamed vectors
      # are matched positionally to the sequences (named by seq_id).
      process_sequence_param(seq_label_params$seq_labels, seqs, "seq_labels",
                             default_value = seqs)
    }
    seq_label_radius <- process_sequence_param(
      seq_label_params$seq_label_radius, seqs, "seq_label_radius", 1)
    seq_label_rotation <- process_sequence_param(
      seq_label_params$seq_label_rotation, seqs, "seq_label_rotation", 0)
    seq_label_size <- process_sequence_param(
      seq_label_params$seq_label_size, seqs, "seq_label_size", 3)
    seq_label_orientation <- seq_label_params$seq_label_orientation %||% "arc"
    seq_label_hjust <- if (is.null(seq_label_params$seq_label_hjust)) {
      if (identical(seq_label_orientation, "arc")) {
        process_sequence_param(-0.2, seqs, "seq_label_hjust", -0.2)
      } else {
        NULL
      }
    } else {
      process_sequence_param(seq_label_params$seq_label_hjust, seqs,
                             "seq_label_hjust", 0.5)
    }
    seq_label_vjust <- if (is.null(seq_label_params$seq_label_vjust)) {
      NULL
    } else {
      process_sequence_param(seq_label_params$seq_label_vjust, seqs,
                             "seq_label_vjust", 0.5)
    }
    seq_label_numeric <- list(
      seq_label_radius = seq_label_radius,
      seq_label_rotation = seq_label_rotation,
      seq_label_size = seq_label_size,
      seq_label_hjust = seq_label_hjust,
      seq_label_vjust = seq_label_vjust
    )
    for (nm in names(seq_label_numeric)) {
      value <- seq_label_numeric[[nm]]
      if (!is.null(value) && (!is.numeric(value) || any(!is.finite(value)))) {
        ggchord_stop(nm, " must contain finite numbers")
      }
    }
    if (any(seq_label_size <= 0)) {
      ggchord_stop("seq_label_size must contain positive numbers")
    }
  }

  # --- Process sequence-region highlight data ---
  region_data <- seq_region_params$regions
  region_fill   <- seq_region_params$region_fill %||% "#F59E0B"
  region_color  <- seq_region_params$region_color %||% "#B45309"
  region_alpha  <- seq_region_params$region_alpha %||% 0.25
  region_width  <- seq_region_params$region_width %||% 0.08
  region_offset <- seq_region_params$region_offset %||% 0
  region_side   <- seq_region_params$region_side %||% "inside"
  if (!is.null(region_data) && !is.data.frame(region_data)) {
    ggchord_stop("geom_seq_region(): regions must be a data.frame")
  }

  # --- Compute ribbon highlight selection (safe, no string evaluation) ---
  ribbon_highlight_rows <- integer(0)
  if (length(ribbon_highlight_params) > 0 &&
      !is.null(data_list$ribbon_data) && nrow(data_list$ribbon_data) > 0) {
    rd <- data_list$ribbon_data
    keep <- rep(TRUE, nrow(rd))

    if (!is.null(ribbon_highlight_params$ribbon_ids)) {
      ids <- ribbon_highlight_params$ribbon_ids
      if (!is.numeric(ids)) {
        ggchord_stop("geom_ribbon_highlight(): ribbon_ids must be numeric row numbers")
      }
      keep <- keep & seq_len(nrow(rd)) %in% as.integer(ids)
    }
    if (!is.null(ribbon_highlight_params$qaccver)) {
      keep <- keep & rd$qaccver %in% ribbon_highlight_params$qaccver
    }
    if (!is.null(ribbon_highlight_params$saccver)) {
      keep <- keep & rd$saccver %in% ribbon_highlight_params$saccver
    }
    if (!is.null(ribbon_highlight_params$min_pident)) {
      keep <- keep & rd$pident >= ribbon_highlight_params$min_pident
    }
    if (!is.null(ribbon_highlight_params$max_pident)) {
      keep <- keep & rd$pident <= ribbon_highlight_params$max_pident
    }
    if (!is.null(ribbon_highlight_params$min_length)) {
      keep <- keep & rd$length >= ribbon_highlight_params$min_length
    }
    if (!is.null(ribbon_highlight_params$max_length)) {
      keep <- keep & rd$length <= ribbon_highlight_params$max_length
    }
    pred <- ribbon_highlight_params$predicate
    if (!is.null(pred)) {
      pred_res <- tryCatch(as.logical(pred(rd)),
                           error = function(e) {
                             ggchord_stop("geom_ribbon_highlight(): predicate failed: ",
                                          conditionMessage(e))
                           })
      if (length(pred_res) != nrow(rd)) {
        ggchord_stop("geom_ribbon_highlight(): predicate must return one logical value per ribbon row")
      }
      keep <- keep & pred_res
    }
    ribbon_highlight_rows <- which(keep)
  }

  # ====================================================================
  # Step 2: compute the layout
  # ====================================================================
  layout <- compute_chord_layout(
    seqs = seqs, lens = lens, seq_labels = seq_labels,
    seq_colors = seq_colors, seqRadius = seqRadius,
    seq_curvature = seq_curvature, orientation = orientation,
    seq_gap = seq_gap,
    seq_group = seq_group,
    seq_group_gap = seq_group_gap,
    seq_group_labels = seq_group_labels,
    seq_group_label_radius = seq_group_label_radius,
    seq_group_colors = seq_group_colors,
    ribbon_data = data_list$ribbon_data, ribbonGap = ribbonGap,
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors = ribbon_colors, ribbon_alpha = ribbon_alpha,
    ribbon_color_by = ribbon_color_by,
    ribbon_color_limits = ribbon_color_limits,
    ribbon_color_breaks = ribbon_color_breaks,
    ribbon_color_name = ribbon_color_name,
    ribbon_alpha_by = ribbon_alpha_by,
    ribbon_alpha_range = ribbon_alpha_range,
    ribbon_outline_by = ribbon_outline_by,
    ribbon_outline_colors = ribbon_outline_colors,
    ribbon_linetype_by = ribbon_linetype_by,
    ribbon_linetypes = ribbon_linetypes,
    ribbon_direction = ribbon_direction,
    ribbon_direction_colors = ribbon_direction_colors,
    ribbon_direction_linetypes = ribbon_direction_linetypes,
    ribbon_direction_alpha = ribbon_direction_alpha,
    ribbon_ctrl_point = ribbon_ctrl_pt,
    region_data = region_data,
    region_fill = region_fill,
    region_color = region_color,
    region_alpha = region_alpha,
    region_width = region_width,
    region_offset = region_offset,
    region_side = region_side,
    ribbon_highlight_rows = ribbon_highlight_rows,
    gene_data = gene_data_override %||% data_list$gene_data,
    geneGap = geneGap, geneWidth = geneWidth,
    geneLabelRadialOffset = geneLabelRadialOffset,
    geneLabelCircumOffset = geneLabelCircumOffset,
    geneLabelCircumLimit = geneLabelCircumLimit,
    geneLabelRotation = geneLabelRotation,
    gene_label_show = gene_ls, gene_label_size = gene_lsz,
    gene_label_wrap = gene_lwrap,
    gene_label_repel_layer = gene_lrepel_layer,
    gene_label_repel_max_overlaps = gene_lrepel_maxov,
    gene_label_repel_box_padding = gene_lrepel_box,
    gene_label_repel_point_padding = gene_lrepel_pt,
    gene_label_repel_min_segment_length = gene_lrepel_minseg,
    gene_label_repel_force = gene_lrepel_force,
    gene_label_repel_seed = gene_lrepel_seed,
    gene_label_orientation = gene_lrepel_orient,
    gene_label_segment = gene_lrepel_seg,
    gene_label_side = gene_lrepel_side,
    gene_label_segment_linetype = gene_lrepel_ltype,
    gene_color_scheme = gene_cs, gene_colors = gene_cols,
    gene_order = gene_ord,
    seq_label_text = seq_label_text,
    seq_label_radius = seq_label_radius,
    seq_label_rotation = seq_label_rotation,
    seq_label_size = seq_label_size,
    seq_label_orientation = seq_label_orientation,
    seq_label_hjust = seq_label_hjust,
    seq_label_vjust = seq_label_vjust,
    axisGap = axisGap, axisMaj = axisMaj, axisMajLen = axisMajLen,
    axisMin = axisMin, axisMinLen = axisMinLen,
    labelSize = labelSize, labelOffset = labelOffset,
    axisLabelOrientation = axisLabelOri,
    axis_label_hide_overlaps = axisLabelHide,
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
              gene_text = integer(0), gene_text_repel = integer(0),
              gene_label_segment = integer(0),
              axis_line = integer(0), axis_seg = integer(0),
              axis_text = integer(0), seq_label = integer(0),
              seq_group_label = integer(0),
              seq_region = integer(0), ribbon_highlight = integer(0))
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
      ggchord_stop("legend_position must be one of 'left', 'right', 'top', 'bottom' or 'inside'")
    }
    if (!is.null(pp$type) && pp$type %in% names(pos)) pos[[pp$type]] <- lp
  }
  pos
}

#' Read the ribbon layer's legend key width/height overrides if set
#' @keywords internal
ggchord_ribbon_key_dims <- function(plot) {
  for (lyr in plot$layers) {
    pp <- lyr$ggchord_params
    if (!is.null(pp) && identical(pp$type, "ribbon")) {
      return(list(width = pp$legend_key_width, height = pp$legend_key_height))
    }
  }
  list(width = NULL, height = NULL)
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
                                positions = list(), legend_key_width = NULL,
                                legend_key_height = NULL) {
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
    if (layout$ribbon_color_scheme %in% c("pident", "value")) {
      # The colorbar follows its effective legend position: vertical and
      # filling the available height at the left/right, horizontal (with a
      # fixed size) at the top/bottom. A "null" key height collapses to zero
      # inside horizontal legend boxes, so a fixed size is used whenever the
      # legend sits at the top/bottom or the box is horizontal.
      ribbon_pos <- positions$ribbon %||% legend_position %||% "right"
      horizontal_legend <- ribbon_pos %in% c("top", "bottom") ||
        identical(legend_box, "horizontal")
      # legend_key_width / legend_key_height control the colourbar key
      # dimensions directly; numbers are interpreted as centimetres.
      key_width <- legend_key_width
      if (!is.null(key_width) && !is.unit(key_width)) key_width <- unit(key_width, "cm")
      key_height <- legend_key_height
      if (!is.null(key_height) && !is.unit(key_height)) key_height <- unit(key_height, "cm")
      value_scheme <- identical(layout$ribbon_color_scheme, "value")
      ribbon_name <- if (value_scheme) {
        layout$ribbon_color_name %||% "value"
      } else {
        "Identity(%)"
      }
      ribbon_limits <- if (value_scheme) {
        lims <- layout$ribbon_color_limits %||% range(layout$ribbon_polys$value, na.rm = TRUE)
        if (lims[1] == lims[2]) lims <- lims + c(-0.5, 0.5)
        lims
      } else {
        c(0, 100)
      }
      ribbon_breaks <- if (value_scheme) {
        layout$ribbon_color_breaks %||% pretty(ribbon_limits, n = 5)
      } else {
        c(0, 50, 80, 90, 95, 100)
      }
      ribbon_fill_scale <- scale_fill_stepsn(
        name    = ribbon_name,
        colours = layout$ribbon_colors,
        limits  = ribbon_limits,
        breaks  = ribbon_breaks,
        guide   = guide_colorbar(
          position = positions$ribbon %||% NULL,
          theme = theme(
            legend.title.position = "top",
            legend.key.height = if (horizontal_legend) {
              key_height %||% unit(1.5, "cm")
            } else {
              key_height %||% unit(1, "null")
            },
            # A horizontal colorbar needs a longer key; the vertical bar keeps
            # the default key width.
            legend.key.width = if (horizontal_legend) {
              key_width %||% unit(4, "cm")
            } else {
              key_width %||% NULL
            }
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

  # The ribbon layer always uses the internal "zfill" aesthetic (never the
  # plain "fill" aesthetic), so that its mapping stays consistent with the
  # ribbon geom's default aesthetics (which are renamed to avoid injecting a
  # plain "fill" default into the ribbon data). This also keeps ggplot2's
  # guide matching working when no gene layer is present, so the Identity(%)
  # colourbar legend is shown even without gene data.
  ribbon_aes <- "fill"
  if (!is.null(ribbon_fill_scale)) {
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
    if (!is.null(gene_fill_scale)) {
      scales[[length(scales) + 1]] <- gene_fill_scale
    }
  } else if (!is.null(gene_fill_scale)) {
    scales[[length(scales) + 1]] <- gene_fill_scale
  }

  # Ribbon alpha is a preset value; use an identity scale so it renders as specified
  if (!is.null(layout$ribbon_polys)) {
    scales[[length(scales) + 1]] <- scale_alpha_identity()
  }
  # Optional per-ribbon outline / linetype mappings use internal aesthetics so
  # they do not collide with the sequence colour scale or the ribbon fill scale.
  if (isTRUE(layout$ribbon_use_outline)) {
    scales[[length(scales) + 1]] <- scale_colour_identity(aesthetics = "zoutline")
  }
  if (isTRUE(layout$ribbon_use_linetype)) {
    scales[[length(scales) + 1]] <- scale_linetype_identity(aesthetics = "zlinetype")
  }
  # Sequence-region bands use their own internal fill aesthetic.
  if (!is.null(layout$region_polys) && nrow(layout$region_polys) > 0) {
    scales[[length(scales) + 1]] <- scale_fill_identity(aesthetics = "zregionfill")
  }
  # Sequence-group label colours use an internal aesthetic so they do not
  # collide with the Seq ID colour scale used by geom_seq().
  if (!is.null(layout$group_labels) && nrow(layout$group_labels) > 0) {
    scales[[length(scales) + 1]] <- scale_colour_identity(aesthetics = "zcolour")
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
  lim <- ggchord_adaptive_limits(layout)
  plot$coordinates <- coord_fixed(
    ratio = 1,
    xlim  = lim$xlim,
    ylim  = lim$ylim,
    clip  = "off"
  )
  plot
}

#' Compute coordinate limits that fit the rendered text boxes
#'
#' The chord geometry is placed in a square, fixed-aspect panel.  Instead of
#' adding one global text-width pad on every side, this helper measures the
#' actual gene/sequence/group/axis label boxes and expands only the sides that
#' need it.  The result is a tighter plot that uses the available panel area.
#' @keywords internal
ggchord_adaptive_limits <- function(layout) {
  ext <- layout$extremes
  if (is.null(ext) || !all(is.finite(c(ext$x_min, ext$x_max,
                                        ext$y_min, ext$y_max)))) {
    return(list(xlim = c(-1, 1), ylim = c(-1, 1)))
  }

  x_range <- ext$x_max - ext$x_min
  y_range <- ext$y_max - ext$y_min
  span <- max(x_range, y_range, 1)
  units_per_inch <- span / 6

  x_lim <- c(ext$x_min, ext$x_max)
  y_lim <- c(ext$y_min, ext$y_max)

  add_boxes <- function(b) {
    if (is.null(b) || nrow(b) == 0) return(invisible(NULL))
    x_lim <<- range(c(x_lim, b$xmin, b$xmax), na.rm = TRUE)
    y_lim <<- range(c(y_lim, b$ymin, b$ymax), na.rm = TRUE)
    invisible(NULL)
  }

  if (nrow(layout$gene_labels) > 0) {
    add_boxes(ggchord_text_boxes(
      layout$gene_labels,
      x_col = "text_x", y_col = "text_y", text_col = "text",
      angle_col = "text_angle", size_col = "size",
      hjust_col = "hjust", vjust_col = "vjust",
      units_per_inch = units_per_inch, box_padding = 0.03
    ))
  }
  if (nrow(layout$seq_labels_df) > 0) {
    add_boxes(ggchord_text_boxes(
      layout$seq_labels_df,
      x_col = "text_x", y_col = "text_y", text_col = "label",
      angle_col = "text_angle", size_col = "size",
      hjust_col = "hjust", vjust_col = "vjust",
      units_per_inch = units_per_inch, box_padding = 0.03
    ))
  }
  if (!is.null(layout$group_labels) && nrow(layout$group_labels) > 0) {
    add_boxes(ggchord_text_boxes(
      layout$group_labels,
      x_col = "text_x", y_col = "text_y", text_col = "label",
      angle_col = "text_angle", size_col = "size",
      hjust_col = "hjust", vjust_col = "vjust",
      units_per_inch = units_per_inch, box_padding = 0.03
    ))
  }
  if (isTRUE(layout$show_axis) && nrow(layout$axis_ticks) > 0) {
    axis_labels <- layout$axis_ticks[!is.na(layout$axis_ticks$label), ,
                                     drop = FALSE]
    if (nrow(axis_labels) > 0) {
      add_boxes(ggchord_text_boxes(
        axis_labels,
        x_col = "label_x", y_col = "label_y", text_col = "label",
        angle_col = "label_angle", size_col = "size",
        hjust_col = "label_hjust", vjust_col = "label_vjust",
        units_per_inch = units_per_inch, box_padding = 0.03
      ))
    }
  }

  # coord_fixed() needs a square data range; centering the square on the
  # rendered content avoids any extra blank strip on one side.
  x_mid <- mean(x_lim)
  y_mid <- mean(y_lim)
  half <- max(diff(range(x_lim)), diff(range(y_lim))) / 2
  half <- half + 0.02 * max(half * 2, 1)

  list(
    xlim = c(x_mid - half, x_mid + half),
    ylim = c(y_mid - half, y_mid + half)
  )
}

#' Estimate the coordinate margin (in data units) needed so that the text
#' labels rendered by the gene/sequence label layers stay inside the figure.
#'
#' Kept for backwards compatibility.  Plot limits are now computed adaptively
#' by \code{\link{ggchord_adaptive_limits}()}, which fits the actual text
#' boxes rather than adding a single conservative margin on every side.
#' @keywords internal
ggchord_label_pad <- function(layout) {
  0
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
                            has_gene = length(cls$gene_poly) > 0,
                            legend_position = plot$theme$legend.position,
                            legend_box = plot$theme$legend.box,
                            positions = ggchord_legend_positions(plot),
                            legend_key_width = ggchord_ribbon_key_dims(plot)$width,
                            legend_key_height = ggchord_ribbon_key_dims(plot)$height)
  plot <- rename_ribbon_layers(plot, cls$ribbon, sc$ribbon_aes, layout)
  plot <- attach_ggchord_scales(plot, sc$scales)
  plot <- set_ggchord_coord(plot, layout)
  plot
}

#' @export
ggplot_build.ggchord <- function(plot, ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  chord <- plot$ggchord
  if (is.null(chord)) {
    ggchord_stop("Not a valid ggchord object: no data stored on the plot. ",
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
  gene_text_repel_indices <- integer(0)
  gene_label_segment_indices <- integer(0)
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
      gene_text_repel = gene_text_repel_indices <- c(gene_text_repel_indices, i),
      gene_label_segment = gene_label_segment_indices <- c(gene_label_segment_indices, i),
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

  # Sequence-group labels are appended at build time (geom_seq() itself keeps
  # returning a single layer for backward compatibility).
  group_label_layer <- ggchord_group_label_layer(layout$group_labels)
  if (!is.null(group_label_layer)) {
    plot$layers[[length(plot$layers) + 1L]] <- group_label_layer
  }


  # ====================================================================
  # Step 4: build and attach scales
  # ====================================================================
  sc <- make_ggchord_scales(layout,
                            has_seq = length(seq_indices) > 0,
                            has_gene = length(gene_poly_indices) > 0,
                            legend_position = plot$theme$legend.position,
                            legend_box = plot$theme$legend.box,
                            positions = ggchord_legend_positions(plot),
                            legend_key_width = ggchord_ribbon_key_dims(plot)$width,
                            legend_key_height = ggchord_ribbon_key_dims(plot)$height)
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
