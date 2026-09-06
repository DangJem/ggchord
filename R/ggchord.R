# ggchord.R - constructor
# Entry point of the layered ggplot2 API for chord diagrams
# v0.6.0: plot objects are self-contained. Data and parameters are stored on
# the plot object itself, and the layout is computed at build time so that
# print(), ggsave(), ggplot_build() and other ggplot2 workflows all work.

# Global variable declarations (to avoid R CMD check NOTEs)
globalVariables(c(
  "x", "y", "group", "pident", "fill", "colour", "strand", "anno", "accver",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "xend", "yend", ".component",
  "label", "label_x", "label_y", "size", "angle",
  "fill_col", "alpha", "label_hjust", "label_vjust", "label_angle",
  "linetype", "zcolour", "zregionfill", "zoutline", "zlinetype",
  "outline_col", "linetype_val", "value", "source_row", "direction",
  "seq_colour", "ribbon_fill", "ribbon_alpha",
  "ribbon_colour", "ribbon_linetype", "gene_fill", "feature_fill",
  "feature_shape", "seq_ring", "bundle_n", "bundle_weight", "density",
  "region_fill", ".bundle_n", ".bundle_weight", ".bundle_density"
))

#' ggchord: layered multi-sequence alignment chord diagrams for ggplot2
#'
#' ggchord visualizes multi-sequence alignment results using ggplot2's layered grammar.
#' The \code{ggchord()} constructor handles data validation and global settings;
#' the \code{geom_*} layers are stacked as needed, each responsible for its own layout parameters and visual rendering.
#' The layout is computed lazily when the plot is built (e.g. via \code{print()},
#' \code{ggsave()}, or \code{ggplot_build()}).
#'
#' @param seq_data data.frame/tibble with accver and length.
#' @param ribbon_data Optional alignment table with qaccver, saccver, qstart,
#'   qend, sstart and send. length and pident are optional unless requested
#'   by a mapping, filter or statistic.
#' @param gene_data Optional table with accver, start, end and strand.
#'   Annotation (anno) may be omitted.
#' @param debug Logical. Whether to output debug information, default FALSE
#' @param validate Character, default \code{"warn"}. How to run the structured
#'   input-data validation (see \code{\link{validate_ggchord_data}}):
#'   \code{"warn"} emits a single summary warning when the data has problems
#'   and caches the full report on the plot object
#'   (\code{p$ggchord$validation}); \code{"error"} stops on severe problems;
#'   \code{"none"} skips the diagnostic validation (the cheap structural
#'   checks that prevent crashes are still performed).
#' @param ... Reserved for clear migration errors from removed constructor
#'   arguments. Use [labs()], [coord_chord()] and [theme()] for presentation.
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
#'   geom_link_ribbon() +
#'   geom_gene()
#' print(p)
#'
#' @importFrom ggplot2 ggplot_build
ggchord <- function(
    seq_data,
    ribbon_data = NULL,
    gene_data = NULL,
    debug = FALSE,
    validate = c("warn", "error", "none"),
    ...
) {
  seq_data <- ggchord_normalize_accver(seq_data)
  gene_data <- ggchord_normalize_accver(gene_data)

  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  ggchord_reject_retired(dots, "ggchord()", c(
    title = "labs(title = ...)",
    rotation = "coord_chord(rotation = ...)",
    panel_margin = "theme(plot.margin = margin(...))",
    show_legend = "theme(legend.position = 'none')"
  ))
  if (length(dots)) {
    ggchord_stop("ggchord(): unused argument(s): ",
                 paste(names(dots), collapse = ", "))
  }

  validate <- match.arg(validate)
  # ====================================================================
  # 1. Validate data
  # ====================================================================
  if (!is.logical(debug) || length(debug) != 1 || is.na(debug)) {
    ggchord_stop("debug must be TRUE or FALSE")
  }
  # One structured validation pass supplies both the always-on crash-prevention
  # checks and the optional diagnostics. `validate = "none"` skips the costly
  # coordinate, duplicate and self-link checks but still rejects unsafe input.
  validation_result <- validate_ggchord_data(
    seq_data, ribbon_data, gene_data, strict = FALSE,
    check_coordinates = validate != "none",
    check_duplicates = validate != "none",
    check_self_links = validate != "none"
  )
  structural_errors <- ggchord_structural_validation_errors(validation_result)
  if (nrow(structural_errors) > 0) {
    ggchord_stop("ggchord(): ", structural_errors$message[1])
  }

  if (!is.null(ribbon_data)) {
    if (nrow(ribbon_data) == 0) warning("No valid alignment data in ribbon_data")
    if (debug) cat("Number of alignment data rows: ", nrow(ribbon_data), "\n")
  }
  if (!is.null(gene_data)) {
    if (nrow(gene_data) == 0) warning("No valid gene annotation data in gene_data")
    if (debug) cat("Number of gene annotation rows: ", nrow(gene_data), "\n")
  }

  validation <- if (validate == "none") NULL else validation_result
  if (!is.null(validation)) {
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
  p <- ggplot2::ggplot() +
    coord_chord(rotation = 45) +
    theme_ggchord()

  p$ggchord <- list(
    data   = list(seq_data = seq_data, ribbon_data = ribbon_data,
                  gene_data = gene_data),
    global = list(rotation = 45, debug = debug, validate = validate),
    validation = validation,
    # Shared reference environment for plot-owned layout caching.
    ref    = new.env(),
    layout = NULL
  )
  class(p) <- c("ggchord", class(p))
  p
}

# ====================================================================
# +.ggchord method
# ====================================================================

#' Combine a ggchord plot with ggplot2 objects
#'
#' Uses ggplot2's standard composition semantics for layers, scales, themes,
#' coordinates and annotations.
#'
#' @param e1 A ggchord object
#' @param e2 A ggplot2 component.
#' @return A ggchord object
#' @export
`+.ggchord` <- function(e1, e2) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  old_settings <- attr(e1$theme, "ggchord.settings")
  new_settings <- if (inherits(e2, "theme")) {
    attr(e2, "ggchord.settings")
  } else NULL
  if (inherits(e2, "theme") && is.null(new_settings) &&
      !is.null(old_settings)) {
    # Preserve chord settings across an ordinary theme(), while recording
    # explicit generic legend fields so role guides inherit those values
    # instead of mistaking ggplot2's built-in line-unit defaults for user
    # choices.
    new_settings <- old_settings
    for (field in names(new_settings$legend)) {
      value <- e2[[paste0("legend.", field)]]
      if (!is.null(value)) new_settings$legend[[field]] <- value
    }
  }
  if (inherits(e2, "Scale")) {
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
  # Invalidate the plot-owned layout after adding any component.
  if (!is.null(e1$ggchord$ref)) e1$ggchord$ref$layout <- NULL
  attr(p$theme, "ggchord.settings") <- new_settings %||% old_settings
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

compute_chord_geometry_single <- function(plot, geometry_cache = NULL) {
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
  seq_data_override <- NULL
  ribbon_data_override <- NULL
  gene_data_override <- NULL
  region_data_override <- NULL
  gene_label_params <- list()
  gene_repel_params <- list()
  axis_params   <- list()
  seq_label_params <- list()
  seq_region_params <- list()
  seq_layer_requested <- FALSE
  seq_ring_mapped <- FALSE
  ribbon_layer_requested <- FALSE
  gene_geometry_layer_requested <- FALSE
  gene_label_layer <- FALSE
  gene_repel_layer <- FALSE
  feature_stack_position <- NULL

  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    pp <- lyr$ggchord_params
    if (is.null(pp)) next
    switch(pp$type,
      seq               = {
        seq_params <- pp
        seq_layer_requested <- TRUE
        seq_ring_mapped <- "seq_ring" %in% names(lyr$ggchord_input_mapping)
        seq_data_override <- ggchord_resolve_layer_input(
          lyr, data_list$seq_data
        )
      },
      ribbon            = {
        ribbon_params <- pp
        ribbon_layer_requested <- TRUE
        ribbon_data_override <- ggchord_resolve_layer_input(
          lyr, data_list$ribbon_data
        )
      },
      gene              = {
        if (isTRUE(lyr$position$ggchord_feature_stack)) {
          pp$feature_stack_position <- lyr$position
          feature_stack_position <- lyr$position
        }
        gene_params <- pp
        gene_geometry_layer_requested <- TRUE
        gene_data_override <- pp$gene_data_override %||%
          ggchord_resolve_layer_input(lyr, data_list$gene_data)

      },
      gene_label        = {
        if (isTRUE(lyr$position$ggchord_feature_stack)) {
          feature_stack_position <- lyr$position
        }
        gene_label_params <- pp
        gene_label_layer <- TRUE
        label_data <- ggchord_resolve_layer_input(lyr, data_list$gene_data)
        if (!is.null(lyr$ggchord_input_data) ||
            length(intersect(names(lyr$ggchord_input_mapping),
                             lyr$ggchord_role_aes)) > 0) {
          gene_data_override <- label_data
        }
      },
      gene_label_repel  = {
        if (isTRUE(lyr$position$ggchord_feature_stack)) {
          feature_stack_position <- lyr$position
        }
        gene_repel_params <- pp
        gene_repel_layer <- TRUE
        label_data <- ggchord_resolve_layer_input(lyr, data_list$gene_data)
        if (!is.null(lyr$ggchord_input_data) ||
            length(intersect(names(lyr$ggchord_input_mapping),
                             lyr$ggchord_role_aes)) > 0) {
          gene_data_override <- label_data
        }
      },
      seq_label         = seq_label_params <- pp,
      seq_region        = {
        seq_region_params <- pp
        region_data_override <- ggchord_resolve_layer_input(lyr, pp$regions)
      }

    )
  }

  # --- Process sequences ---
  seq_data <- seq_data_override %||% data_list$seq_data
  seqs     <- seq_data$accver
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

  # Rings are explicit input roles. Their scale values are radii, so no ring
  # count or spacing is guessed from the data. This deliberately keeps the
  # existing seq_radius interface unchanged for plots without a ring mapping.
  seq_ring <- NULL
  if (isTRUE(seq_ring_mapped)) {
    if (!is.null(seq_params$seq_radius)) {
      ggchord_stop(
        "geom_seq(): `seq_radius` cannot be combined with a `seq_ring` ",
        "mapping; set radii with scale_seq_ring_manual(values = ...)"
      )
    }
    ring_scale <- plot$scales$get_scales("seq_ring")
    if (is.null(ring_scale)) {
      ggchord_stop(
        "A mapped `seq_ring` requires scale_seq_ring_manual(values = ...)"
      )
    }
    raw_ring <- as.character(seq_data$seq_ring[match(seqs, seq_data$accver)])
    if (anyNA(raw_ring) || any(!nzchar(raw_ring))) {
      ggchord_stop("Mapped `seq_ring` values must be non-missing")
    }
    trained_ring_scale <- ring_scale$clone()
    trained_ring_scale$train(raw_ring)
    mapped_radius <- suppressWarnings(as.numeric(
      trained_ring_scale$map(raw_ring)
    ))
    if (anyNA(mapped_radius) || any(!is.finite(mapped_radius)) ||
        any(mapped_radius <= 0)) {
      ggchord_stop(
        "scale_seq_ring_manual(): every used ring must map to a finite ",
        "positive radius"
      )
    }
    seq_ring <- stats::setNames(raw_ring, seqs)
    seqRadius <- stats::setNames(mapped_radius, seqs)
  }

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

  # --- Process ribbons ---
  ribbon_gap_auto <- is.null(ribbon_params$ribbon_gap)
  ribbonGap  <- process_sequence_param(ribbon_params$ribbon_gap %||% 0.15,
                                       seqs, "ribbon_gap", 0.15)
  ribbon_color_scheme <- ribbon_params$ribbon_color_scheme %||% "pident"
  ribbon_color_by     <- ribbon_params$ribbon_color_by
  ribbon_color_limits <- ribbon_params$ribbon_color_limits
  ribbon_color_breaks <- ribbon_params$ribbon_color_breaks
  ribbon_color_name   <- ribbon_params$ribbon_color_name
  ribbon_alpha    <- ribbon_params$ribbon_alpha %||% 0.42
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
  ribbon_data <- if (ribbon_layer_requested) {
    ribbon_data_override %||%
      data_list$ribbon_data
  } else {
    NULL
  }
  ribbon_stat_report <- NULL
  ribbon_stat_data <- NULL
  ribbon_stat <- ribbon_params$ribbon_stat
  if (!is.null(ribbon_stat) && !is.null(ribbon_data)) {
    if (identical(ribbon_stat$type, "bundle")) {
      computed <- bundle_ggchord_ribbons(
        ribbon_data = ribbon_data,
        seq_data = seq_data,
        bins = ribbon_stat$bins,
        min_bundle = ribbon_stat$min_bundle,
        weight = ribbon_stat$weight,
        group_by = ribbon_stat$group_by
      )
      computed$data$bundle_n <- computed$data$.bundle_n
      computed$data$bundle_weight <- computed$data$.bundle_weight
      computed$data$density <- computed$data$.bundle_density
    } else if (identical(ribbon_stat$type, "density")) {
      computed <- ggchord_ribbon_density(
        ribbon_data = ribbon_data,
        seq_data = seq_data,
        bins = ribbon_stat$bins,
        weight = ribbon_stat$weight,
        group_by = ribbon_stat$group_by,
        caller = "stat_ribbon_density()"
      )
    } else {
      ggchord_stop("Unknown ggchord ribbon stat: ", ribbon_stat$type)
    }
    ribbon_data <- computed$data
    ribbon_stat_data <- computed$data
    ribbon_stat_report <- computed$report
  }
  has_ribbon_data <- !is.null(ribbon_data) && nrow(ribbon_data) > 0

  if (has_ribbon_data) {
    if (is.null(ribbon_colors)) {
      ribbon_colors <- switch(ribbon_color_scheme,
        single = "steelblue",
        query  = {
          mix <- 0.5
          sapply(seq_colors, function(col) {
            cols <- grDevices::col2rgb(col)
            light_cols <- cols + (255 - cols) * mix
            grDevices::rgb(light_cols[1,], light_cols[2,], light_cols[3,],
                maxColorValue = 255)
          })
        },
        subject = {
          mix <- 0.5
          sapply(seq_colors, function(col) {
            cols <- grDevices::col2rgb(col)
            light_cols <- cols + (255 - cols) * mix
            grDevices::rgb(light_cols[1,], light_cols[2,], light_cols[3,],
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
      if (!ribbon_color_by %in% colnames(ribbon_data)) {
        ggchord_stop("ribbon_color_by column '", ribbon_color_by, "' not found in ribbon_data")
      }
      if (!is.numeric(ribbon_data[[ribbon_color_by]]) ||
          any(!is.finite(ribbon_data[[ribbon_color_by]]))) {
        ggchord_stop("ribbon_color_by column '", ribbon_color_by, "' must be numeric and finite")
      }
    }
    if (!is.null(ribbon_alpha_by)) {
      if (!ribbon_alpha_by %in% colnames(ribbon_data)) {
        ggchord_stop("ribbon_alpha_by column '", ribbon_alpha_by, "' not found in ribbon_data")
      }
      if (!is.numeric(ribbon_data[[ribbon_alpha_by]]) ||
          any(!is.finite(ribbon_data[[ribbon_alpha_by]]))) {
        ggchord_stop("ribbon_alpha_by column '", ribbon_alpha_by, "' must be numeric and finite")
      }
    }
    if (!is.null(ribbon_outline_by) &&
        !ribbon_outline_by %in% colnames(ribbon_data)) {
      ggchord_stop("ribbon_outline_by column '", ribbon_outline_by, "' not found in ribbon_data")
    }
    if (!is.null(ribbon_linetype_by) &&
        !ribbon_linetype_by %in% colnames(ribbon_data)) {
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
    gene_params$gene_label_size %||%
    ggchord_theme_text_size(plot, "ggchord.gene.label", 2.5)
  # Repelled labels now use mode-owned deterministic positioning. Manual
  # rotation and offsets remain available through geom_gene_label(), but are
  # intentionally not inherited by geom_gene_label_repel().
  if (gene_repel_layer) {
    gene_lr <- gene_lro <- gene_lco <- 0
    gene_lcl <- TRUE
  } else {
    gene_lr <- lbl$gene_label_rotation %||%
      gene_params$gene_label_rotation %||% 0
    gene_lro <- lbl$gene_label_radial_offset %||%
      gene_params$gene_label_radial_offset %||% 0
    gene_lco <- lbl$gene_label_circum_offset %||%
      gene_params$gene_label_circum_offset %||% 0
    gene_lcl <- if (is.null(lbl$gene_label_circum_limit)) {
      if (is.null(gene_params$gene_label_circum_limit)) TRUE
      else gene_params$gene_label_circum_limit
    } else lbl$gene_label_circum_limit
  }
  gene_lwrap  <- lbl$gene_label_wrap %||% gene_params$gene_label_wrap
  gene_lorientation <- if (gene_repel_layer) {
    "radial"
  } else {
    lbl$gene_label_orientation %||% "horizontal"
  }
  gene_loverlap <- if (gene_repel_layer) {
    "allow"
  } else {
    lbl$gene_label_overlap %||% "hide"
  }
  gene_lrepel_layer <- gene_repel_layer
  gene_lrepel_maxov <- gene_repel_params$max_overlaps %||% Inf
  gene_lrepel_layout <- gene_repel_params$gene_label_layout %||% "radial"
  gene_lrepel_fit    <- gene_repel_params$gene_label_fit %||% "wrap"
  gene_lrepel_lines  <- gene_repel_params$gene_label_max_lines %||% 2L
  gene_lrepel_side   <- if (gene_repel_layer) {
    gene_repel_params$gene_label_side %||% "outside"
  } else {
    lbl$gene_label_side %||% "outside"
  }
  gene_lrepel_segment_overlap <-
    gene_repel_params$gene_label_segment_overlap %||% "fade"
  gene_lrepel_segment_overlap_alpha <-
    gene_repel_params$gene_label_segment_overlap_alpha %||% 0.18
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
  if (!is.numeric(gene_lrepel_maxov) || length(gene_lrepel_maxov) != 1 ||
      is.na(gene_lrepel_maxov) || gene_lrepel_maxov < 0) {
    ggchord_stop("max_overlaps must be a non-negative number or Inf")
  }

  geneGap    <- process_gene_param(gene_off, seqs, "gene_offset", 0.1, FALSE)
  geneWidth  <- process_gene_param(gene_w, seqs, "gene_width", 0.05, FALSE)

  # Feature geometry is resolved before coordinate generation because these
  # values change the actual polygon, not only its appearance. A user-supplied
  # feature-shape scale is cloned and trained here so the layout and legend use
  # exactly the same category-to-shape mapping.
  gene_data_layout <- if (gene_geometry_layer_requested || gene_label_layer ||
      gene_repel_layer) {
    gene_data_override %||% data_list$gene_data
  } else {
    NULL
  }
  ribbon_obstacles <- plot$ggchord$obstacles %||% ggchord_empty_obstacles()
  if (!is.null(gene_data_layout) && !"anno" %in% names(gene_data_layout)) {
    gene_data_layout$anno <- rep(NA_character_, nrow(gene_data_layout))
  }
  feature_shape_pal <- NULL
  feature_shape_order <- NULL
  if (isTRUE(gene_params$is_feature) && !is.null(gene_data_layout) &&
      nrow(gene_data_layout) > 0L) {
    raw_shape <- as.character(
      gene_data_layout$.feature_shape_raw %||%
        rep(gene_params$feature_shape %||% "arrow", nrow(gene_data_layout))
    )
    feature_shape_order <- unique(raw_shape)
    if (isTRUE(gene_params$feature_shape_mapped)) {
      shape_scale <- plot$scales$get_scales("feature_shape")
      if (!is.null(shape_scale)) {
        shape_scale <- shape_scale$clone()
        shape_scale$train(raw_shape)
        mapped_shape <- as.character(shape_scale$map(raw_shape))
      } else {
        allowed_shape <- c("arrow", "block", "chevron", "lollipop")
        if (all(feature_shape_order %in% allowed_shape)) {
          feature_shape_pal <- stats::setNames(
            feature_shape_order, feature_shape_order
          )
        } else {
          feature_shape_pal <- stats::setNames(
            rep(allowed_shape, length.out = length(feature_shape_order)),
            feature_shape_order
          )
        }
        mapped_shape <- unname(feature_shape_pal[raw_shape])
      }
    } else {
      mapped_shape <- raw_shape
    }
    allowed_shape <- c("arrow", "block", "chevron", "lollipop")
    mapped_shape[is.na(mapped_shape)] <- "arrow"
    if (any(!mapped_shape %in% allowed_shape)) {
      ggchord_stop(
        "feature_shape scale values must use 'arrow', 'block', ",
        "'chevron', or 'lollipop'"
      )
    }
    if (is.null(feature_shape_pal) &&
        isTRUE(gene_params$feature_shape_mapped)) {
      feature_shape_pal <- stats::setNames(
        mapped_shape[match(feature_shape_order, raw_shape)],
        feature_shape_order
      )
    }
    gene_data_layout <- as.data.frame(gene_data_layout, stringsAsFactors = FALSE)
    gene_data_layout$.feature_shape <- mapped_shape
  }
  if (!is.null(feature_stack_position) && !is.null(gene_data_layout)) {
    gene_data_layout <- ggchord_stack_feature_tracks(
      gene_data_layout, feature_stack_position
    )
  }
  geneLabelRadialOffset <- process_gene_param(gene_lro, seqs,
                                              "gene_label_radial_offset", 0, FALSE)
  geneLabelCircumOffset <- process_gene_param(gene_lco, seqs,
                                              "gene_label_circum_offset", 0, FALSE)
  geneLabelCircumLimit  <- process_gene_param(gene_lcl, seqs,
                                              "gene_label_circum_limit", TRUE, TRUE)
  geneLabelRotation     <- process_gene_param(gene_lr, seqs,
                                              "gene_label_rotation", 0, FALSE)

  # --- Process axes ---
  axis_theme <- ggchord_plot_settings(plot)$axis
  show_axis <- !isTRUE(axis_theme$hidden)
  # Pass physical distances as inches. compute_chord_layout() converts them
  # against the actual curved sequence span once that geometry is available.
  axis_unit_data <- ggchord_unit_inches
  axisGap    <- process_sequence_param(axis_unit_data(axis_theme$gap),
                                       seqs, "axis.gap", 0.04)
  axisMaj    <- process_sequence_param(axis_params$axis_tick_major_number %||% 3,
                                       seqs, "axis_tick_major_number", 3)
  axisMajLen <- process_sequence_param(axis_unit_data(axis_theme$ticks.length),
                                       seqs, "axis.ticks.length", 0.02)
  axisMin    <- process_sequence_param(axis_params$axis_tick_minor_number %||% 4,
                                       seqs, "axis_tick_minor_number", 4)
  axisMinLen <- process_sequence_param(
    axis_unit_data(axis_theme$minor.ticks.length),
    seqs, "axis.minor.ticks.length", 0.01
  )
  axis_theme_size <- ggchord_theme_text_size(
    plot, "ggchord.axis.text", 3
  )
  labelSize  <- process_sequence_param(
    axis_theme_size,
    seqs, "axis_label_size", axis_theme_size
  )
  labelOffset <- process_sequence_param(axis_unit_data(axis_theme$text.offset),
                                        seqs, "axis.text.offset", 0.02)
  axisLabelHide <- isTRUE(axis_theme$text.check.overlap)
  axisLabelOri <- process_axis_orientation(
    axis_theme$text.orientation, seqs
  )
  if (!is.logical(show_axis) || length(show_axis) != 1 || is.na(show_axis)) {
    ggchord_stop("show_axis must be TRUE or FALSE")
  }
  axis_numeric <- list(
    axis_gap = axisGap,
    axis_tick_major_length = axisMajLen,
    axis_tick_minor_length = axisMinLen,
    axis_label_size = labelSize,
    axis_text_offset = labelOffset
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
  axis_breaks <- axis_minor_breaks <- axis_labels <- NULL
  position_scale <- plot$scales$get_scales("seq_position")
  if (!is.null(position_scale)) {
    axis_breaks <- axis_minor_breaks <- axis_labels <- setNames(
      vector("list", length(seqs)), seqs
    )
    for (id in seqs) {
      sc <- position_scale$clone()
      sc$train(c(0, lens[[id]]))
      br <- sc$get_breaks()
      br <- br[is.finite(br) & br >= 0 & br <= lens[[id]]]
      axis_breaks[[id]] <- br
      minor <- sc$get_breaks_minor()
      axis_minor_breaks[[id]] <- minor[
        is.finite(minor) & minor >= 0 & minor <= lens[[id]]
      ]
      axis_labels[[id]] <- sc$get_labels(br)
    }
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
      # are matched positionally to the sequences (named by accver).
      process_sequence_param(seq_label_params$seq_labels, seqs, "seq_labels",
                             default_value = seqs)
    }
    seq_label_radius <- process_sequence_param(
      seq_label_params$seq_label_radius, seqs, "seq_label_radius", 1)
    seq_label_rotation <- process_sequence_param(
      seq_label_params$seq_label_rotation, seqs, "seq_label_rotation", 0)
    seq_theme_size <- ggchord_theme_text_size(
      plot, "ggchord.seq.label", 3
    )
    seq_label_size <- process_sequence_param(
      seq_label_params$seq_label_size, seqs, "seq_label_size", seq_theme_size
    )
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
  region_data <- region_data_override %||% seq_region_params$regions
  region_fill   <- seq_region_params$region_fill %||% "#F59E0B"
  region_color  <- seq_region_params$region_color %||% "#B45309"
  region_alpha  <- seq_region_params$region_alpha %||% 0.25
  region_width  <- seq_region_params$region_width %||% 0.08
  region_offset <- seq_region_params$region_offset %||% 0
  region_side   <- seq_region_params$region_side %||% "inside"
  if (!is.null(region_data) && !is.data.frame(region_data)) {
    ggchord_stop("geom_seq_region(): regions must be a data.frame")
  }

  # ====================================================================
  # Step 2: compute the layout
  # ====================================================================
  coord_rotation <- if (isTRUE(plot$coordinates$ggchord_coord)) {
    plot$coordinates$rotation
  } else {
    global$rotation
  }

  layout <- compute_chord_layout(
    seqs = seqs, lens = lens, seq_labels = seq_labels,
    seq_colors = seq_colors, seqRadius = seqRadius,
    seq_curvature = seq_curvature, orientation = orientation,
    seq_gap = seq_gap,
    ribbon_data = ribbon_data, ribbonGap = ribbonGap,
    ribbon_gap_auto = ribbon_gap_auto,
    link_avoid = ribbon_params$link_avoid %||% "none",
    ribbon_obstacles = ribbon_obstacles,
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
    gene_data = gene_data_layout,
    draw_gene_geometry = gene_geometry_layer_requested,
    geneGap = geneGap, geneWidth = geneWidth,
    geneLabelRadialOffset = geneLabelRadialOffset,
    geneLabelCircumOffset = geneLabelCircumOffset,
    geneLabelCircumLimit = geneLabelCircumLimit,
    geneLabelRotation = geneLabelRotation,
    gene_label_show = gene_ls, gene_label_size = gene_lsz,
    gene_label_wrap = gene_lwrap,
    gene_label_fit = gene_lrepel_fit,
    gene_label_max_lines = gene_lrepel_lines,
    gene_label_orientation = gene_lorientation,
    gene_label_overlap = gene_loverlap,
    gene_label_repel_layer = gene_lrepel_layer,
    gene_label_repel_max_overlaps = gene_lrepel_maxov,
    gene_label_layout = gene_lrepel_layout,
    gene_label_side = gene_lrepel_side,
    gene_label_segment_overlap = gene_lrepel_segment_overlap,
    gene_label_segment_overlap_alpha = gene_lrepel_segment_overlap_alpha,
    gene_label_segment_linetype = gene_lrepel_ltype,
    gene_color_scheme = gene_cs, gene_colors = gene_cols,
    gene_order = gene_ord,
    feature_shape_pal = feature_shape_pal,
    feature_shape_order = feature_shape_order,
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
    axis_breaks = axis_breaks,
    axis_minor_breaks = axis_minor_breaks,
    axis_labels = axis_labels,
    axis_label_hide_overlaps = axisLabelHide,
    show_axis = show_axis,
    rotation = coord_rotation, debug = global$debug,
    geometry_cache = geometry_cache
  )

  layout$seq_ring <- seq_ring
  layout$seq_ring_radius <- if (is.null(seq_ring)) NULL else seqRadius
  layout$ribbon_stat_data <- ribbon_stat_data
  layout$ribbon_stat_report <- ribbon_stat_report
  if (!is.null(seq_ring) && length(layout$seq_arcs)) {
    layout$seq_arcs <- lapply(names(layout$seq_arcs), function(id) {
      arc <- layout$seq_arcs[[id]]
      arc$seq_ring <- unname(seq_ring[[id]])
      arc
    }) |>
      stats::setNames(names(layout$seq_arcs))
  }

  layout
}

#' Extract one drawable component from a computed layout
#' @noRd
ggchord_layout_component <- function(layout, type, fallback = data.frame()) {
  switch(type,
    seq = if (length(layout$seq_arcs) > 0) do.call(rbind, layout$seq_arcs) else fallback,
    ribbon = layout$ribbon_polys %||% fallback,
    link = layout$link_lines %||% fallback,
    gene_poly = layout$gene_polys %||% fallback,
    gene_text = layout$gene_labels %||% fallback,
    gene_text_repel = layout$gene_labels %||% fallback,
    gene_label_segment = layout$gene_label_segments %||% fallback,
    gene_label_repel = ggchord_repel_geometry(layout),
    seq_label = layout$seq_labels_df %||% fallback,
    seq_region = layout$region_polys %||% fallback,
    axis_line = layout$axis_lines %||% fallback,
    axis_seg = layout$axis_ticks %||% fallback,
    axis_text = {
      d <- layout$axis_ticks %||% fallback
      if (nrow(d) > 0 && "label" %in% names(d)) d[!is.na(d$label), , drop = FALSE] else d
    },
    axis = ggchord_axis_geometry(layout),
    fallback
  )
}

#' Compute a plot-owned, per-layer geometry registry
#' @noRd
compute_chord_geometry <- function(plot) {
  chord <- plot$ggchord
  if (is.null(chord)) {
    ggchord_stop("Not a valid ggchord object: no data stored on the plot")
  }
  # Layers added through ordinary ggplot2 mechanisms may not have passed the
  # list branch of +.ggchord. Assign deterministic IDs before grouping.
  next_id <- 1L
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (is.null(lyr$ggchord_type)) next
    if (is.null(lyr$ggchord_layer_id)) {
      lyr$ggchord_layer_id <- sprintf("layer-%04d", next_id)
      plot$layers[[i]] <- lyr
    }
    next_id <- next_id + 1L
  }

  geometry_cache <- new.env(parent = emptyenv())
  avoid_requested <- any(vapply(plot$layers, function(layer) {
    params <- layer$ggchord_params
    !is.null(params) && params$type %in% c("ribbon", "link") &&
      !identical(params$link_avoid %||% "none", "none") &&
      is.null(params$ribbon_gap %||% params$link_gap)
  }, logical(1)))
  plot$ggchord$obstacles <- if (avoid_requested) {
    ggchord_collect_obstacles(plot, geometry_cache)
  } else {
    ggchord_empty_obstacles()
  }
  primary <- compute_chord_geometry_single(plot, geometry_cache)
  primary$obstacles <- plot$ggchord$obstacles
  ids <- vapply(plot$layers, function(x) x$ggchord_layer_id %||% "", character(1))
  groups <- split(which(nzchar(ids)), ids[nzchar(ids)])

  group_type <- vapply(groups, function(idx) {
    types <- vapply(idx, function(i) plot$layers[[i]]$ggchord_params$type %||% "",
                    character(1))
    main <- types[types %in% c(
      "seq", "ribbon", "link", "gene", "gene_label", "gene_label_repel", "axis",
      "seq_label", "seq_region"
    )]
    if (length(main)) main[length(main)] else ""
  }, character(1))
  type_counts <- table(group_type[nzchar(group_type)])

  first_group <- function(type) {
    hit <- names(group_type)[group_type == type]
    if (length(hit)) groups[[hit[1]]] else integer(0)
  }
  seq_dep <- first_group("seq")
  gene_dep <- first_group("gene")
  ribbon_dep <- first_group("ribbon")
  gene_geometry_dep <- unlist(
    groups[names(group_type)[group_type == "gene"]], use.names = FALSE
  )

  registry <- list()
  inputs <- list()
  layouts <- list()
  for (id in names(groups)) {
    idx <- groups[[id]]
    main_type <- group_type[[id]]
    needs_own <- nzchar(main_type) && type_counts[[main_type]] > 1
    sub_layout <- primary
    if (isTRUE(needs_own)) {
      deps <- seq_dep
      if (main_type %in% c("gene_label", "gene_label_repel")) {
        deps <- c(deps, gene_dep)
      }
      if (main_type == "ribbon") deps <- c(deps, gene_geometry_dep)
      sub_plot <- plot
      sub_plot$layers <- plot$layers[sort(unique(c(deps, idx)))]
      sub_layout <- compute_chord_geometry_single(sub_plot, geometry_cache)
    }
    if (main_type == "link") {
      link_layer <- plot$layers[[idx[1L]]]
      link_input <- ggchord_resolve_layer_input(link_layer)
      sub_layout$link_lines <- ggchord_attach_input_columns(
        ggchord_link_geometry(link_input, link_layer$ggchord_params, primary), link_input)
    }
    layouts[[id]] <- sub_layout
    registry[[id]] <- list()
    inputs[[id]] <- list()
    for (i in idx) {
      lyr <- plot$layers[[i]]
      component <- lyr$ggchord_type
      registry[[id]][[component]] <- ggchord_layout_component(
        sub_layout, component, lyr$ggchord_placeholder %||% data.frame()
      )
      fallback <- switch(component,
        seq = chord$data$seq_data,
        axis_line = chord$data$seq_data,
        axis_seg = chord$data$seq_data,
        axis_text = chord$data$seq_data,
        axis = chord$data$seq_data,
        seq_label = chord$data$seq_data,
        ribbon = chord$data$ribbon_data,
        gene_poly = chord$data$gene_data,
        gene_text = chord$data$gene_data,
        gene_text_repel = chord$data$gene_data,
        gene_label_segment = chord$data$gene_data,
        gene_label_repel = chord$data$gene_data,
        seq_region = lyr$ggchord_params$regions,
        NULL
      )
      inputs[[id]][[component]] <- if (
          component == "ribbon" && !is.null(sub_layout$ribbon_stat_data)) {
        sub_layout$ribbon_stat_data
      } else {
        ggchord_resolve_layer_input(lyr, fallback)
      }
      registry[[id]][[component]] <- ggchord_attach_input_columns(
        registry[[id]][[component]], inputs[[id]][[component]]
      )
    }
  }
  primary$layer_geometry <- registry
  primary$layer_inputs <- inputs
  primary$layer_layouts <- layouts

  # Plot limits must see all independent layers, not only the compatibility
  # fields in the primary layout.
  collect <- function(component) {
    values <- lapply(registry, `[[`, component)
    values <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, values)
    if (length(values)) ggchord_rbind_fill(values) else data.frame()
  }
  for (pair in list(
    c("ribbon_polys", "ribbon"), c("link_lines", "link"), c("gene_polys", "gene_poly"),
    c("gene_labels", "gene_text"), c("gene_label_segments", "gene_label_segment"),
    c("seq_labels_df", "seq_label"), c("region_polys", "seq_region"),
    c("axis_lines", "axis_line"), c("axis_ticks", "axis_seg")
  )) {
    combined <- collect(pair[2])
    if (nrow(combined) > 0) primary[[pair[1]]] <- combined
  }
  repel_labels <- collect("gene_text_repel")
  fixed_labels <- collect("gene_text")
  all_gene_labels <- Filter(function(x) nrow(x) > 0,
                            list(fixed_labels, repel_labels))
  if (length(all_gene_labels)) {
    primary$gene_labels <- ggchord_rbind_fill(all_gene_labels)
  }

  palettes <- lapply(layouts, function(x) x$gene_pal)
  palettes <- Filter(function(x) !is.null(x) && length(x) > 0, palettes)
  if (length(palettes)) {
    pal <- do.call(c, unname(palettes))
    primary$gene_pal <- pal[!duplicated(names(pal), fromLast = TRUE)]
    orders <- unlist(lapply(layouts, function(x) x$final_gene_order),
                     use.names = FALSE)
    primary$final_gene_order <- unique(orders)
  }
  shape_palettes <- lapply(layouts, function(x) x$feature_shape_pal)
  shape_palettes <- Filter(
    function(x) !is.null(x) && length(x) > 0L, shape_palettes
  )
  if (length(shape_palettes)) {
    shape_pal <- do.call(c, unname(shape_palettes))
    primary$feature_shape_pal <- shape_pal[
      !duplicated(names(shape_pal), fromLast = TRUE)
    ]
    primary$feature_shape_order <- unique(unlist(lapply(
      layouts, function(x) x$feature_shape_order
    ), use.names = FALSE))
  }
  primary$extremes <- get_plot_extremes(
    allRibbon = ggchord_rbind_fill(Filter(Negate(is.null), list(primary$ribbon_polys, primary$link_lines))),
    seqArcs = primary$seq_arcs,
    axisLines = primary$axis_lines,
    axisTicks = primary$axis_ticks,
    gene_polys = primary$gene_polys,
    gene_arrows = primary$gene_labels,
    seq_labels = primary$seq_labels_df,
    show_axis = primary$show_axis
  )

  if (!is.null(plot$ggchord$ref)) plot$ggchord$ref$layout <- primary
  plot$ggchord$layout <- primary
  primary
}


# ====================================================================
# Shared helpers used by ggplot_build.ggchord() and layout preparation.
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
    mapping = mapping %||% ggchord_effective_mapping(lyr), position = lyr$position,
    params = params,
    inherit.aes = lyr$inherit.aes,
    show.legend = lyr$show.legend,
    check.aes = FALSE
  )
  # Preserve the ggchord custom fields on the reconstructed layer
  for (fld in c(
    "ggchord_type", "ggchord_params", "ggchord_placeholder",
    "ggchord_layer_id", "ggchord_input_data", "ggchord_input_mapping",
    "ggchord_role_aes", "ggchord_resolved_input", "ggchord_theme_element",
    "ggchord_theme_components", "ggchord_input_transform", "ggchord_obstacle_provider"
  )) {
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
              axis_text = integer(0), axis = integer(0),
              gene_label_repel = integer(0), seq_label = integer(0),
              seq_region = integer(0))
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    type <- lyr$ggchord_type %||% ""
    if (type %in% names(idx)) idx[[type]] <- c(idx[[type]], i)
  }
  idx
}

#' Build the list of scales for a computed layout
#'
#' @keywords internal
make_ggchord_scales <- function(layout, has_seq = FALSE, has_gene = FALSE,
                                has_feature = FALSE,
                                has_feature_shape = FALSE,
                                plot = NULL,
                                legend_scale = ggchord_device_scale(),
                                legend_text_size = 8,
                                legend_title_size = 9) {
  scales <- list()
  responsive_legend_theme <- ggplot2::theme(
    legend.text = ggplot2::element_text(
      size = legend_text_size * legend_scale
    ),
    legend.title = ggplot2::element_text(
      size = legend_title_size * legend_scale
    )
  )
  role_guide <- function(role, colourbar = FALSE, order = 0,
                         override.aes = list()) {
    if (!is.null(plot)) {
      return(ggchord_role_guide(
        plot, role, colourbar = colourbar, order = order,
        override.aes = override.aes
      ))
    }
    if (colourbar) {
      guide_ggchord_colourbar(order = order, theme = responsive_legend_theme,
                              size_scale = legend_scale)
    } else {
      guide_ggchord_legend(order = order, theme = responsive_legend_theme,
                           size_scale = legend_scale,
                           override.aes = override.aes)
    }
  }

  if (has_seq) {
    scales[[length(scales) + 1]] <- scale_seq_colour_manual(
      name   = "Seq ID",
      values = layout$seq_colors,
      labels = layout$seq_labels,
      breaks = layout$seqs,
      guide  = role_guide("seq", order = 1)
    )
  }

  ribbon_fill_scale <- NULL
  if (!is.null(layout$ribbon_polys)) {
    if (layout$ribbon_color_scheme %in% c("pident", "value")) {
      value_scheme <- identical(layout$ribbon_color_scheme, "value")
      ribbon_name <- if (value_scheme) {
        layout$ribbon_color_name %||% "value"
      } else {
        "Identity (%)"
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
        c(0, 50, 80, 90, 100)
      }
      ribbon_fill_scale <- scale_ribbon_fill_stepsn(
        name    = ribbon_name,
        colours = layout$ribbon_colors,
        limits  = ribbon_limits,
        breaks  = ribbon_breaks,
        guide   = role_guide("ribbon", colourbar = TRUE, order = 2)
      )
    } else {
      ribbon_fill_scale <- scale_ribbon_fill_identity()
    }
  }

  gene_fill_scale <- NULL
  if (has_gene) {
    if (layout$gene_color_scheme == "strand") {
      strand_breaks <- intersect(
        c("+", "-"), unique(as.character(layout$gene_polys$strand))
      )
      if (length(strand_breaks) == 0L) strand_breaks <- c("+", "-")
      gene_fill_scale <- scale_gene_fill_manual(
        name   = "Strand",
        breaks = strand_breaks,
        values = layout$gene_pal,
        guide  = role_guide("gene", order = 3,
          override.aes = list(strand = strand_breaks))
      )
    } else {
      gene_fill_scale <- scale_gene_fill_manual(
        name   = "Gene Annotation",
        breaks = layout$final_gene_order,
        values = layout$gene_pal,
        guide  = role_guide("gene", order = 3)
      )
    }
  }

  # The ribbon layer always uses the internal "zfill" aesthetic (never the
  # plain "fill" aesthetic), so that its mapping stays consistent with the
  # ribbon geom's default aesthetics (which are renamed to avoid injecting a
  # plain "fill" default into the ribbon data). This also keeps ggplot2's
  # guide matching working when no gene layer is present, so the Identity(%)
  # colourbar legend is shown even without gene data.
  feature_fill_scale <- NULL
  merge_feature_guides <- FALSE
  if (has_feature) {
    merge_feature_guides <- has_feature_shape &&
      !is.null(layout$feature_shape_pal) &&
      identical(
        as.character(layout$final_gene_order),
        as.character(layout$feature_shape_order)
      )
    feature_override <- if (merge_feature_guides) {
      list(feature_shape = unname(
        layout$feature_shape_pal[layout$final_gene_order]
      ))
    } else {
      list()
    }
    feature_fill_scale <- scale_feature_fill_manual(
      name = "Feature", breaks = layout$final_gene_order,
      values = layout$gene_pal,
      guide = role_guide("feature", order = 3,
        override.aes = feature_override)
    )
  }
  feature_shape_scale <- NULL
  if (has_feature_shape && !is.null(layout$feature_shape_pal)) {
    feature_shape_scale <- scale_feature_shape_manual(
      name = "Feature",
      breaks = layout$feature_shape_order,
      values = layout$feature_shape_pal,
      guide = if (isTRUE(merge_feature_guides)) "none" else {
        role_guide("feature", order = 3)
      }
    )
  }

  ribbon_aes <- "ribbon_fill"
  if (!is.null(ribbon_fill_scale)) {
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
  if (!is.null(feature_fill_scale)) {
    scales[[length(scales) + 1]] <- feature_fill_scale
  }
  if (!is.null(feature_shape_scale)) {
    scales[[length(scales) + 1]] <- feature_shape_scale
  }

  # Ribbon alpha is a preset value; use an identity scale so it renders as specified
  if (!is.null(layout$ribbon_polys)) {
    scales[[length(scales) + 1]] <- ggplot2::scale_alpha_identity(
      aesthetics = "ribbon_alpha"
    )
  }
  # Optional per-ribbon outline / linetype mappings use internal aesthetics so
  # they do not collide with the sequence colour scale or the ribbon fill scale.
  if (isTRUE(layout$ribbon_use_outline)) {
    scales[[length(scales) + 1]] <- ggplot2::scale_colour_identity(aesthetics = "ribbon_colour")
  }
  if (isTRUE(layout$ribbon_use_linetype)) {
    scales[[length(scales) + 1]] <- ggplot2::scale_linetype_identity(aesthetics = "ribbon_linetype")
  }
  # Sequence-region bands use their own internal fill aesthetic.
  if (!is.null(layout$region_polys) && nrow(layout$region_polys) > 0) {
    scales[[length(scales) + 1]] <- ggplot2::scale_fill_identity(aesthetics = "region_fill")
  }
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

#' Inspect explicit user mappings for default-scale inference
#' @noRd
ggchord_mapping_info <- function(plot, aesthetic) {
  infos <- list()
  for (lyr in plot$layers) {
    mapping <- lyr$ggchord_input_mapping
    if (is.null(mapping) || !aesthetic %in% names(mapping)) next
    expr <- rlang::quo_get_expr(mapping[[aesthetic]])
    label <- rlang::as_label(expr)
    staged <- grepl("after_stat\\s*\\(", label)
    value <- if (staged) NULL else tryCatch(
      rlang::eval_tidy(mapping[[aesthetic]], data = lyr$data),
      error = function(e) NULL
    )
    kind <- if (staged) {
      if (aesthetic %in% c("ribbon_linetype", "feature_shape")) {
        "discrete"
      } else {
        "continuous"
      }
    } else if (is.numeric(value) && !is.factor(value)) {
      "continuous"
    } else {
      "discrete"
    }
    infos[[length(infos) + 1L]] <- list(
      kind = kind, value = value, label = label
    )
  }
  if (!length(infos)) return(NULL)
  kinds <- unique(vapply(infos, `[[`, character(1), "kind"))
  if (length(kinds) > 1L) {
    ggchord_stop(
      "Incompatible continuous and discrete mappings were supplied for `",
      aesthetic, "`"
    )
  }
  values <- unlist(lapply(infos, `[[`, "value"), use.names = FALSE)
  list(kind = kinds, values = values, label = infos[[1L]]$label)
}

#' Replace managed defaults when an explicit mapping needs another scale type
#' @noRd
ggchord_infer_visual_scales <- function(plot, layout, scales) {
  aesthetics <- c(
    "seq_colour", "ribbon_fill", "ribbon_alpha", "ribbon_colour",
    "ribbon_linetype", "gene_fill", "feature_fill", "region_fill"
  )
  for (aesthetic in aesthetics) {
    if (plot$scales$has_scale(aesthetic)) next
    info <- ggchord_mapping_info(plot, aesthetic)
    if (is.null(info)) next
    scales <- Filter(function(s) !aesthetic %in% s$aesthetics, scales)
    if (info$kind == "discrete") {
      levels <- unique(as.character(info$values))
      levels <- levels[!is.na(levels)]
      values <- chord_default_palette(max(length(levels), 1L))
      names(values) <- levels
      scale <- switch(aesthetic,
        seq_colour = scale_seq_colour_manual(
          name = info$label, values = values
        ),
        ribbon_fill = scale_ribbon_fill_manual(
          name = info$label, values = values
        ),
        ribbon_alpha = scale_ribbon_alpha_manual(
          name = info$label,
          values = stats::setNames(seq(0.35, 0.9, length.out = max(length(levels), 1L)), levels)
        ),
        ribbon_colour = scale_ribbon_colour_manual(
          name = info$label, values = values
        ),
        ribbon_linetype = scale_ribbon_linetype_manual(
          name = info$label,
          values = stats::setNames(rep(c(1, 2, 3, 4, 5, 6), length.out = max(length(levels), 1L)), levels)
        ),
        gene_fill = scale_gene_fill_manual(
          name = info$label, values = values
        ),
        feature_fill = scale_feature_fill_manual(
          name = info$label, values = values
        ),
        region_fill = scale_region_fill_manual(
          name = info$label, values = values
        )
      )
    } else {
      colours <- layout$ribbon_colors %||%
        c("#34457E", "#2FA96B", "#8BD925", "#F0E51B")
      scale <- switch(aesthetic,
        ribbon_fill = scale_ribbon_fill_gradientn(
          name = info$label, colours = colours
        ),
        ribbon_alpha = scale_ribbon_alpha_continuous(name = info$label),
        ribbon_colour = ggplot2::scale_colour_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        seq_colour = ggplot2::scale_colour_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        gene_fill = ggplot2::scale_fill_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        feature_fill = ggplot2::scale_fill_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        region_fill = ggplot2::scale_fill_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        ggchord_stop("A continuous mapping is not supported for `", aesthetic, "`")
      )
    }
    role <- switch(aesthetic,
      seq_colour = "seq", ribbon_fill = "ribbon",
      ribbon_alpha = "ribbon", ribbon_colour = "ribbon",
      ribbon_linetype = "ribbon", gene_fill = "gene",
      feature_fill = "feature", region_fill = "region")
    if (!is.null(role)) {
      scale$guide <- ggchord_role_guide(
        plot, role, colourbar = identical(aesthetic, "ribbon_fill") &&
          identical(info$kind, "continuous")
      )
    }
    scales[[length(scales) + 1L]] <- scale
  }
  scales
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
      plot$layers[[idx]] <- reconstruct_layer(lyr, lyr$data, mapping = mp)
    }
  }
  plot
}

#' Set the fixed coordinate system from the layout extremes
#' @keywords internal
set_ggchord_coord <- function(plot, layout) {
  coord <- plot$coordinates
  if (!isTRUE(coord$ggchord_coord)) return(plot)

  lim <- switch(
    coord$fit %||% "labels",
    labels = ggchord_adaptive_limits(layout),
    geometry = ggchord_geometry_limits(layout),
    manual = list(xlim = coord$user_xlim, ylim = coord$user_ylim)
  )
  xlim <- coord$user_xlim %||% lim$xlim
  ylim <- coord$user_ylim %||% lim$ylim

  resolved <- new_coord_chord(
    rotation = coord$rotation %||% 45,
    ratio = coord$ratio %||% 1,
    xlim = xlim,
    ylim = ylim,
    expand = coord$expand %||% TRUE,
    clip = coord$clip %||% "off",
    fit = coord$fit %||% "labels",
    user_xlim = coord$user_xlim,
    user_ylim = coord$user_ylim
  )
  plot$coordinates <- resolved
  plot
}

#' Compute tight coordinate limits for geometry only
#' @noRd
ggchord_geometry_limits <- function(layout) {
  ext <- layout$extremes
  if (is.null(ext) || !all(is.finite(c(ext$x_min, ext$x_max,
                                        ext$y_min, ext$y_max)))) {
    return(list(xlim = c(-1, 1), ylim = c(-1, 1)))
  }
  x_pad <- 0.02 * max(ext$x_max - ext$x_min, 1)
  y_pad <- 0.02 * max(ext$y_max - ext$y_min, 1)
  list(
    xlim = c(ext$x_min - x_pad, ext$x_max + x_pad),
    ylim = c(ext$y_min - y_pad, ext$y_max + y_pad)
  )
}

#' Compute coordinate limits that fit the rendered text boxes
#'
#' Instead of adding one global text-width pad on every side, this helper
#' measures the actual gene, sequence and axis label boxes and expands only
#' the sides that need it. x and y are fitted independently: `coord_fixed()`
#' preserves equal physical units without requiring a square data range. This
#' lets wide or tall rendered content use the available panel more efficiently.
#' @keywords internal
ggchord_adaptive_limits <- function(layout) {
  ext <- layout$extremes
  if (is.null(ext) || !all(is.finite(c(ext$x_min, ext$x_max,
                                        ext$y_min, ext$y_max)))) {
    return(list(xlim = c(-1, 1), ylim = c(-1, 1)))
  }

  units_per_inch <- layout$text_units_per_inch
  if (is.null(units_per_inch) || length(units_per_inch) != 1L ||
      !is.finite(units_per_inch) || units_per_inch <= 0) {
    units_per_inch <- ggchord_device_units_per_inch(
      c(ext$x_min, ext$x_max), c(ext$y_min, ext$y_max)
    )
  }

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

  x_pad <- 0.02 * max(diff(x_lim), 1)
  y_pad <- 0.01 * max(diff(y_lim), 1)

  list(
    xlim = c(x_lim[1] - x_pad, x_lim[2] + x_pad),
    ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad)
  )
}

#' Fully prepare a ggchord plot and return it (compute layout, rename ribbon
#' mappings, attach scales, set coordinates). The layout is cached on the plot
#' (and on the shared reference environment) during preparation. Used by the
#' callers that need a fully prepared ggplot object.
#' @keywords internal
prepare_ggchord_plot <- function(plot) {
  plot$scales$scales <- Filter(function(s) is.null(s$ggchord_managed),
                               plot$scales$scales)
  layout <- compute_chord_geometry(plot)
  cls <- classify_ggchord_layers(plot)
  new_layers <- plot$layers
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (is.null(lyr$ggchord_type)) next
    lyr$ggchord_resolved_input <- layout$layer_inputs[[lyr$ggchord_layer_id]][[lyr$ggchord_type]]
    new_layers[[i]] <- reconstruct_layer(
      lyr, extract_ggchord_layer_data(lyr, layout)
    )
  }
  plot$layers <- new_layers
  has_feature <- any(vapply(plot$layers, function(x) {
    "feature_fill" %in% names(x$mapping)
  }, logical(1)))
  has_feature_shape <- any(vapply(plot$layers, function(x) {
    "feature_shape" %in% names(ggchord_effective_mapping(x))
  }, logical(1)))
  has_gene <- any(vapply(plot$layers, function(x) {
    "gene_fill" %in% names(x$mapping)
  }, logical(1)))
  sc <- make_ggchord_scales(layout,
                            has_seq = length(cls$seq) > 0,
                            has_gene = has_gene,
                            has_feature = has_feature,
                            has_feature_shape = has_feature_shape,
                            plot = plot,
                            legend_text_size = ggchord_theme_point_size(
                              plot, "legend.text", 8
                            ),
                            legend_title_size = ggchord_theme_point_size(
                              plot, "legend.title", 9
                            ))
  sc$scales <- ggchord_infer_visual_scales(plot, layout, sc$scales)
  plot <- rename_ribbon_layers(plot, cls$ribbon, sc$ribbon_aes, layout)
  plot <- attach_ggchord_scales(plot, sc$scales)
  plot <- ggchord_add_link_scales(plot)
  if (!isTRUE(ggchord_plot_settings(plot)$axis$hidden) &&
      nrow(layout$axis_lines %||% data.frame()) > 0L) {
    plot$layers[[length(plot$layers) + 1L]] <- ggchord_axis_layer(layout)
  }
  plot <- set_ggchord_coord(plot, layout)
  plot <- ggchord_apply_theme_styles(plot)
  plot
}

#' @export
ggplot_build.ggchord <- function(plot, ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (is.null(plot$ggchord)) {
    ggchord_stop("Not a valid ggchord object: no data stored on the plot. ",
         "Please build the plot with ggchord().")
  }
  plot <- prepare_ggchord_plot(plot)
  class(plot) <- setdiff(class(plot), "ggchord")
  ggchord_branch_built(ggplot2::ggplot_build(plot))
}
