#' Set the fixed coordinate system from the layout extremes
#' @noRd
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
  if (isTRUE(coord$ggchord_circular)) {
    resolved <- new_ggchord_coord(
      CoordCircular,
      rotation = coord$rotation %||% 0,
      ratio = coord$ratio %||% 1,
      xlim = xlim, ylim = ylim,
      expand = coord$expand %||% FALSE,
      clip = coord$clip %||% "off",
      fit = coord$fit %||% "labels",
      user_xlim = coord$user_xlim,
      user_ylim = coord$user_ylim
    )
    resolved$ggchord_circular <- TRUE
    resolved$circular_gap <- coord$circular_gap
    resolved$circular_direction <- coord$circular_direction
  }
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
#' @noRd
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
  restriction <- layout$restriction_sites
  if (!is.null(restriction) && nrow(restriction) > 0L &&
      all(c(".component", "label") %in% names(restriction))) {
    restriction <- restriction[
      restriction$.component == "label" & !is.na(restriction$label),,
      drop = FALSE
    ]
    if (nrow(restriction)) {
      add_boxes(ggchord_text_boxes(
        restriction,
        x_col = "x", y_col = "y", text_col = "label",
        angle_col = "angle", size_col = "size",
        hjust_col = "hjust", vjust_col = "vjust",
        units_per_inch = units_per_inch, box_padding = 0.03
      ))
    }
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
#' @noRd
prepare_ggchord_plot <- function(plot) {
  plot$scales$scales <- Filter(function(s) is.null(s$ggchord_managed),
                               plot$scales$scales)
  layout <- compute_chord_geometry(plot)
  cls <- classify_ggchord_layers(plot)
  new_layers <- plot$layers
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (is.null(lyr$ggchord_type)) next
    if (isTRUE(plot$coordinates$ggchord_circular) &&
        identical(lyr$ggchord_type, "seq") &&
        !isTRUE(lyr$ggchord_params$seq_arrow_supplied)) {
      lyr$geom_params$arrow <- NULL
    }
    if (isTRUE(plot$coordinates$ggchord_circular) &&
        identical(lyr$ggchord_type, "seq") &&
        !isTRUE(lyr$ggchord_params$seq_legend_supplied)) {
      lyr$show.legend <- FALSE
    }
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
