#' Set the fixed coordinate system from the layout extremes
#' @noRd
set_ggchord_coord <- function(plot, layout, axis_labels = NULL) {
  coord <- plot$coordinates
  if (!isTRUE(coord$ggchord_coord) || !length(layout$seqs)) return(plot)

  lim <- switch(
    coord$fit %||% "labels",
    labels = ggchord_adaptive_limits(layout, axis_labels),
    geometry = ggchord_geometry_limits(layout),
    manual = list(xlim = coord$user_xlim, ylim = coord$user_ylim)
  )
  # Allow room for the final feature outlines, text boxes and leaders. Add a
  # modest margin along the paths so terminal annotations have breathing room.
  if (!is.null(layout$reference_spacing) && !identical(coord$fit, "manual")) {
    panels <- layout$panel_references %||% list(layout)
    along <- if (identical(layout$seq_axis, "vertical")) 2L else 1L
    for (i in 1:2) {
      key <- c("xlim", "ylim")[i]
      boxes <- unlist(lapply(panels, function(p) p$reference_spacing[[i]]$limits))
      lim[[key]] <- range(lim[[key]], boxes)
      if (i == along && isTRUE(plot$ggchord$sequence_layout$auto_row_gap))
        lim[[key]] <- lim[[key]] + c(-.04, .04) * diff(lim[[key]])
    }
  }
  # On a linear layout the transverse spacing is set in millimetres.
  # Once everything fits, keep the spare panel space: fitting this axis
  # tightly again only rescales its arbitrary native units on every rebuild.
  if (inherits(coord, "CoordCollinear") && is.null(coord$ratio) &&
      coord$collinear_seq_axis %in% c("horizontal", "vertical", "sankey") &&
      !identical(coord$fit, "manual")) {
    vertical <- identical(coord$collinear_seq_axis, "vertical") || identical(coord$sankey_axis, "y")
    axis <- if (vertical) "x" else "y"
    key <- paste0(axis, "lim")
    context <- ggchord_render_context()
    span <- unname(context$units_per_mm[axis] * context$panel_mm[if (vertical) "width" else "height"])
    if (length(span) == 1L && is.finite(span) && span > diff(lim[[key]]))
      lim[[key]] <- mean(lim[[key]]) + c(-.5, .5) * span
  }
  if (isTRUE(ggchord_render_context()$contain_panel) && isTRUE(coord$ggchord_circular)) {
    context <- ggchord_render_context()
    span <- unname(context$units_per_mm * context$panel_mm)
    if (span[1L] > diff(lim$xlim)) lim$xlim <- mean(lim$xlim) + c(-.5, .5) * span[1L]
    if (span[2L] > diff(lim$ylim)) lim$ylim <- mean(lim$ylim) + c(-.5, .5) * span[2L]
  }
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
      rotation = coord$rotation %||% 90,
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
  if (isTRUE(coord$ggchord_collinear)) {
    resolved <- new_ggchord_coord(
      CoordCollinear, rotation = 0, ratio = coord$ratio,
      xlim = xlim, ylim = ylim, expand = coord$expand %||% FALSE,
      clip = coord$clip %||% "off", fit = coord$fit %||% "labels",
      user_xlim = coord$user_xlim, user_ylim = coord$user_ylim
    )
    resolved$ggchord_collinear <- TRUE
    for (field in c("collinear_align", "collinear_anchor", "collinear_orientation",
        "collinear_min_support", "collinear_segments", "collinear_segment_gap_mm", "collinear_seq_axis",
        "ggchord_sankey", "sankey_column_gap_mm", "sankey_self_links", "sankey_axis",
        "ggchord_radial", "radial_start_angle", "radial_clockwise", "radial_inner_radius")) resolved[[field]] <- coord[[field]]
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
ggchord_adaptive_limits <- function(layout, axis_labels = NULL) {
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
    if (is.null(axis_labels)) {
      ticks <- layout$axis_ticks
      axis_labels <- ggchord_axis_label_geometry(layout, ticks[!is.na(ticks$label), , drop = FALSE])
    }
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

  for (owner in layout$layer_geometry) {
    labels <- owner$domain_label
    if (!NROW(labels)) next
    labels <- labels[labels$.component == "text", , drop = FALSE]
    if (all(c("bbox_xmin", "bbox_xmax", "bbox_ymin", "bbox_ymax") %in% names(labels))) {
      add_boxes(data.frame(xmin = labels$bbox_xmin, xmax = labels$bbox_xmax,
        ymin = labels$bbox_ymin, ymax = labels$bbox_ymax))
    }
  }

  resources <- layout$circular_annotation_registry
  for (entry in list(list(data = layout$track_axes, size = 2.2),
      list(data = layout$sequence_scales, size = 3))) {
    text <- entry$data
    if (!is.data.frame(text) || !nrow(text) || !"label" %in% names(text)) next
    text <- text[!is.na(text$label) & nzchar(as.character(text$label)), , drop = FALSE]
    if (!nrow(text)) next
    if (".axis_hidden" %in% names(text)) text <- text[!text$.axis_hidden, , drop = FALSE]
    for (name in c("size", "family", "fontface", "angle", "lineheight", "hjust", "vjust"))
      if (paste0(".axis_", name) %in% names(text)) text[[name]] <- text[[paste0(".axis_", name)]]
    text$text_angle <- text$angle %||% 0
    text$size <- text$size %||% entry$size
    text$hjust <- text$hjust %||% .5; text$vjust <- text$vjust %||% .5
    add_boxes(ggchord_text_boxes(text, x_col = "x", y_col = "y", text_col = "label",
      units_per_inch = units_per_inch, box_padding = .03))
  }
  if (is.data.frame(resources) && nrow(resources) &&
      all(c("side", "bbox_xmin", "bbox_xmax", "bbox_ymin",
        "bbox_ymax") %in% names(resources))) {
    outer <- resources[resources$side %in% "outer" &
      is.finite(resources$bbox_xmin) &
      is.finite(resources$bbox_xmax) &
      is.finite(resources$bbox_ymin) &
      is.finite(resources$bbox_ymax), , drop = FALSE]
    if (nrow(outer)) add_boxes(data.frame(
      xmin = outer$bbox_xmin, xmax = outer$bbox_xmax,
      ymin = outer$bbox_ymin, ymax = outer$bbox_ymax
    ))
  }

  # Rounded feature callouts and bold enzyme fragments extend slightly beyond
  # the bare text metrics. Keep a modest physical-looking safety margin so a
  # dense annotation perimeter shrinks the map instead of clipping its last row.
  x_pad <- 0.05 * max(diff(x_lim), 1)
  y_pad <- 0.035 * max(diff(y_lim), 1)

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
  if (isTRUE(plot$coordinates$ggchord_dotplot)) return(ggchord_prepare_dotplot(plot))
  plot <- ggchord_prepare_geometry_inputs(plot)
  layout <- compute_chord_geometry_prepared(plot)
  new_layers <- plot$layers
  track_axes <- list()
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (is.null(lyr$ggchord_type)) next
    lyr$ggchord_resolved_input <- ggchord_restore_instance_data(
      layout$layer_inputs[[lyr$ggchord_layer_id]][[lyr$ggchord_type]], layout$display_instances)
    geometry <- extract_ggchord_layer_data(lyr, layout)
    if (isTRUE(layout$circular) && identical(lyr$ggchord_type, "domain_label")) {
      lyr$ggchord_params$label_segment_layer <- "below"
      lyr$geom_params$bold_callout <- ggchord_domain_bold_callout(lyr, plot)
    }
    if (identical(lyr$ggchord_type, "track")) {
      axis <- geometry$.component %in% c("track_axis_line", "track_axis_tick", "track_axis_text", "track_axis_title")
      if (any(axis)) {
        track_axes[[lyr$ggchord_layer_id]] <- geometry[axis, , drop = FALSE]
        geometry <- geometry[!axis, , drop = FALSE]
      }
    }
    if (isTRUE(lyr$ggchord_params$user_feature_fill)) geometry$feature_fill_explicit <- NULL
    flow <- ggchord_ribbon_flow(lyr, geometry, layout, plot)
    lyr <- flow$layer; geometry <- flow$data
    # An absent component can have no column prototype in the public layout.
    # Draw its empty build-local layer without evaluating nonexistent x/text
    # mappings. Keep the source layer and exported component unchanged.
    empty <- is.data.frame(geometry) && !nrow(geometry)
    new_layers[[i]] <- if (empty) reconstruct_layer(lyr, geometry, mapping = aes()) else
      ggchord_draw_layer(lyr, geometry, lyr$ggchord_resolved_input)
    if (empty) new_layers[[i]]$inherit.aes <- FALSE
    # Domain geometry already owns its final panel, including explicit repeats.
    if (!is.null(layout$facet_keys)) new_layers[[i]]$layout <- NULL

  }
  plot$layers <- new_layers
  if (!is.null(layout$facet_keys)) {
    plot$facet <- ggchord_precomputed_facet(plot$facet, layout$facet_layout)
  }
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
                            has_seq = any(vapply(plot$layers, function(x)
                              "seq_colour" %in% names(x$mapping), logical(1))),
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
  plot <- attach_ggchord_scales(plot, sc$scales)
  plot <- ggchord_add_link_scales(plot)
  axis_labels <- NULL
  if (!isTRUE(ggchord_plot_settings(plot)$axis$hidden) &&
      (nrow(layout$axis_lines %||% data.frame()) > 0L ||
       nrow(layout$axis_ticks %||% data.frame()) > 0L || NROW(layout$axis_heading) > 0L)) {
    axis_data <- ggchord_axis_geometry(layout)
    axis_labels <- axis_data[axis_data$.component %in% c("text", "title"), , drop = FALSE]
    feature_layers <- which(vapply(plot$layers, function(layer)
      identical(layer$ggchord_type, "gene_poly"), logical(1)))
    if (length(feature_layers)) {
      # Automatic marks are background decoration. Keep the user's layers in
      # their declared order, but let features cover ticks that meet them.
      marks <- ggchord_axis_layer(layout, data = axis_data, components = c("line", "major_tick", "minor_tick"))
      if (nrow(marks$data)) plot$layers <- append(
        plot$layers, list(marks), after = feature_layers[1L] - 1L)
      # Text has its own measured clearance and remains in the foreground.
      text <- ggchord_axis_layer(layout, c("text", "title"), data = axis_data)
      if (nrow(text$data)) plot$layers[[length(plot$layers) + 1L]] <- text
    } else {
      plot$layers[[length(plot$layers) + 1L]] <- ggchord_axis_layer(layout, data = axis_data)
    }
  }
  plot <- set_ggchord_coord(plot, layout, axis_labels)
  plot$theme <- plot$theme + theme(axis.title.x = element_blank(),
    axis.title.y = element_blank())
  sequence_theme <- plot$scales$get_scales("seq_position")$ggchord_axis_guide$theme
  if (!is.null(sequence_theme)) {
    local <- setNames(as.list(sequence_theme), paste0("ggchord.", names(sequence_theme)))
    axes_only <- plot
    axes_only$theme <- plot$theme + do.call(theme, local)
    plot$layers <- ggchord_apply_theme_styles(axes_only)$layers
  } else {
    plot <- ggchord_apply_theme_styles(plot)
  }
  # Axes use the owner's resolved theme and draw immediately above its marks.
  # They contain no observations with which to train data scales or guides.
  if (length(track_axes)) {
    for (i in rev(seq_along(plot$layers))) {
      owner <- plot$layers[[i]]
      id <- owner$ggchord_layer_id
      if (is.null(id) || is.null(track_axes[[id]])) next
      axis_data <- track_axes[[id]]
      if (is.null(axis_data$xend)) axis_data$xend <- axis_data$yend <- NA_real_
      axis <- ggplot2::layer(data = axis_data,
        mapping = aes(x = x, y = y, xend = xend, yend = yend,
          group = group, label = label, .component = I(.component)),
        geom = owner$geom, stat = "identity", position = "identity", params = owner$geom_params,
        layout = owner$layout, inherit.aes = FALSE, show.legend = FALSE, check.aes = FALSE)
      plot$layers <- append(plot$layers, list(axis), after = i)
    }
  }
  layout$render_context <- ggchord_render_context() %||% ggchord_estimated_render_context(layout)
  plot$ggchord$ref$layout <- layout
  plot
}

#' @export
print.ggchord <- function(x, ...) {
  ggchord_local_console_errors()
  NextMethod()
}

#' @export
ggplot_build.ggchord <- function(plot, ...) {
  ggchord_local_console_errors()
  if (is.null(plot$ggchord)) {
    ggchord_stop("Not a valid ggchord object: no data stored on the plot. ",
         "Please build the plot with ggchord().")
  }
  source <- ggchord_render_track_source(ggchord_render_inputs(plot))
  plot <- ggchord_apply_adaptive_styles(source)
  plot$ggchord$ref <- new.env(parent = emptyenv())
  plot <- prepare_ggchord_plot(plot)
  plot$ggchord$ref$layout$adaptive_layout <- plot$ggchord$adaptive_report
  scopes <- ggchord_drawing_scopes(plot)
  user_scales <- Filter(function(scale) !isTRUE(scale$ggchord_managed), plot$scales$scales)
  user_scopes <- Filter(function(scope) !isTRUE(scope$automatic), scopes)
  explicit_aes <- c(unlist(lapply(user_scales, `[[`, "aesthetics"), use.names = FALSE),
    unlist(lapply(user_scopes, `[[`, "aliases"), use.names = FALSE))
  plot <- ggchord_apply_scoped_scales(plot, scopes)
  plot <- ggchord_track_mapping_titles(plot)
  plot <- ggchord_add_visual_scales(plot, explicit_aes)
  class(plot) <- setdiff(class(plot), "ggchord")
  built <- ggchord_branch_built(ggplot_build(plot))
  ggchord_cache_layout(source, built$plot$ggchord$ref$layout)
  attr(built, "ggchord_source_plot") <- source
  class(built) <- unique(c("ggchord_built", class(built)))
  built
}
