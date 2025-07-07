#' @title Core Plotting Function Collection
#' @description Contains core functions for integrating graphic elements and generating the final ggplot2 object.
#' @import ggplot2
#' @import ggnewscale
#' @import grid
#' @name chordPlotFunc
#' @title Core Function for Chord Plotting
#' @description Integrates all graphic elements to generate the final chord plot
#' @param allRibbon Ribbon data (data.frame)
#' @param ribbon_alpha Ribbon transparency
#' @param ribbon_color_scheme Ribbon color scheme
#' @param ribbon_colors Ribbon color parameters
#' @param show_legend Whether to display the legend
#' @param gene_polys Gene arrow polygon data
#' @param gene_pal Gene color vector
#' @param gene_color_scheme Gene color scheme
#' @param final_gene_order Order of genes in the legend
#' @param seqArcs Sequence arc data
#' @param gene_arrows Gene label data
#' @param gene_label_show Whether to display gene labels
#' @param gene_label_size Font size of gene labels
#' @param show_axis Whether to display axes
#' @param axisLines Axis line data
#' @param axisTicks Tick mark data
#' @param axisLabelOrientation Orientation of axis labels
#' @param seq_colors Sequence color vector
#' @param seq_labels Sequence labels
#' @param seqs List of sequence IDs
#' @param extremes Plot extrema (min/max coordinates)
#' @param panel_margin Panel margins
#' @param title Plot title
#' @param process_panel_margin Function to process panel margins (internal)
#' @return ggplot2 object
#' @keywords internal
chordPlotFunc <- function(allRibbon, ribbon_alpha, ribbon_color_scheme, ribbon_colors, show_legend, gene_polys, gene_pal, gene_color_scheme, final_gene_order, seqArcs, gene_arrows, gene_label_show, gene_label_size, show_axis, axisLines, axisTicks, axisLabelOrientation, seq_colors, seq_labels, seqs, extremes, panel_margin, title) {
  ggplot() +
    # 1. Draw ribbons and set the first fill scale (only if valid ribbon data exists)
    { if (!is.null(allRibbon))
      geom_polygon(
        data = allRibbon,
        aes(
          x, y, group = group,
          fill = if (ribbon_color_scheme == "pident") pident else fill
        ),
        alpha = ribbon_alpha,
        color = "grey40",
        linewidth = 0.2
      )
    } +
    # Ribbon color scale (only if valid ribbon data exists)
    { if (!is.null(allRibbon)) {
      if (ribbon_color_scheme == "pident") {
        scale_fill_stepsn(
          name = "Identity(%)",
          colours = ribbon_colors,
          limits = c(0, 100),
          breaks = c(0, 50, 80, 90, 95, 100),
          guide = if (show_legend) guide_colorbar(theme = theme(legend.title.position = "top", legend.key.height = unit(120, "mm")), order = 1, position = "left") else "none"
        )
      } else {
        scale_fill_identity(guide = "none")
      }
    }
    } +
    # 2. Reset fill scale (core: use ggnewscale, only if gene data exists)
    { if (nrow(gene_polys) > 0) new_scale_fill() } +

    # 3. Draw gene arrows and set the second fill scale (map strand or anno based on mode)
    { if (nrow(gene_polys) > 0)
      geom_polygon(
        data = gene_polys,
        aes(x = x, y = y, group = group,
            fill = if (gene_color_scheme == "strand") strand else anno), # Dynamically map fill variable
        color = "black", key_glyph = draw_key_gene_arrow
      )
    } +
    # Gene arrow color scale (only if gene data exists)
    { if (nrow(gene_polys) > 0)
      scale_fill_manual(
        name = if (gene_color_scheme == "strand") "Strand" else "Gene Annotation",
        breaks = if (gene_color_scheme == "strand") c("+", "-") else final_gene_order,
        values = gene_pal,
        guide = if (show_legend) guide_legend(order = 3) else "none"
      )
    } +
    # 4. Draw other elements
    # Sequence arcs
    geom_path(
      data = do.call(rbind, seqArcs),
      aes(x, y, group = seq_id, color = seq_id),
      linewidth = 1.2,
      arrow = arrow(angle = 25, length = unit(2, "mm"), type = "closed")
    ) +
    # Gene labels (only if display is enabled and gene data exists)
    { if (nrow(gene_arrows) > 0 && gene_label_show)
      geom_text(
        data = gene_arrows,
        aes(x = text_x, y = text_y, label = text, angle = text_angle,
            hjust = hjust, vjust = vjust), # Use precomputed alignment parameters
        size = gene_label_size,
        color = "black",
        inherit.aes = FALSE
      )
    } +
    # Axis lines (only drawn if show_axis is TRUE)
    { if (show_axis)
      geom_path(
        data = axisLines,
        aes(x, y, group = seq_id),
        color = "black",
        linewidth = 0.3,
        inherit.aes = FALSE
      )
    } +
    # Tick marks (only drawn if show_axis is TRUE)
    { if (show_axis)
      geom_segment(
        data = axisTicks,
        aes(x = x0, y = y0, xend = x1, yend = y1),
        color = "black",
        linewidth = 0.3,
        inherit.aes = FALSE
      )
    } +
    # Tick labels (only drawn if show_axis is TRUE)
    { if (show_axis) {
      # First create a subset of tick data with labels
      label_data <- subset(axisTicks, !is.na(label))
      # Calculate angle for each label, ensuring length matches
      label_angles <- ifelse(
        axisLabelOrientation[label_data$seq_id] == "horizontal",
        0,
        atan2(label_data$label_y - label_data$y0, label_data$label_x - label_data$x0) * 180 / pi + 90 +
          as.numeric(axisLabelOrientation[label_data$seq_id])
      )

      geom_text(
        data = label_data,
        aes(x = label_x, y = label_y, label = label, size = size, group = seq_id),
        inherit.aes = FALSE,
        color = "black",
        angle = label_angles # Use angle vector matching data length
      )
    }
    } +
    scale_size_identity() +
    # Sequence colors
    scale_color_manual(
      name = "Seq ID",
      values = seq_colors,
      labels = seq_labels,
      breaks = seqs,
      guide = if (show_legend) guide_legend(order = 2) else "none"
    ) +

    # Theme settings
    coord_equal(clip = "off", xlim = c(extremes$x_min - process_panel_margin(panel_margin)$l, extremes$x_max + process_panel_margin(panel_margin)$r), ylim = c(extremes$y_min - process_panel_margin(panel_margin)$b, extremes$y_max + process_panel_margin(panel_margin)$t)) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
      legend.background = element_blank(),
      legend.box.spacing = unit(10, "mm"),
      legend.spacing = unit(5, "mm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.background = element_blank()
    )
}
