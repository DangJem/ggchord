# ggplotly.R - plotly conversion
#
# plotly::ggplotly() rebuilds a ggplot with its own pipeline. That pipeline
# cannot cope with ggchord's two independent fill scales (ribbon + gene) or the
# custom internal aesthetics the ribbon layer uses (they raise "Discrete value
# supplied to a continuous scale" / "Continuous value supplied to a discrete
# scale"). To make interactive charts work, ggchord provides its own
# ggplotly.ggchord() method: it converts the computed chord geometry into an
# equivalent standard ggplot2 plot whose colors are pre-mapped (identity
# scales), lets plotly convert that plot, and then appends legend-only traces
# so the interactive chart still shows the Seq ID / Strand / Identity legends.

#' Convert a ggchord plot to a plotly object
#'
#' ggchord provides an S3 method for \code{plotly::ggplotly()} so that chord
#' diagrams built with \code{ggchord()} (including plots that combine the
#' ribbon and the gene layers) can be converted to interactive \pkg{plotly}
#' charts. The conversion first renders the computed geometry with standard
#' ggplot2 layers whose colors are pre-mapped, then delegates to
#' \code{plotly::ggplotly()}.
#'
#' @param p A ggchord plot object created with \code{ggchord()}.
#' @param ... Additional arguments passed to \code{plotly::ggplotly()}.
#'
#' @return A plotly object.
#' @exportS3Method plotly::ggplotly
ggplotly.ggchord <- function(p, ...) {
  out <- ggchord_plotly_ggplot(p)
  layout <- out$layout
  pl <- plotly::ggplotly(out$plot, ...)
  if (isTRUE(p$ggchord$global$show_legend)) {
    pl <- tryCatch(ggchord_plotly_legend(pl, layout, out$colors),
                   error = function(e) pl)
  }
  # The standard plot uses identity scales with show.legend = FALSE, so plotly
  # disables the legend at the layout level. Turn it back on (the legend-only
  # traces added above provide the Seq ID / Strand / Identity entries).
  pl$x$layout$showlegend <- TRUE
  # geom_seq draws arrowheads that plotly cannot reproduce on line traces;
  # add them as plotly annotations so the arc direction stays visible.
  tryCatch(pl <- ggchord_plotly_arrows(pl, layout, out$colors),
           error = function(e) NULL)
  pl
}

#' Build a plotly-ready standard ggplot2 plot from a ggchord plot
#'
#' The returned plot uses standard geoms with pre-mapped colors and identity
#' scales, so plotly::ggplotly() can convert every element (sequence arcs,
#' ribbons, gene arrows, axis) without scale conflicts. User-supplied color
#' scales (added with \code{+}) are honored.
#' @keywords internal
ggchord_plotly_ggplot <- function(p) {
  if (is.null(p$ggchord$ref) || is.null(p$ggchord$ref$layout)) {
    p <- prepare_ggchord_plot(p)
  }
  layout <- p$ggchord$ref$layout %||% p$ggchord$layout
  cols <- ggchord_plotly_colors(p, layout)

  std <- ggplot2::ggplot() +
    p$theme +
    p$coordinates +
    ggplot2::labs(title = p$labels$title %||% NULL)

  # ---- sequence arcs (one layer per sequence; fixed colour per arc so no
  # "colour" data column is left over for plotly to complain about) ----
  if (length(layout$seq_arcs) > 0) {
    for (i in seq_along(layout$seq_arcs)) {
      arc <- layout$seq_arcs[[i]]
      sid <- unique(arc$seq_id)[1]
      std <- std + ggplot2::geom_path(
        data = arc,
        mapping = ggplot2::aes(x = x, y = y, group = seq_id),
        inherit.aes = FALSE,
        show.legend = FALSE,
        colour = unname(cols$seq[sid]),
        linewidth = 1.2
      )
    }
  }

  # ---- alignment ribbons ----
  if (!is.null(layout$ribbon_polys)) {
    rp <- layout$ribbon_polys
    rp$fill_col <- ggchord_ribbon_fill(layout, cols$ribbon)
    outline <- ggchord_layer_params(p, "ribbon")
    std <- std + ggplot2::geom_polygon(
      data = rp,
      mapping = ggplot2::aes(x = x, y = y, group = group,
                             fill = fill_col, alpha = alpha),
      inherit.aes = FALSE,
      show.legend = FALSE,
      colour = outline$ribbon_outline_color %||% "black",
      linewidth = outline$ribbon_outline_width %||% 0.05,
      linetype = outline$ribbon_outline_linetype %||% 1
    )
  }

  # ---- gene arrows ----
  if (nrow(layout$gene_polys) > 0) {
    gp <- layout$gene_polys
    gp$fill_col <- ggchord_gene_fill(layout, cols$gene)
    std <- std + ggplot2::geom_polygon(
      data = gp,
      mapping = ggplot2::aes(x = x, y = y, group = group, fill = fill_col),
      inherit.aes = FALSE,
      show.legend = FALSE,
      colour = "black"
    )
  }

  # ---- gene labels ----
  if (nrow(layout$gene_labels) > 0) {
    # leader lines for repelled labels
    if (nrow(layout$gene_label_segments) > 0) {
      std <- std + ggplot2::geom_segment(
        data = layout$gene_label_segments,
        mapping = ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1,
                               group = group),
        inherit.aes = FALSE,
        colour = "grey50",
        linewidth = 0.3
      )
    }
    std <- std + ggplot2::geom_text(
      data = layout$gene_labels,
      mapping = ggplot2::aes(x = text_x, y = text_y, label = text,
                             angle = text_angle, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  # ---- sequence labels ----
  if (nrow(layout$seq_labels_df) > 0) {
    std <- std + ggplot2::geom_text(
      data = layout$seq_labels_df,
      mapping = ggplot2::aes(x = text_x, y = text_y, label = label,
                             angle = text_angle, size = size,
                             hjust = hjust, vjust = vjust),
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  # ---- axis ----
  if (nrow(layout$axis_lines) > 0) {
    std <- std + ggplot2::geom_path(
      data = layout$axis_lines,
      mapping = ggplot2::aes(x = x, y = y, group = seq_id),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.3
    )
  }
  if (nrow(layout$axis_ticks) > 0) {
    std <- std + ggplot2::geom_segment(
      data = layout$axis_ticks,
      mapping = ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.3
    )
    axis_lbl <- layout$axis_ticks[!is.na(layout$axis_ticks$label), , drop = FALSE]
    if (nrow(axis_lbl) > 0) {
      std <- std + ggplot2::geom_text(
        data = axis_lbl,
        mapping = ggplot2::aes(x = label_x, y = label_y,
                               label = label, size = size),
        inherit.aes = FALSE,
        colour = "black",
        angle = 0
      )
    }
  }

  # Identity scales pass the pre-mapped colors straight through to plotly.
  list(plot = std + ggplot2::scale_fill_identity() +
         ggplot2::scale_colour_identity() +
         ggplot2::scale_alpha_identity() +
         ggplot2::scale_size_identity(),
       colors = cols,
       layout = layout)
}

#' Resolve the colors used by the plotly conversion (defaults from the layout,
#' overridden by any user-supplied scale added with \code{+}).
#' @keywords internal
ggchord_plotly_colors <- function(p, layout) {
  seq_cols <- layout$seq_colors
  gene_cols <- layout$gene_pal
  ribbon_cols <- layout$ribbon_colors

  for (sc in p$scales$scales) {
    if (!is.null(sc$ggchord_managed)) next
    aes <- sc$aesthetics[1]
    if (is.null(aes)) next
    if (aes == "colour") {
      tryCatch({
        sc$train(layout$seqs)
        mapped <- sc$map(layout$seqs)
        if (all(!is.na(mapped))) seq_cols <- stats::setNames(mapped, layout$seqs)
      }, error = function(e) NULL)
    } else if (aes == "fill") {
      keys <- if (identical(layout$gene_color_scheme, "strand")) {
        c("+", "-")
      } else {
        layout$final_gene_order
      }
      tryCatch({
        sc$train(keys)
        mapped <- sc$map(keys)
        if (all(!is.na(mapped))) gene_cols <- stats::setNames(mapped, keys)
      }, error = function(e) NULL)
    }
  }

  list(seq = seq_cols, gene = gene_cols, ribbon = ribbon_cols)
}

#' Pre-map the ribbon polygon colors
#' @keywords internal
ggchord_ribbon_fill <- function(layout, ribbon_colors = NULL) {
  rp <- layout$ribbon_polys
  if (layout$ribbon_color_scheme == "pident") {
    # Interpolate the pident value (0-100) across the ribbon color gradient.
    ramp <- grDevices::colorRamp(ribbon_colors %||% layout$ribbon_colors)
    t <- pmin(1, pmax(0, rp$pident / 100))
    m <- ramp(t)
    grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  } else {
    rp$fill
  }
}

#' Pre-map the gene arrow colors
#' @keywords internal
ggchord_gene_fill <- function(layout, gene_colors = NULL) {
  gp <- layout$gene_polys
  pal <- gene_colors %||% layout$gene_pal
  key <- if (layout$gene_color_scheme == "strand") as.character(gp$strand) else as.character(gp$anno)
  unname(pal[key])
}

#' Extract stored parameters of the first layer of a given ggchord type
#' @keywords internal
ggchord_layer_params <- function(p, type) {
  for (lyr in p$layers) {
    pp <- lyr$ggchord_params
    if (!is.null(pp) && identical(pp$type, type)) return(pp)
  }
  list()
}

#' Add arrowhead annotations at the tip of every sequence arc
#'
#' plotly's scatter traces cannot draw line arrowheads, so the directional
#' arrows used by \code{geom_seq()} are reproduced as plotly annotations.
#' @keywords internal
ggchord_plotly_arrows <- function(pl, layout, colors = NULL) {
  if (length(layout$seq_arcs) == 0) return(pl)
  seq_cols <- colors$seq %||% layout$seq_colors
  ext <- layout$extremes
  span <- max(ext$x_max - ext$x_min, ext$y_max - ext$y_min, 1)
  arrow_len <- 0.05 * span

  ann <- pl$x$layout$annotations
  if (is.null(ann)) ann <- list()

  for (arc in layout$seq_arcs) {
    n <- nrow(arc)
    if (n < 3) next
    sid <- unique(arc$seq_id)[1]
    tip_x <- arc$x[n]
    tip_y <- arc$y[n]
    dx <- tip_x - arc$x[n - 1]
    dy <- tip_y - arc$y[n - 1]
    len <- sqrt(dx^2 + dy^2)
    if (len == 0) next
    # arrow tail sits behind the tip, pointing along the arc's direction
    tail_x <- tip_x - dx / len * arrow_len
    tail_y <- tip_y - dy / len * arrow_len
    ann[[length(ann) + 1]] <- list(
      x = tip_x, y = tip_y, xref = "x", yref = "y",
      ax = tail_x, ay = tail_y, axref = "x", ayref = "y",
      showarrow = TRUE, arrowhead = 2, arrowsize = 1,
      arrowcolor = unname(seq_cols[sid]),
      text = "", hoverinfo = "none"
    )
  }

  pl$x$layout$annotations <- ann
  pl
}

#' Append legend-only traces (Seq ID, Strand/Annotation, Identity) to a plotly
#' object built from a ggchord plot.
#' @keywords internal
ggchord_plotly_legend <- function(pl, layout, colors = NULL) {
  seq_cols <- colors$seq %||% layout$seq_colors
  gene_cols <- colors$gene %||% layout$gene_pal
  ribbon_cols <- colors$ribbon %||% layout$ribbon_colors

  # Seq ID
  if (length(layout$seq_arcs) > 0) {
    for (sid in layout$seqs) {
      pl <- plotly::add_trace(
        pl, type = "scatter", mode = "markers", x = 0, y = 0,
        marker = list(color = unname(seq_cols[sid]), size = 12),
        name = unname(layout$seq_labels[sid]),
        legendgroup = "Seq ID",
        visible = "legendonly", showlegend = TRUE, hoverinfo = "none"
      )
    }
  }

  # Strand / Gene Annotation
  if (nrow(layout$gene_polys) > 0) {
    if (layout$gene_color_scheme == "strand") {
      items <- c("+", "-")
      grp <- "Strand"
    } else {
      items <- layout$final_gene_order
      grp <- "Gene Annotation"
    }
    for (it in items) {
      pl <- plotly::add_trace(
        pl, type = "scatter", mode = "markers", x = 0, y = 0,
        marker = list(color = unname(gene_cols[it]), size = 12),
        name = it,
        legendgroup = grp,
        visible = "legendonly", showlegend = TRUE, hoverinfo = "none"
      )
    }
  }

  # Identity(%) colorbar
  if (!is.null(layout$ribbon_polys) &&
      identical(layout$ribbon_color_scheme, "pident")) {
    cols <- ribbon_cols
    n <- length(cols)
    stops <- lapply(seq_len(n), function(i) {
      list((i - 1) / max(n - 1, 1), cols[i])
    })
    pl <- plotly::add_trace(
      pl, type = "scatter", mode = "markers", x = 0, y = 0,
      marker = list(
        color = seq(0, 100, length.out = 100),
        colorscale = stops,
        colorbar = list(title = "Identity(%)", len = 0.6)
      ),
      visible = "legendonly", showlegend = FALSE,
      name = "Identity(%)", hoverinfo = "none"
    )
  }

  pl
}
