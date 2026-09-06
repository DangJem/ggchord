#' Build the list of scales for a computed layout
#'
#' @noRd
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
#' @noRd
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
#' @noRd
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

