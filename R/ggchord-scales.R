ggchord_uses_plasmid_feature_scale <- function(circular, levels) {
  isTRUE(circular) && all(levels %in% c(
    "CDS", "promoter", "rep_origin", "terminator", "protein_bind", "RBS",
    "repeat_region", "misc_feature"
  ))
}

ggchord_default_feature_shapes <- function(levels) {
  shapes <- c(
    "arrow", "compact_arrow", "promoter_arrow", "primer_arrow", "marker",
    "block", "capped_line", "primer_arc", "chevron", "lollipop"
  )
  values <- if (all(levels %in% shapes)) levels else rep(shapes, length.out = length(levels))
  setNames(values, levels)
}

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
  feature_fill_order <- layout$feature_fill_order %||% layout$final_gene_order
  feature_fill_pal <- layout$feature_fill_pal %||% layout$gene_pal
  responsive_legend_theme <- theme(
    legend.text = element_text(
      size = legend_text_size * legend_scale
    ),
    legend.title = element_text(
      size = legend_title_size * legend_scale
    )
  )
  role_guide <- function(role, colourbar = FALSE, order = 0,
                         override.aes = list()) {
    labels <- switch(role, seq=layout$legend_seq_labels %||% layout$seq_labels,
      gene=names(layout$gene_pal), feature=feature_fill_order, region=names(layout$region_pal), NULL)
    if (!is.null(plot)) {
      return(ggchord_role_guide(
        plot, role, colourbar = colourbar, order = order,
        override.aes = override.aes, labels=labels
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
      values = layout$legend_seq_colors %||% layout$seq_colors,
      labels = layout$legend_seq_labels %||% layout$seq_labels,
      breaks = layout$legend_seqs %||% layout$seqs,
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
        c(0, 50, 80, 90, 95, 100)
      }
      ribbon_fill_scale <- scale_ribbon_fill_stepsn(
        name    = ribbon_name,
        colours = layout$ribbon_colors,
        limits  = ribbon_limits,
        breaks  = ribbon_breaks,
        guide   = if (!value_scheme) ggchord_identity_guide(plot) else role_guide("ribbon", colourbar = TRUE, order = 2)
      )
    } else {
      ribbon_fill_scale <- scale_ribbon_fill_identity()
    }
  }

  gene_fill_scale <- NULL
  if (has_gene) {
    if (layout$gene_color_scheme == "strand") {
      strand_breaks <- c("+", "-")
      gene_fill_scale <- scale_gene_fill_manual(
        name   = "Strand",
        limits = strand_breaks,
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

  feature_fill_scale <- NULL
  merge_feature_guides <- FALSE
  if (has_feature) {
    merge_feature_guides <- has_feature_shape &&
      !is.null(layout$feature_shape_pal) &&
      identical(
        as.character(feature_fill_order),
        as.character(layout$feature_shape_order)
      )
    feature_override <- if (merge_feature_guides) {
      list(feature_shape = unname(
        layout$feature_shape_pal[feature_fill_order]
      ))
    } else {
      list()
    }
    feature_fill_scale <- if (ggchord_uses_plasmid_feature_scale(layout$circular, feature_fill_order)) {
      scale_feature_fill_plasmid(
        name = "Feature", breaks = feature_fill_order,
        guide = role_guide("feature", order = 3,
          override.aes = feature_override)
      )
    } else scale_feature_fill_manual(
        name = "Feature", breaks = feature_fill_order,
        values = feature_fill_pal,
        guide = role_guide("feature", order = 3,
          override.aes = feature_override)
      )
  }
  feature_shape_scale <- NULL
  if (has_feature_shape && !is.null(layout$feature_shape_pal)) {
    feature_shape_scale <- if (ggchord_uses_plasmid_feature_scale(layout$circular, layout$feature_shape_order)) {
      scale_feature_shape_plasmid(
        name = "Feature", breaks = layout$feature_shape_order,
        guide = if (isTRUE(merge_feature_guides)) "none" else
          role_guide("feature", order = 3)
      )
    } else scale_feature_shape_manual(
        name = "Feature", breaks = layout$feature_shape_order,
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
    scales[[length(scales) + 1]] <- scale_fill_identity(aesthetics = "region_fill")
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
    if (aesthetic %in% names(lyr$aes_params)) next
    mapping <- lyr$ggchord_input_mapping
    if (is.null(mapping) || !aesthetic %in% names(mapping)) next
    expr <- get_expr(mapping[[aesthetic]])
    label <- as_label(expr)
    staged <- ggchord_has_stage(expr, "after_stat")
    columns <- lyr$ggchord_params$visual_columns
    value <- if (staged) {
      NULL
    } else if (aesthetic %in% names(columns)) {
      lyr$data[[columns[[aesthetic]]]]
    } else {
      tryCatch(eval_tidy(mapping[[aesthetic]], data = lyr$data),
        error = function(e) NULL)
    }
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
  list(kind = kinds, values = values, label = infos[[1L]]$label,
    labels = unique(vapply(infos, `[[`, character(1), "label")))
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
    if (identical(aesthetic, "ribbon_fill") && info$kind == "continuous" &&
        identical(layout$ribbon_color_scheme, "pident") &&
        all(info$labels == "pident")) next
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
          values = setNames(seq(0.35, 0.9, length.out = max(length(levels), 1L)), levels)
        ),
        ribbon_colour = scale_ribbon_colour_manual(
          name = info$label, values = values
        ),
        ribbon_linetype = scale_ribbon_linetype_manual(
          name = info$label,
          values = setNames(rep(c(1, 2, 3, 4, 5, 6), length.out = max(length(levels), 1L)), levels)
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
        ribbon_colour = scale_colour_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        seq_colour = scale_colour_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        gene_fill = scale_fill_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        feature_fill = scale_fill_gradientn(
          name = info$label, colours = colours, aesthetics = aesthetic
        ),
        region_fill = scale_fill_gradientn(
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
          identical(info$kind, "continuous"),
        labels=if(identical(info$kind,"discrete")) unique(as.character(info$values)) else NULL
      )
    }
    if (identical(info$kind, "continuous") && all(info$labels == "pident")) {
      scale$name <- "Identity (%)"
      scale$limits <- c(0, 100)
      scale$breaks <- c(0, 50, 80, 90, 95, 100)
    }
    scales[[length(scales) + 1L]] <- scale
  }
  scales
}

#' Train a visual scale on selected layers
#'
#' Apply one visual scale to the selected layers, overriding both the
#' plot-wide scale and their independent defaults. See [auto_id()].
#'
#' @details Unscoped backbones and features can use the sequence palette when
#' mapped directly to `accver`. Features mapped to `strand` use red (+),
#' blue (-), purple (+/-) and grey (unknown). Circular features mapped to
#' `type` use the plasmid palette with stable colours for additional types.
#' An explicit scale takes precedence.
#'
#' @param layers One or more explicitly named layer IDs.
#' @param scale An ordinary ggplot2 visual scale for those layers.
#' @return A scale declaration for a ggchord plot or a Fusion root.
#' @export
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' gene_colours <- c(CDS = "#187B77", tRNA = "#CC7D3A")
#'
#' ggchord(data = seqs) + geom_seq() +
#'   geom_gene(data = genes, aes(fill = type), id = "genes") +
#'   geom_link_ribbon(data = links, aes(fill = pident), id = "matches") +
#'   scale_for("genes", scale_fill_manual(values = gene_colours, name = "Type")) +
#'   scale_for("matches", scale_fill_gradient(limits = c(0, 100), name = "Identity (%)")) +
#'   coord_chord()
scale_for <- function(layers, scale) {
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_id_vector(layers, "layers")
  if (!inherits(scale, "Scale") || !length(scale$aesthetics))
    ggchord_stop("scale_for(): scale must be an ordinary ggplot2 visual scale")
  visual <- c("fill", "colour", "alpha", "linewidth", "linetype", "shape",
    "size", "stroke", "fill_from", "fill_to")
  if (any(!scale$aesthetics %in% visual))
    ggchord_stop("scale_for(): position and domain scales cannot be scoped to layers")
  structure(list(layers = layers, scale = scale), class = "ggchord_scoped_scale")
}

ggchord_scale_scopes <- function(plot, scopes = plot$ggchord$scoped_scales) {
  if (!length(scopes)) return(list())
  ids <- vapply(plot$layers, function(layer) layer$ggchord_layer_id %||% "", character(1))
  claimed <- character()
  lapply(seq_along(scopes), function(i) {
    scope <- scopes[[i]]
    selected <- match(scope$layers, ids)
    if (anyNA(selected)) ggchord_stop("scale_for(): unknown layer '",
      scope$layers[which(is.na(selected))[1L]], "'")
    aes <- scope$scale$aesthetics
    for (j in selected) {
      layer <- plot$layers[[j]]
      native <- layer$ggchord_params$native_geom
      geom <- if (is.null(native)) layer$geom else
        ggplot2::layer(geom = native, stat = "identity", position = "identity")$geom
      unsupported <- setdiff(aes, geom$aesthetics())
      if (length(unsupported)) ggchord_stop("scale_for(): layer '", ids[j],
        "' does not support ", unsupported[1L])
    }
    keys <- as.vector(outer(selected, aes, paste, sep = "/"))
    overlap <- match(TRUE, keys %in% claimed)
    if (!is.na(overlap)) {
      layer_index <- selected[(overlap - 1L) %% length(selected) + 1L]
      aesthetic <- aes[(overlap - 1L) %/% length(selected) + 1L]
      ggchord_stop("scale_for(): layer '", ids[layer_index],
        "' already has a scoped scale for `", aesthetic, "`")
    }
    claimed <<- c(claimed, keys)
    list(selected = selected, aesthetics = aes, scale = scope$scale,
      aliases = setNames(paste0("ggchord_scope_", i, "_", aes), aes))
  })
}

# Only competing, unclaimed mappings need a private scale. Keep an ordinary
# scale addressable under its public name whenever a single layer uses it.
ggchord_independent_layers <- function(plot, aesthetic, mapped, explicit) {
  if (plot$scales$has_scale(aesthetic)) return(integer())
  for (scope in explicit) if (aesthetic %in% scope$aesthetics)
    mapped <- setdiff(mapped, scope$selected)
  if (length(mapped) < 2L) return(integer())
  Filter(function(i) !ggchord_shared_id(plot$layers[[i]]$ggchord_params$id), mapped)
}

ggchord_trains_visual_scale <- function(mapping) {
  start <- ggchord_start_mapping(mapping)
  (!is.null(start) && !is_call(get_expr(start), "I")) ||
    ggchord_has_stage(mapping, "after_stat")
}

ggchord_drawing_scopes <- function(plot, scopes = plot$ggchord$scoped_scales) {
  specs <- ggchord_scale_scopes(plot, scopes)
  visual <- c("colour", "fill", "alpha", "size", "linewidth", "linetype", "shape", "stroke")
  mappings <- lapply(plot$layers, function(layer) {
    mapping <- layer$mapping
    if (isTRUE(layer$inherit.aes)) {
      for (key in setdiff(names(plot$mapping), names(mapping))) mapping[[key]] <- plot$mapping[[key]]
    }
    mapping <- mapping[intersect(setdiff(names(mapping), names(layer$aes_params)), visual)]
    Filter(ggchord_trains_visual_scale, mapping)
  })
  for (key in visual) {
    mapped <- which(vapply(mappings, function(m) key %in% names(m), logical(1)))
    selected <- ggchord_independent_layers(plot, key, mapped, specs)
    for (i in selected) {
      role <- as_label(mappings[[i]][[key]])
      scale <- if (key %in% c("colour", "fill")) ggchord_mapped_scale(plot, plot$layers[i], key, role)
      index <- length(specs) + 1L
      specs[[index]] <- list(selected = i, aesthetics = key, scale = scale, automatic = TRUE,
        aliases = setNames(paste0("ggchord_scope_", index, "_", key), key))
    }
  }
  specs
}

ggchord_alias_aes <- function(x, aliases) {
  if (is.null(x)) return(NULL)
  fields <- names(x)
  selected <- fields %in% names(aliases)
  fields[selected] <- unname(aliases[fields[selected]])
  names(x) <- fields
  x
}

ggchord_scope_layer <- function(layer, aliases, plot_mapping) {
  geom <- layer$geom
  inverse <- setNames(names(aliases), aliases)
  changed <- function(x) {
    at <- x %in% names(aliases)
    x[at] <- unname(aliases[x[at]])
    x
  }
  copy <- ggproto(NULL, layer)
  mapping <- layer$mapping %||% aes()
  if (isTRUE(layer$inherit.aes)) {
    for (key in setdiff(names(plot_mapping), names(mapping))) mapping[[key]] <- plot_mapping[[key]]
  }
  copy$inherit.aes <- FALSE
  copy$mapping <- ggchord_alias_aes(mapping, aliases)
  columns <- layer$ggchord_params$visual_columns
  if (length(columns)) copy$ggchord_params$visual_columns <- ggchord_alias_aes(columns, aliases)
  copy$aes_params <- ggchord_alias_aes(layer$aes_params, aliases)
  copy$show.legend <- ggchord_alias_aes(layer$show.legend, aliases)
  copy$geom <- ggproto(NULL, geom,
    required_aes = changed(geom$required_aes),
    non_missing_aes = changed(geom$non_missing_aes),
    optional_aes = changed(geom$optional_aes),
    default_aes = ggchord_alias_aes(geom$default_aes, aliases),
    setup_data = function(data, params) {
      data <- geom$setup_data(ggchord_alias_aes(data, inverse), params)
      ggchord_alias_aes(data, aliases)
    },
    use_defaults = function(data, params = list(), modifiers = aes(), default_aes = NULL,
        theme = NULL, ...) {
      data <- geom$use_defaults(ggchord_alias_aes(data, inverse),
        ggchord_alias_aes(params, inverse), ggchord_alias_aes(modifiers, inverse),
        ggchord_alias_aes(default_aes, inverse), theme = theme, ...)
      ggchord_alias_aes(data, aliases)
    },
    draw_layer = function(data, params, layout, coord) {
      geom$draw_layer(ggchord_alias_aes(data, inverse), params, layout, coord)
    },
    draw_key = function(data, params, size) {
      geom$draw_key(ggchord_alias_aes(data, inverse), params, size)
    })
  copy
}

ggchord_scope_guide <- function(guide, aliases) {
  if (is.character(guide) && length(guide) == 1L &&
      guide %in% c("colourbar", "colorbar", "coloursteps", "colorsteps"))
    guide <- getExportedValue("ggplot2", paste0("guide_", guide))()
  if (inherits(guide, "Guide") && !is.null(guide$available_aes) && !"any" %in% guide$available_aes) {
    available <- aliases[intersect(names(aliases), guide$available_aes)]
    original <- guide
    guide <- ggproto(NULL, original, available_aes = union(guide$available_aes, unname(available)))
  }
  guide
}

ggchord_scope_scale <- function(scale, aliases, theme) {
  if (is.null(scale$palette)) {
    palette <- ggplot()$scales
    palette$add(scale)
    palette$set_palettes(theme)
  }
  scale$aesthetics <- unname(aliases)
  scale$guide <- ggchord_scope_guide(scale$guide, aliases)
  scale
}

# ggplot2 knows the mapping's type when it evaluates the start or Stat stage.
# Ask its ordinary scale registry for exactly that default using the public
# aesthetic, then install it under the layer's private aesthetic.
ggchord_defer_scale_defaults <- function(parent, channels, theme) {
  force(parent)
  force(channels)
  force(theme)
  ggproto(NULL, parent, add_defaults = function(self, data, env) {
    for (alias in intersect(setdiff(names(data), self$input()), names(channels))) {
      key <- unname(channels[[alias]])
      defaults <- ggproto(NULL, parent, scales = list())
      defaults$add_defaults(setNames(list(data[[alias]]), key), env)
      scale <- defaults$get_scales(key)
      if (!is.null(scale)) self$add(ggchord_scope_scale(scale, setNames(alias, key), theme))
    }
    ggproto_parent(parent, self)$add_defaults(data, env)
  })
}

ggchord_apply_scoped_scales <- function(plot, specs = ggchord_drawing_scopes(plot)) {
  if (!length(specs)) return(plot)
  plot$scales <- plot$scales$clone()
  aliases <- vector("list", length(plot$layers))
  deferred <- character()
  theme <- ggplot2::theme_get() + plot$theme
  guides <- plot$guides$guides
  for (spec in specs) {
    if (is.null(spec$scale)) {
      deferred[spec$aliases] <- spec$aesthetics
    } else {
      scale <- ggchord_prepare_domain_scale(spec$scale, plot)$clone()
      plot$scales$add(ggchord_scope_scale(scale, spec$aliases, theme))
    }
    # Preserve public guide settings and titles across the internal rename.
    for (aesthetic in spec$aesthetics) {
      for (i in sort(spec$selected)) {
        layer <- plot$layers[[i]]
        if (!is.null(layer$aes_params[[aesthetic]])) next
        mapping <- layer$mapping[[aesthetic]] %||%
          if (isTRUE(layer$inherit.aes)) plot$mapping[[aesthetic]]
        if (is.null(mapping)) next
        key <- spec$aliases[[aesthetic]]
        label <- layer$ggchord_params$track_mapping_labels[[aesthetic]]
        if (!is.null(label)) plot$labels[[key]] <- label
        if (isTRUE(spec$automatic)) {
          if (aesthetic %in% names(plot$labels)) plot$labels[key] <- plot$labels[aesthetic]
          guides[[key]] <- ggchord_scope_guide(guides[[aesthetic]], spec$aliases)
        }
        break
      }
    }
    for (i in spec$selected) {
      aliases[[i]] <- c(aliases[[i]], spec$aliases)
    }
  }
  if (length(deferred)) plot$scales <- ggchord_defer_scale_defaults(plot$scales, deferred, theme)
  for (i in which(lengths(aliases) > 0L))
    plot$layers[[i]] <- ggchord_scope_layer(plot$layers[[i]], aliases[[i]], plot$mapping)
  automatic <- Filter(function(spec) isTRUE(spec$automatic), specs)
  copied <- intersect(unlist(lapply(automatic, `[[`, "aesthetics")),
    union(names(guides), names(plot$labels)))
  for (key in copied) {
    shared <- any(vapply(plot$layers, function(layer) {
      !is.null(layer$mapping[[key]]) ||
        (isTRUE(layer$inherit.aes) && !is.null(plot$mapping[[key]]))
    }, logical(1)))
    if (!shared) {
      plot$labels[key] <- NULL
      guides[[key]] <- NULL
    }
  }
  parent_guides <- plot$guides
  plot$guides <- ggproto(NULL, parent_guides, guides = guides)
  plot
}

# Obtain the already-evaluated column without running user mappings twice.
ggchord_scale_column <- function(layer, aesthetic, role) {
  columns <- layer$ggchord_params$visual_columns
  column <- if (aesthetic %in% names(columns)) unname(columns[aesthetic]) else role
  if (is.na(column) || !column %in% names(layer$data)) return(NULL)
  layer$data[[column]]
}

ggchord_mapped_scale <- function(plot, layers, aesthetic, roles) {
  role <- unique(roles)
  if (length(role) != 1L) return(NULL)
  kinds <- vapply(layers, function(l) l$ggchord_type %||% "", character(1))
  feature <- kinds == "gene_poly" | vapply(layers, function(l)
    identical(l$ggchord_params$label_kind, "feature"), logical(1))
  if (role %in% c("accver", "strand", "type")) {
    discrete <- vapply(layers, function(layer) {
      value <- ggchord_scale_column(layer, aesthetic, role)
      is.character(value) || is.factor(value)
    }, logical(1))
    if (!all(discrete)) return(NULL)
    if (role == "accver" && all(kinds == "seq" | feature))
      return(scale_colour_sequence(aesthetics = aesthetic))
    if (role == "strand" && all(feature))
      return(scale_fill_strand(aesthetics = aesthetic))
    if (role == "type" && all(feature) && isTRUE(plot$coordinates$ggchord_circular))
      return(scale_fill_feature_type(aesthetics = aesthetic))
  }
  types <- if (aesthetic == "fill") c("ribbon", "fusion_ribbon") else c("link", "alignment")
  if (!identical(role, "pident") || !all(kinds %in% types)) return(NULL)
  name <- if (aesthetic %in% names(plot$labels)) waiver() else "Identity (%)"
  if (aesthetic == "fill") return(scale_fill_pident(name = name))
  key_width <- plot$ggchord$theme_explicit$legend.key.width %||%
    plot$ggchord$theme_explicit$legend.key.size %||% unit(5, "mm")
  scale_colour_gradientn(name = name, colours = c("#34457E", "#2FA96B"),
    limits = c(0, 100), breaks = c(0, 50, 80, 90, 95, 100),
    guide = guide_legend(order = 2, theme = theme(legend.key.width = key_width)))
}

ggchord_check_mixed_scale <- function(layers, aesthetic, roles) {
  kinds <- vapply(seq_along(layers), function(i) {
    value <- ggchord_scale_column(layers[[i]], aesthetic, roles[i])
    if (!length(value) || inherits(value, "AsIs")) return("")
    if (is.numeric(value)) return("continuous")
    if (is.character(value) || is.factor(value) || is.logical(value)) return("discrete")
    ""
  }, character(1))
  pair <- match(c("continuous", "discrete"), kinds)
  if (anyNA(pair)) return(invisible(NULL))
  describe <- function(i) {
    id <- layers[[i]]$ggchord_layer_id %||% layers[[i]]$ggchord_type %||% "unnamed"
    paste0("'", id, "' maps `", roles[i], "` (", kinds[i], ")")
  }
  id <- ggchord_declared_id(layers[[pair[2L]]]$ggchord_params$id) %||% "features"
  ggchord_stop("Incompatible mappings for shared `", aesthetic, "`: ",
    describe(pair[1L]), " and ", describe(pair[2L]), ". Give a layer an explicit id and use ",
    "scale_for(\"", id, "\", scale_", aesthetic, "_discrete()) for its independent scale.")
}

# GuideLegend already knows the raw break behind each swatch. Pass that
# strand to just the feature keys, after breaks, labels and reversal resolve.
ggchord_strand_guide <- function(guide) {
  if (identical(guide, "legend") || is_waiver(guide)) guide <- guide_legend()
  if (!inherits(guide, "GuideLegend")) return(guide)
  ggproto(NULL, guide, get_layer_key = function(params, layers, data, theme = NULL) {
    out <- guide$get_layer_key(params, layers, data, theme)
    if (!is.null(params$override.aes$strand)) return(out)
    for (i in seq_along(out$decor)) {
      layer <- layers[[i]]
      if (!identical(layer$ggchord_type, "gene_poly")) next
      mapped <- layer$computed_mapping[intersect(names(layer$computed_mapping), params$aesthetic)]
      if (any(vapply(mapped, function(q) identical(as_label(q), "strand"), logical(1)))) {
        out$decor[[i]]$data$strand <- as.character(params$key$.value)
      }
    }
    out
  })
}

ggchord_add_visual_scales <- function(plot) {
  # Use local copies for the plot-dependent reference palette and identity guide.
  plot$scales <- plot$scales$clone()
  for (i in seq_along(plot$scales$scales)) {
    s <- plot$scales$scales[[i]]
    if (!any(startsWith(s$aesthetics, "ggchord_scope_")) && !is.null(s$ggchord_palette))
      plot$scales$scales[[i]] <- ggchord_prepare_domain_scale(s, plot)
  }
  for (aesthetic in c("fill", "colour")) {
    if (plot$scales$has_scale(aesthetic)) next
    layers <- Filter(function(layer) {
      if (aesthetic %in% names(layer$aes_params)) return(FALSE)
      mapping <- layer$mapping[[aesthetic]]
      if (is.null(mapping) && isTRUE(layer$inherit.aes)) mapping <- plot$mapping[[aesthetic]]
      !is.null(mapping) && ggchord_trains_visual_scale(mapping)
    }, plot$layers)
    if (!length(layers)) next
    roles <- vapply(layers, function(layer) {
      mapping <- layer$mapping[[aesthetic]] %||% plot$mapping[[aesthetic]]
      as_label(mapping)
    }, character(1))
    if (length(unique(roles)) > 1L) ggchord_check_mixed_scale(layers, aesthetic, roles)
    scale <- ggchord_mapped_scale(plot, layers, aesthetic, roles)
    if (is.null(scale)) next
    scale <- ggchord_prepare_domain_scale(scale, plot)
    scale$ggchord_managed <- TRUE
    plot$scales$add(scale)
  }
  strand_aes <- unique(unlist(lapply(plot$layers, function(layer) {
    if (!identical(layer$ggchord_type, "gene_poly")) return(NULL)
    mapped <- names(Filter(function(q) identical(as_label(q), "strand"), layer$mapping))
    grep("(^|_)(colour|fill)$", mapped, value = TRUE)
  }), use.names = FALSE))
  if (length(strand_aes)) {
    guides <- plot$guides$guides
    for (key in strand_aes) {
      scale <- plot$scales$get_scales(key)
      guide <- guides[[key]] %||% scale$guide %||% "legend"
      guides[[key]] <- ggchord_strand_guide(guide)
    }
    parent <- plot$guides
    plot$guides <- ggproto(NULL, parent, guides = guides)
  }
  plot
}
