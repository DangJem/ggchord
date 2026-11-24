# Cross-view visual scopes train through the same native scale that will
# map the drawing layers. The recorder sees values after ggplot2's transform
# and Stat stages, before palette, OOB and after_scale.
ggchord_fusion_scale_probe <- function(scale, record) {
  out <- scale$clone()
  parent <- out
  ggproto(NULL, parent, train = function(self, x) {
    if (length(x)) record$values[[length(record$values) + 1L]] <- x
    ggproto_parent(parent, self)$train(x)
  })
}

ggchord_fusion_freeze_scale <- function(scale, aesthetics = scale$aesthetics) {
  if (!inherits(scale, c("ScaleDiscrete", "ScaleContinuous")) || inherits(scale, "ScaleBinned"))
    ggchord_stop("Fusion: this trained scale has no supported fixed domain")
  result <- scale$clone()
  result$aesthetics <- aesthetics
  result$limits <- scale$get_limits()
  result
}

ggchord_fusion_get <- function(plot, path, caller = "Fusion scale_for()") {
  for (id in path) {
    key <- ggchord_fusion_path_key(id)
    if (!inherits(plot, "ggchord_fusion") || !key %in% names(plot$panels))
      ggchord_stop(caller, ": unknown view path ", ggchord_fusion_path_key(path),
        "; candidates: ", paste(names(plot$panels), collapse = ", "))
    plot <- plot$panels[[key]]$plot
  }
  plot
}

ggchord_fusion_modify <- function(plot, path, modify) {
  if (!length(path)) return(modify(plot))
  key <- ggchord_fusion_path_key(path[1L])
  plot$panels[[key]]$plot <- ggchord_fusion_modify(plot$panels[[key]]$plot, path[-1L], modify)
  plot
}

ggchord_fusion_scale_targets <- function(plot, targets, scale) {
  paths <- lapply(targets, `[[`, "panel_path")
  keys <- vapply(paths, ggchord_fusion_path_key, character(1))
  ids <- vapply(targets, `[[`, "layer_id", FUN.VALUE = character(1))
  if (anyDuplicated(data.frame(keys, ids)))
    ggchord_stop("Fusion scale_for(): duplicate view and layer targets")
  groups <- split(seq_along(keys), match(keys, unique(keys)))
  units <- list()
  out <- lapply(groups, function(at) {
    path <- paths[[at[1L]]]
    view <- ggchord_fusion_get(plot, path)
    selected <- ids[at]
    if (inherits(view, "ggchord_fusion")) {
      layers <- view$ribbons[selected]
      if (any(!selected %in% names(view$ribbons)))
        ggchord_stop("Fusion scale_for(): unknown root Link in ", ggchord_fusion_path_key(path),
          ": ", paste(setdiff(selected, names(view$ribbons)), collapse = ", "))
      layers <- lapply(layers, function(layer)
        list(mapping = layer$mapping, data = layer$data))
    } else {
      genomic <- inherits(view, "ggchord")
      if (genomic) view <- ggchord_register_layers(view)
      if (!genomic && is.null(view$ggchord_alignment))
        ggchord_stop("Fusion scale_for(): ordinary ggplot does not declare a supported source layer")
      known <- vapply(view$layers, function(layer) layer$ggchord_layer_id %||% "", character(1))
      if (any(!selected %in% known))
        ggchord_stop("Fusion scale_for(): unknown layer: ", paste(setdiff(selected, known), collapse = ", "))
      layers <- lapply(view$layers[match(selected, known)], function(layer) {
        if (is.null(layer$ggchord_type) && !inherits(layer, "LayerResidue"))
          ggchord_stop("Fusion scale_for(): unsupported native layer '", layer$ggchord_layer_id, "'")
        mapping <- layer$ggchord_input_mapping %||% layer$mapping
        if (isTRUE(layer$inherit.aes))
          mapping <- c(mapping, view$mapping[setdiff(names(view$mapping), names(mapping))])
        data <- layer$ggchord_input_data
        if (is.null(data) && inherits(layer, "LayerResidue")) {
          data <- if (is.data.frame(layer$data)) layer$data else view$data
        }
        list(mapping = mapping, data = data)
      })
    }
    for (layer in layers) for (key in intersect(names(layer$mapping), scale$aesthetics)) {
      field <- all.vars(get_expr(layer$mapping[[key]]))
      if (length(field) == 1L && is.data.frame(layer$data)) {
        column <- layer$data[[field]]
        units[[key]] <<- union(units[[key]], as.character(attr(column, "unit") %||% attr(column, "units")))
      }
    }
    list(path = path, ids = selected)
  })
  if (any(lengths(units) > 1L)) ggchord_stop("Fusion scale_for(): selected mappings declare conflicting units")
  unname(out)
}

ggchord_fusion_add_scope <- function(plot, groups, scale) {
  for (group in groups) {
    scope <- scale_for(scale, layers = group$ids)
    plot <- ggchord_fusion_modify(plot, group$path, function(view) {
      if (inherits(view, "ggchord_fusion")) {
        view$scope_scales <- c(view$scope_scales, list(scope))
      } else if (!inherits(view, "ggchord")) {
        view$ggchord_fusion_native_scopes <- c(view$ggchord_fusion_native_scopes, list(scope))
      } else {
        view$ggchord$scoped_scales <- c(view$ggchord$scoped_scales, list(scope))
        view$ggchord_fusion_level_lock <- TRUE
      }
      view
    })
  }
  plot
}

ggchord_fusion_train_scales <- function(plot, on_error, path = character()) {
  declare <- function(scope) {
    groups <- ggchord_fusion_scale_targets(plot, scope$targets, scope$scale)
    scale <- scope$scale$clone()
    if (!inherits(scale, c("ScaleContinuous", "ScaleDiscrete")) || inherits(scale, "ScaleBinned"))
      ggchord_stop("Fusion scale_for(): shared visual scale must be native continuous or discrete")
    if (is.null(scale$palette)) {
      styles <- Filter(function(x) inherits(x, "theme"), plot$styles)
      theme <- Reduce(`+`, styles, init = ggplot2::theme_get())
      registry <- ggplot()$scales
      registry$add(scale)
      registry$set_palettes(theme)
    }
    record <- new.env(parent = emptyenv())
    probe <- ggchord_fusion_add_scope(plot, groups, ggchord_fusion_scale_probe(scale, record))
    for (group in groups) {
      view <- ggchord_fusion_get(probe, group$path)
      if (inherits(view, "ggchord_fusion")) ggchord_fusion_build_styles(view) else ggchord_fusion_build_leaf(view)
    }
    for (values in record$values) scale$train(values)
    locked <- ggchord_fusion_freeze_scale(scale)
    plot <<- ggchord_fusion_add_scope(plot, groups, locked)
    plot$fusion_trained_scales[[length(plot$fusion_trained_scales) + 1L]] <<-
      list(panel_path = character(), targets = scope$targets, aesthetics = scale$aesthetics,
        limits = locked$limits, transform = if (!scale$is_discrete()) scale$get_transformation()$name,
        discrete = scale$is_discrete(), observations = sum(lengths(record$values)))
  }
  for (scope in plot$fusion_shared_scales) tryCatch(declare(scope), error = function(e) {
    if (is.null(on_error)) ggchord_rethrow(e)
    on_error(path, e)
  })
  for (key in names(plot$panels)) {
    child <- plot$panels[[key]]$plot
    if (inherits(child, "ggchord_fusion"))
      plot$panels[[key]]$plot <- ggchord_fusion_train_scales(child, on_error,
        c(path, ggchord_fusion_path_ids(key)))
  }
  plot
}

ggchord_fusion_build_leaf <- function(plot) {
  scopes <- plot$ggchord_fusion_native_scopes
  if (length(scopes))
    plot <- ggchord_apply_scoped_scales(plot, ggchord_scale_scopes(plot, scopes))
  ggplot_build(plot)
}

# Cache source work for the lifetime of a composition build. In a derived
# window, copy only the prepared fields of unchanged layers, keeping the
# detail's own export capture and geometry/guide state.
ggchord_fusion_prepare_sources <- function(plot, on_error, path = character()) {
  for (id in names(plot$panels)) {
    panel <- plot$panels[[id]]
    child_path <- c(path, ggchord_fusion_path_ids(id))
    view <- panel$plot
    prepared <- tryCatch({
      if (inherits(view, "ggchord_fusion")) {
        ggchord_fusion_prepare_sources(view, on_error, child_path)
      } else if (inherits(view, "ggchord")) {
        source <- if (length(panel$derived_from) == 1L)
          plot$panels[[ggchord_fusion_path_key(panel$derived_from)]]$plot
        view <- ggchord_render_plot_copy(view)
        if (!is.null(source)) view <- ggchord_magnify_reuse_source(view, source)
        ggchord_render_track_source(ggchord_render_inputs(view))
      } else view
    }, error = function(e) {
      if (is.null(on_error)) ggchord_rethrow(e)
      on_error(child_path, e)
      view
    })
    plot$panels[[id]]$plot <- prepared
  }
  plot
}

ggchord_magnify_reuse_source <- function(detail, source) {
  before <- detail$ggchord$magnify_parent
  if (is.null(before) || !identical(before$mapping, detail$mapping)) return(detail)
  for (i in seq_along(detail$layers)) {
    layer <- detail$layers[[i]]
    id <- layer$ggchord_layer_id
    if (!id %in% before$unchanged_layers) next
    from <- Filter(function(l) identical(l$ggchord_layer_id, id), source$layers)
    if (length(from) != 1L) next
    from <- from[[1L]]
    copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    stat <- from$ggchord_params$source_stat_cache
    if (!is.null(stat)) copy$ggchord_params$source_stat_cache <- stat
    if (!is.null(from$ggchord_params$link_stat_cache))
      copy$ggchord_params$link_stat_cache <- from$ggchord_params$link_stat_cache
    reader <- from$ggchord_params$render_source_reader
    if (is.function(reader)) {
      copy$ggchord_params$render_source_reader <- reader
      copy$ggchord_input_data <- reader
    }
    same_y <- identical(before$position_scales, detail$ggchord$position_scales) &&
      identical(Filter(function(x) isTRUE(x$y), before$scopes),
        Filter(function(x) isTRUE(x$y), detail$ggchord$scoped_scales))
    if (identical(layer$ggchord_type, "annotation") ||
        (isTRUE(from$ggchord_params$native_prepared) && same_y)) {
      capture <- layer$ggchord_params$stage_capture
      copy <- reconstruct_layer(from, from$data, mapping = from$mapping)
      copy$ggchord_params$stage_capture <- capture
    }
    detail$layers[[i]] <- copy
  }
  detail
}

ggchord_fusion_auto_scales <- function(plot, on_error, path = character()) {
  prepared <- list()
  for (id in names(plot$panels)) {
    panel <- plot$panels[[id]]
    view <- panel$plot
    child_path <- c(path, ggchord_fusion_path_ids(id))
    if (inherits(view, "ggchord_fusion")) {
      view <- ggchord_fusion_auto_scales(view, on_error, child_path)
    } else if (length(panel$derived_from) == 1L) {
      parent_id <- ggchord_fusion_path_key(panel$derived_from)
      view <- tryCatch({
        if (is.null(prepared[[parent_id]]))
          prepared[[parent_id]] <- ggplot_build(plot$panels[[parent_id]]$plot)
        ggchord_magnify_inherit_scales(view, prepared[[parent_id]])
      }, error = function(e) {
        if (is.null(on_error)) ggchord_rethrow(e)
        on_error(child_path, e)
        view$ggchord$magnify_scale_failure <- conditionMessage(e)
        view
      })
    }
    plot$panels[[id]]$plot <- view
  }
  plot
}

ggchord_magnify_inherit_scales <- function(detail, built) {
  before <- detail$ggchord$magnify_parent
  scopes <- detail$ggchord$scoped_scales
  scope_unchanged <- vapply(scopes, function(s)
    any(vapply(before$scopes, identical, logical(1), y = s)), logical(1))
  detail_layers <- vapply(detail$layers, function(l) l$ggchord_layer_id %||% "", character(1))
  visuals <- c("fill", "colour", "alpha", "linewidth", "linetype", "shape", "size", "stroke", "fill_from", "fill_to")
  ids_for <- function(scope) {
    if (!is.null(scope$layers)) return(scope$layers)
    detail_layers[vapply(detail$layers, function(layer)
      (layer$ggchord_params$track_id %||% "") %in% scope$tracks, logical(1))]
  }
  add <- list()
  status <- list()
  for (scale in built$plot$scales$scales) {
    for (aesthetic in scale$aesthetics) {
      key <- sub("^ggchord_scope_[0-9]+_", "", aesthetic)
      if (!key %in% visuals) next
      ids <- vapply(Filter(function(layer) {
        trained <- layer$mapping[[aesthetic]]
        !is.null(trained) && ggchord_trains_visual_scale(trained) &&
          !aesthetic %in% names(layer$aes_params) && layer$ggchord_layer_id %in% detail_layers
      }, built$plot$layers), function(layer) layer$ggchord_layer_id, character(1))
      if (!length(ids)) next
      same_scale <- !key %in% before$scale_changed
      same_preset <- identical(detail$ggchord$palette_preset, before$preset)
      explicit <- which(vapply(scopes, function(s) !isTRUE(s$y) && key %in% s$scale$aesthetics, logical(1)))
      overridden <- unique(unlist(lapply(scopes[explicit[!scope_unchanged[explicit]]], ids_for)))
      inherit <- intersect(setdiff(unique(ids), overridden), before$unchanged_layers)
      if (!same_scale || !same_preset || !identical(before$mapping, detail$mapping)) inherit <- character()
      status[[length(status) + 1L]] <- list(aesthetic = key, layers = ids,
        inherited = inherit, overridden = setdiff(ids, inherit), limits = scale$get_limits())
      if (!length(inherit)) next
      locked <- ggchord_fusion_freeze_scale(scale, key)
      # A trained alias uses the original visual aesthetic in the receiving
      # plot; the scope machinery applies its new private alias.
      guide <- locked$guide
      if (inherits(guide, "Guide") && !is.null(guide$available_aes))
        locked$guide <- ggproto(NULL, guide, available_aes = union(guide$available_aes, key))
      add[[length(add) + 1L]] <- scale_for(locked, layers = inherit)
    }
  }
  # Keep declarations for members outside a newly inherited domain.
  for (i in which(scope_unchanged)) {
    scope <- scopes[[i]]
    if (isTRUE(scope$y)) next
    affected <- Filter(function(s) length(intersect(s$scale$aesthetics, scope$scale$aesthetics)), add)
    if (!length(affected)) next
    if (length(scope$scale$aesthetics) > 1L)
      ggchord_stop("Auto Magnify: multi-aesthetic scoped scales require Manual")
    ids <- setdiff(ids_for(scope), unique(unlist(lapply(affected, `[[`, "layers"))))
    scopes[i] <- list(if (length(ids)) scale_for(scope$scale, layers = ids))
  }
  detail$ggchord$scoped_scales <- c(Filter(Negate(is.null), scopes), add)
  layout <- built$plot$ggchord$ref$layout
  same_y <- identical(before$position_scales, detail$ggchord$position_scales) &&
    identical(Filter(function(x) isTRUE(x$y), before$scopes),
      Filter(function(x) isTRUE(x$y), detail$ggchord$scoped_scales))
  if (same_y) for (id in names(layout$track_scales)) {
    selected <- layout$track_registry$track_id == id
    trained <- any(layout$track_registry$domain_status[selected] == "trained")
    scale <- layout$track_scales[[id]]
    members <- vapply(Filter(function(layer) identical(layer$ggchord_params$track_id, id),
      detail$layers), function(layer) layer$ggchord_layer_id, character(1))
    if (all(members %in% before$unchanged_layers))
      detail$ggchord$shared_track_scales[[id]] <- list(scale = scale, has_domain = trained)
  }
  detail$ggchord$magnify_scale_status <- list(visual = status,
    track_y = if (same_y) "inherited" else "overridden", statistic = "source_domain")
  detail
}
