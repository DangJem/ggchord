# A discrete legend is already trained when this wrapper sees its layers.
ggchord_legend_role <- function(layer) {
  switch(layer$ggchord_type %||% "",
    seq = "seq", gene_poly = layer$ggchord_params$feature_role %||% "feature",
    restriction_site = "site", link = "link", ribbon = "ribbon", "")
}

ggchord_legend_variable <- function(mapping, columns) {
  expr <- get_expr(ggchord_start_mapping(mapping))
  field <- if (is_symbol(expr)) as_string(expr) else
    if (is_call(expr, c("$", "[["), n = 2L) && is_symbol(expr[[2L]], ".data")) {
      key <- expr[[3L]]
      if (is_string(key)) key else if (is_call(expr, "$") && is_symbol(key)) as_string(key)
    }
  if (length(field) && field %in% columns) field else NA_character_
}

ggchord_legend_default_group <- function(layer, mapping, matched, default_breaks) {
  if (!all(matched %in% names(default_breaks))) return(NULL)
  data <- layer$ggchord_params$visual_source$data %||% layer$data
  fields <- unique(vapply(mapping[matched], ggchord_legend_variable, character(1), names(data)))
  breaks <- unique(unname(default_breaks[matched]))
  if (length(fields) == 1L && !is.na(fields) && length(breaks) == 1L)
    rlang::hash(list(fields, breaks[[1L]]))
}

ggchord_legend_layer <- function(layer, matched, sources, default_breaks) {
  role <- ggchord_legend_role(layer)
  default_geom <- switch(role, seq = GeomChordSeq, gene = feature_geom,
    feature = feature_geom, primer = feature_geom, site = GeomSequenceSite,
    link = GeomChordLinkLine, ribbon = GeomRibbonFlow)
  geom <- layer$ggchord_key_geom %||% layer$geom
  # Compare the underlying method. Plot preparation can wrap the same Geom
  # to supply defaults, and an explicit key_glyph installs its own method.
  builtin <- !is.null(default_geom) &&
    identical(environment(geom$draw_key)$f, environment(default_geom$draw_key)$f)
  mapping <- c(layer$computed_mapping, layer$stat$default_aes)
  matched <- matched[vapply(mapping[matched], ggchord_trains_visual_scale, logical(1))]
  ids <- sources[matched]
  source <- if (length(ids) && all(lengths(ids) == 1L)) unique(unlist(ids, use.names = FALSE))
  default_group <- if (builtin && length(source))
    ggchord_legend_default_group(layer, mapping, matched, default_breaks)
  head_length <- if ("head_length" %in% layer$ggchord_params$head_args)
    layer$ggchord_params$head_length_mm
  aliases <- layer$ggchord_key_aliases
  styled <- c(matched, names(layer$aes_params))
  styled <- c(styled, names(aliases)[aliases %in% styled])
  list(role = role, builtin = builtin, source = source, default_group = default_group,
    aliases = aliases, route = layer$ggchord_params$route$route,
    alpha_fixed = any(c("alpha", "ribbon_alpha", "feature_alpha") %in% styled),
    site_shape_fixed = isTRUE(layer$ggchord_params$fixed_shape),
    head_length = head_length)
}

ggchord_composite_guide <- function(parent, sources, explicit, default_breaks) {
  force(parent)
  force(sources)
  force(explicit)
  force(default_breaks)
  ggproto(NULL, parent,
    get_layer_key = function(params, layers, data, theme = NULL) {
      domain <- vapply(layers, function(layer) nzchar(ggchord_legend_role(layer)), logical(1))
      if (!any(domain)) return(parent$get_layer_key(params, layers, data, theme))
      # At this point ggplot2 passes only the matched trained aesthetics and
      # .id to compute_geom_2, before adding defaults and fixed aesthetics.
      matched <- vector("list", length(layers))
      observed <- lapply(seq_along(layers), function(i) {
        layer <- layers[[i]]
        ggproto(NULL, layer, compute_geom_2 = function(data, ...) {
          matched[[i]] <<- setdiff(names(data), ".id")
          layer$compute_geom_2(data, ...)
        })
      })
      names(observed) <- names(layers)
      out <- parent$get_layer_key(params, observed, data, theme)
      if (!length(out$decor) || length(out$decor) != length(layers)) return(out)
      contributions <- Map(function(layer, aes)
        ggchord_legend_layer(layer, aes, sources, default_breaks),
        layers, matched)
      for (i in seq_along(out$decor)) {
        aliases <- contributions[[i]]$aliases
        # A public override must replace the private column consumed by
        # the draw_key wrapper, even when the guide merges several aliases.
        for (field in intersect(names(aliases), names(out$override.aes))) {
          alias <- aliases[[field]]
          if (!alias %in% names(out$override.aes))
            out$decor[[i]]$data[[alias]] <- out$decor[[i]]$data[[field]]
        }
      }
      draws <- lapply(out$decor, function(part) part$data$.draw %||% rep(TRUE, nrow(out$key)))
      builtin <- vapply(contributions, `[[`, logical(1), "builtin")
      out$ggchord_composite <- (any(builtin) && any(Reduce(`|`, draws[builtin]), na.rm = TRUE)) ||
        any(Reduce(`+`, draws) > 1L & Reduce(`|`, draws[domain]), na.rm = TRUE)
      if (out$ggchord_composite) {
        for (i in seq_along(contributions)) {
          aliases <- contributions[[i]]$aliases
          override <- names(out$override.aes)
          contributions[[i]]$alpha_fixed <- contributions[[i]]$alpha_fixed ||
            any(c("alpha", "ribbon_alpha", aliases[c("alpha", "ribbon_alpha")]) %in% override)
          contributions[[i]]$site_shape_fixed <- contributions[[i]]$site_shape_fixed ||
            any(c("shape", aliases["shape"]) %in% override)
        }
        out$ggchord_key_layers <- contributions
      }
      out
    },
    build_decor = function(decor, grobs, elements, params) {
      if (!isTRUE(params$ggchord_composite))
        return(parent$build_decor(decor, grobs, elements, params))
      ggchord_composite_keys(decor, elements, params, explicit)
    },
    measure_grobs = function(grobs, params, elements) {
      if (isTRUE(params$ggchord_composite)) {
        elements$width_cm <- attr(grobs$decor[[1L]], "width")
        elements$height_cm <- attr(grobs$decor[[1L]], "height")
      }
      parent$measure_grobs(grobs, params, elements)
    })
}

# Compare only package glyphs. Grid-generated names carry no visual information.
ggchord_legend_signature <- function(glyph) {
  glyph$name <- NULL
  for (field in c("col", "fill")) {
    colour <- glyph$gp[[field]]
    if (is.character(colour)) glyph$gp[[field]] <- unname(grDevices::col2rgb(colour, alpha = TRUE))
  }
  if (inherits(glyph, "gTree")) {
    if (length(glyph$children))
      glyph$children <- lapply(unname(glyph$children), ggchord_legend_signature)
    glyph$childrenOrder <- NULL
  }
  glyph
}

ggchord_legend_parts <- function(decor, layers, row) {
  parts <- list()
  signatures <- new.env(parent = emptyenv())
  for (i in seq_along(decor)) {
    dec <- decor[[i]]
    data <- dec$data[row, , drop = FALSE]
    if (!isTRUE(data$.draw %||% TRUE)) next
    part <- layers[[i]]
    part$draw <- dec$draw_key
    part$params <- dec$params
    part$data <- data
    if (part$builtin) {
      part$style <- ggchord_alias_aes(data, setNames(names(part$aliases), part$aliases))
      part$site_marker <- part$role == "site" && !part$site_shape_fixed &&
        identical(part$params$site_key_shape, "tick") &&
        !any(part$style$shape %in% c("tick", "point", "cut"))
      if (part$role %in% c("gene", "feature", "primer")) {
        part$shape <- ggchord_feature_key_shape(part$style, part$params)
        if (identical(part$shape, "arrow"))
          part$params$arrow_key$length <- min(.55, .16 * (part$head_length %||% 2) / 2)
      }
      part$glyph <- part$draw(data, part$params, ggchord_legend_artwork_mm)
    }
    signature <- if (part$builtin && length(part$source) == 1L)
      list(part$source, part$role, part$route, part$shape, isTRUE(part$params$flow),
        if (part$site_marker) part$style$size %||% 2,
        part$alpha_fixed && part$role %in% c("gene", "feature", "ribbon"),
        ggchord_legend_signature(part$glyph))
    if (!is.null(signature)) {
      hash <- rlang::hash(signature)
      possible <- signatures[[hash]]
      if (any(vapply(possible, identical, logical(1), signature))) next
      signatures[[hash]] <- c(possible, list(signature))
    }
    parts[[length(parts) + 1L]] <- part
  }
  parts
}

ggchord_legend_requested_size <- function(params, explicit) {
  fields <- c("legend.key.width", "legend.key.height")
  guide_theme <- params$theme
  local <- lapply(fields, function(name) guide_theme[[name]] %||% guide_theme$legend.key.size)
  global <- lapply(fields, function(name) explicit[[name]] %||% explicit$legend.key.size)
  Map(`%||%`, local, global)
}

ggchord_legend_key_size <- function(params, explicit) {
  requested <- ggchord_legend_requested_size(params, explicit)
  stretch <- vapply(requested, function(x) !is.null(x) && grid::unitType(x) == "null", logical(1))
  size <- c(8, 6.4)
  measure <- list(convertWidth, convertHeight)
  for (i in which(lengths(requested) > 0L & !stretch))
    size[i] <- measure[[i]](requested[[i]], "mm", TRUE)
  list(size = size, unlocked = lengths(requested) == 0L, stretch = stretch)
}

ggchord_composite_keys <- function(decor, elements, params, explicit) {
  dimensions <- ggchord_legend_key_size(params, explicit)
  size <- dimensions$size
  parts <- lapply(seq_len(params$n_breaks), function(i)
    ggchord_legend_parts(decor, params$ggchord_key_layers, i))
  build <- function(size) lapply(parts, function(parts) {
    plan <- ggchord_legend_plan(parts, size)
    glyph <- ggchord_legend_render(plan)
    list(glyph = glyph, crowded = plan$crowded || ggchord_legend_crowded(glyph))
  })
  layouts <- build(size)
  expanded <- ifelse(dimensions$unlocked, c(9.5, 7.6), size)
  if (is.null(params$ggchord_key_expand)) {
    canvas <- ggchord_render_dimensions()
    allow_expansion <- all(is.finite(canvas)) && all(expanded <= .2 * canvas)
  } else allow_expansion <- isTRUE(params$ggchord_key_expand)
  if (allow_expansion && any(expanded > size) &&
      any(vapply(layouts, `[[`, logical(1), "crowded"))) {
    larger <- build(expanded)
    if (sum(vapply(larger, `[[`, logical(1), "crowded")) <
        sum(vapply(layouts, `[[`, logical(1), "crowded"))) {
      layouts <- larger
      size <- expanded
    }
  }
  just <- elements$key_just %||% c(.5, .5)
  lapply(layouts, function(layout) {
    key <- grobTree(elements$key, layout$glyph, vp = viewport(x = just[1], y = just[2],
      just = just, width = if (dimensions$stretch[1]) unit(1, "npc") else unit(size[1], "mm"),
      height = if (dimensions$stretch[2]) unit(1, "npc") else unit(size[2], "mm"),
      clip = "on"))
    # These measurements describe the allocated viewport, not the nominal
    # size of its children. The label cell may be taller than the visible key.
    attr(key, "width") <- size[1] / 10
    attr(key, "height") <- size[2] / 10
    attr(key, "ggchord_key_crowded") <- layout$crowded
    key
  })
}

ggchord_legend_crowded <- function(grob) {
  isTRUE(attr(grob, "ggchord_key_crowded")) ||
    any(vapply(grob$grobs %||% grob$children, ggchord_legend_crowded, logical(1)))
}

ggchord_warn_composite_keys <- function(table) {
  boxes <- startsWith(table$layout$name, "guide-box")
  if (any(vapply(table$grobs[boxes], ggchord_legend_crowded, logical(1)))) ggchord_legend_warning()
}

ggchord_legend_warning <- function() {
  warning(
    "Composite legend key is too small to display all glyphs clearly; increase legend.key.width/height or simplify the legend.",
    call. = FALSE)
}
