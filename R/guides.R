# guides.R - thin ggplot2 guide wrappers

#' Scale legend furniture relative to the active export device
#' @noRd
ggchord_device_scale <- function(device_inches = NULL,
                                 reference = c(width = 8, height = 6),
                                 limits = c(0.5, 1.6)) {
  if (is.null(device_inches)) {
    available <- ggchord_render_context()$available_canvas_mm
    device_inches <- if (length(available) == 2L && all(is.finite(available))) {
      available / 25.4
    } else if (dev.cur() == 1L) {
      c(NA_real_, NA_real_)
    } else tryCatch(
      dev.size("in"), error = function(e) c(NA_real_, NA_real_)
    )
  }
  if (!is.numeric(device_inches) || length(device_inches) != 2L ||
      any(!is.finite(device_inches)) || any(device_inches <= 0)) {
    return(1)
  }
  value <- min(device_inches / unname(reference))
  min(limits[2], max(limits[1], value))
}

ggchord_check_guide_scale <- function(size_scale, caller) {
  if (!is.numeric(size_scale) || length(size_scale) != 1L ||
      !is.finite(size_scale) || size_scale <= 0) {
    ggchord_stop(caller, ": size_scale must be one positive finite number")
  }
  size_scale
}

#' @noRd
guide_ggchord_legend <- function(
    title = ggplot2::waiver(), theme = NULL, position = NULL, direction = NULL,
    override.aes = list(), nrow = NULL, ncol = NULL, reverse = FALSE,
    order = 0, size_scale = 1, ...) {
  size_scale <- ggchord_check_guide_scale(
    size_scale, "guide_ggchord_legend()"
  )
  compact <- ggplot2::theme(
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.height = unit(4 * size_scale, "mm"),
    legend.key.width = unit(5.5 * size_scale, "mm"),
    legend.text = element_text(
      size = 8 * size_scale,
      margin = margin(l = 1.2 * size_scale)
    ),
    legend.title = element_text(size = 9 * size_scale),
    legend.margin = margin(
      t = size_scale, b = size_scale, unit = "mm"
    )
  )
  if (!is.null(theme)) compact <- compact + theme
  guide_legend(
    title = title, theme = compact, position = position,
    direction = direction, override.aes = override.aes,
    nrow = nrow, ncol = ncol, reverse = reverse, order = order, ...
  )
}

#' @noRd
guide_ggchord_colourbar <- function(
    title = ggplot2::waiver(), theme = NULL, nbin = NULL, display = "raster",
    alpha = NA, draw.ulim = TRUE, draw.llim = TRUE, angle = NULL,
    position = NULL, direction = NULL, reverse = FALSE, order = 0,
    available_aes = c(
      "colour", "color", "fill", "ribbon_fill", "gene_fill",
      "feature_fill", "region_fill"
    ), size_scale = 1, ...) {
  size_scale <- ggchord_check_guide_scale(
    size_scale, "guide_ggchord_colourbar()"
  )
  horizontal <- identical(direction, "horizontal") ||
    (!is.null(position) && position %in% c("top", "bottom"))
  compact <- ggplot2::theme(
    legend.title.position = "top",
    legend.key.width = unit(
      (if (horizontal) 50 else 3.6) * size_scale, "mm"
    ),
    legend.key.height = unit(
      (if (horizontal) 3.6 else 50) * size_scale, "mm"
    ),
    legend.ticks.length = unit(size_scale, "mm"),
    legend.text = element_text(
      size = 8 * size_scale,
      margin = margin(l = 1.2 * size_scale)
    ),
    legend.title = element_text(
      size = 9 * size_scale,
      margin = margin(b = 2 * size_scale, unit = "mm")
    ),
    legend.margin = margin(
      t = size_scale, b = size_scale, unit = "mm"
    )
  )
  if (!is.null(theme)) compact <- compact + theme
  native <- guide_colourbar(
    title = title, theme = compact, nbin = nbin, display = display,
    alpha = alpha, draw.ulim = draw.ulim, draw.llim = draw.llim,
    angle = angle, position = position, direction = direction,
    reverse = reverse, order = order, available_aes = available_aes, ...
  )
  ggchord_colourbar_dimensions(native, theme, size_scale)
}

# Resolve automatic dimensions after ggplot2 knows the final guide direction.
# A constructor cannot infer a later theme(legend.position = "bottom").
ggchord_colourbar_dimensions <- function(native, user_theme, size_scale) {
  override <- native$override_elements
  ggproto(NULL, native,
    override_elements = function(params, elements, theme) {
      elements <- override(params, elements, theme)
      horizontal <- identical(params$direction, "horizontal")
      long <- if(horizontal) "width" else "height"
      short <- if(horizontal) "height" else "width"
      measured <- elements$text
      labels <- params$key$.label
      positions <- as.numeric(params$key$.value)
      if(length(positions) && any(abs(positions)>1,na.rm=TRUE))
        positions <- positions / params$nbin
      required <- 0
      if(length(labels)>1 && !inherits(measured,"element_blank")) {
        extents <- vapply(seq_along(labels),function(i) {
          g <- textGrob(labels[i],rot=measured@angle %||% 0,
            gp=gpar(fontsize=measured@size %||% 8,
              fontfamily=measured@family %||% "",fontface=measured@face %||% "plain"))
          if(horizontal) convertWidth(grobWidth(g),"mm",TRUE) else
            convertHeight(grobHeight(g),"mm",TRUE)
        },numeric(1))
        order <- order(positions); gaps <- diff(positions[order]); valid <- is.finite(gaps) & gaps>0
        needs <- ((head(extents[order],-1)+tail(extents[order],-1))/2+1.2)/gaps
        if(any(valid))required <- max(needs[valid],na.rm=TRUE)
      }
      for(axis in c(long,short)) {
        requested <- user_theme[[paste0("legend.key.",axis)]] %||% user_theme$legend.key.size
        if(is.null(requested)) {
          value <- if(axis==long) max(50*size_scale,required) else 3.6*size_scale
          elements[[paste0("key_",axis)]] <- unit(value,"mm")
          elements[[paste0(axis,"_cm")]] <- value/10
          elements[[paste0("stretch_",if(axis=="width")"x" else "y")]] <- FALSE
        }
      }
      elements
    })
}

#' @noRd
guide_ggchord_colorbar <- guide_ggchord_colourbar

#' Resolve a guide theme and placement for one chord role
#' @noRd
ggchord_role_guide_spec <- function(plot, role, colourbar = FALSE) {
  # An explicitly placed guide or an ordinary plot theme may set the outer
  # side. Keep the default Identity placement and device-aware sizing.
  global_position <- plot$theme$legend.position
  if (identical(global_position, "none")) return(list(hidden = TRUE))
  default_position <- if (identical(role, "ribbon")) "left" else "right"
  resolved_position <- global_position %||% default_position
  resolved_direction <- plot$theme$legend.direction
  size_scale <- if (is.null(plot$ggchord$layout_policy)) ggchord_device_scale() else
    plot$ggchord$adaptive_report$spacing_factor %||% 1
  explicit <- plot$ggchord$theme_explicit
  vals <- list(
    legend.key.width = explicit$legend.key.width %||% explicit$legend.key.size %||%
      plot$theme$legend.key.width,
    legend.key.height = explicit$legend.key.height %||% explicit$legend.key.size %||%
      plot$theme$legend.key.height
  )
  horizontal <- colourbar && (
    identical(resolved_direction, "horizontal") ||
      resolved_position %in% c("top", "bottom") ||
      identical(plot$theme$legend.box, "horizontal")
  )
  if (is.null(vals$legend.key.width))
    vals$legend.key.width <- unit(
      if (colourbar && horizontal) 54 else if (colourbar) 3.8 else 5.5, "mm"
    ) * size_scale
  if (is.null(vals$legend.key.height))
    vals$legend.key.height <- unit(
      if (colourbar && !horizontal) 80 else if (colourbar) 3.8 else 4, "mm"
    ) * size_scale
  # Preserve inherited typography while resolving relative sizes for the
  # device measurement. Absolute sizes and blank elements remain explicit.
  text_fields <- c(legend.text = 8, legend.title = 9)
  for (field in names(text_fields)) {
    style <- plot$theme[[field]] %||% element_text()
    if (inherits(style, "element_text") &&
        (is.null(style@size) || inherits(style@size, "rel"))) {
      style <- ggchord_theme_element(plot, field) %||% style
      if (inherits(style, "element_text"))
        style@size <- ggchord_theme_point_size(plot, field, text_fields[[field]]) * size_scale
    }
    vals[[field]] <- style
  }
  list(
    hidden = FALSE,
    position = resolved_position,
    direction = resolved_direction %||% NULL,
    theme = if (length(vals)) do.call(theme, vals) else NULL,
    size_scale = size_scale,
    height_budget = ggchord_guide_height_budget(plot,resolved_position)
  )
}

#' Construct the default guide for one chord role
#' @noRd
ggchord_role_guide <- function(plot, role, colourbar = FALSE, order = 0,
                               override.aes = list(), labels = NULL) {
  spec <- ggchord_role_guide_spec(plot, role, colourbar)
  if (isTRUE(spec$hidden)) return("none")
  if (colourbar) {
    guide_ggchord_colourbar(
      position = spec$position, direction = spec$direction,
      theme = spec$theme, size_scale = spec$size_scale, order = order,
      available_aes = "ribbon_fill"
    )
  } else {
    columns <- ggchord_guide_columns(labels, spec)
    guide_ggchord_legend(
      ncol = columns, position = spec$position, direction = spec$direction,
      theme = spec$theme, size_scale = spec$size_scale, order = order,
      override.aes = override.aes
    )
  }
}

ggchord_identity_bar_width <- function() unit(2.8, "mm")

ggchord_identity_title_gap <- function(guide) {
  ggproto(NULL, guide, setup_elements = function(params, elements, theme) {
    elements <- guide$setup_elements(params, elements, theme)
    title <- elements$title
    supplied <- theme$legend.title
    if (!inherits(title, "element_text") ||
        (inherits(supplied, "element_text") && !is.null(supplied@margin)))
      return(elements)
    # Endpoint labels are centred on the edge of the bar and extend past it.
    position <- elements$title_position
    facing <- match(position, c("bottom", "left", "top", "right"))
    along <- (params$direction == "horizontal") == (position %in% c("left", "right"))
    text <- elements$text
    endpoint <- 0
    if (along && inherits(text, "element_text") && length(params$key$.label)) {
      label <- params$key$.label[c(1L, length(params$key$.label))]
      measure <- if (params$direction == "horizontal") grobWidth else grobHeight
      endpoint <- max(vapply(seq_along(label), function(i) {
        glyph <- textGrob(label[i], rot = params$angle %||% text@angle %||% 0,
          gp = gpar(fontsize = text@size, fontfamily = text@family, fontface = text@face))
        convertUnit(measure(glyph), "mm", valueOnly = TRUE)
      }, numeric(1)))
    }
    gap <- convertUnit(title@margin %||% margin(), "mm", valueOnly = TRUE)
    gap[facing] <- max(gap[facing], endpoint / 2 + 1.2)
    title@margin <- margin(gap[1L], gap[2L], gap[3L], gap[4L], unit = "mm")
    elements$title <- title
    elements
  })
}

# Unequal numeric intervals occupy their actual share of the 0..100 domain.
ggchord_identity_guide <- function(plot = NULL, aesthetic = "ribbon_fill") {
  spec <- if (!is.null(plot)) ggchord_role_guide_spec(plot,"ribbon",TRUE) else list()
  if (isTRUE(spec$hidden)) return("none")
  if (identical(aesthetic, "fill") && inherits(plot, "ggchord") &&
      identical(spec$position, "right") && is.null(plot$ggchord$theme_explicit[["legend.position"]]))
    spec$position <- "left"
  # Generic legend key sizes describe categorical swatches, not the full
  # numerical ruler. Reserve enough physical distance for its closest labels.
  if (!is.null(plot)) {
    horizontal <- identical(spec$direction, "horizontal") ||
      spec$position %in% c("top", "bottom")
    style <- spec$theme$legend.text
    points <- if(inherits(style,"element_text") && is.numeric(style@size)) style@size else 8
    family <- if(inherits(style,"element_text")) style@family %||% "" else ""
    face <- if(inherits(style,"element_text")) style@face %||% "plain" else "plain"
    values <- c(0,50,80,90,95,100)
    context <- ggchord_render_context()
    dimensions <- if (inherits(style, "element_blank")) rep(0, length(values)) else
      if (!is.null(context) && dev.cur() != 1L) {
        vapply(as.character(values), function(label) {
          g <- textGrob(label, gp = gpar(fontsize = points, fontfamily = family, fontface = face))
          if (horizontal) convertWidth(grobWidth(g), "mm", TRUE) else
            convertHeight(grobHeight(g), "mm", TRUE)
        }, numeric(1))
      } else {
        points * 25.4 / 72 * if (horizontal) nchar(as.character(values)) * .6 else
          rep(1, length(values))
      }
    required <- max(((head(dimensions,-1L)+tail(dimensions,-1L))/2 + 1.2) /
      (diff(values)/100))
    field <- if (horizontal) "key.width" else "key.height"
    thickness_field <- if (horizontal) "key.height" else "key.width"
    explicit_theme <- plot$ggchord$theme_explicit
    thickness <- explicit_theme[[paste0("legend.", thickness_field)]] %||%
      explicit_theme$legend.key.size %||% ggchord_identity_bar_width()
    spec$theme <- spec$theme + do.call(theme,
      setNames(list(thickness), paste0("legend.", thickness_field)))
    explicit <- explicit_theme[[paste0("legend.", field)]] %||%
      explicit_theme$legend.key.size
    if(is.null(explicit)) {
      length_mm <- max(required, (if(horizontal) 100 else 80) * spec$size_scale)
      spec$theme <- spec$theme + do.call(theme,
        setNames(list(unit(length_mm,"mm")),paste0("legend.",field)))
    } else if(dev.cur() != 1L && is.unit(explicit) &&
        convertUnit(explicit,"mm",valueOnly=TRUE) < required) {
      warning(paste0("Identity guide labels need at least ",round(required,1),
        " mm along the bar at this font size; the explicit guide size is retained"),call.=FALSE)
    }
  }
  if (is.null(plot))
    spec$theme <- theme(legend.key.width = ggchord_identity_bar_width())
  ggchord_identity_title_gap(guide_coloursteps(available_aes=aesthetic, even.steps=FALSE,
    show.limits=TRUE, position=spec$position, direction=spec$direction,
    theme=spec$theme, order=2))
}

# Measure actual label envelopes before choosing a categorical guide grid.
# The render loop subsequently measures the combined guide/panel gtable.
ggchord_guide_columns <- function(labels, spec) {
  if (!length(labels)) return(NULL)
  context <- ggchord_render_context()
  canvas <- context$canvas_mm %||% ggchord_render_dimensions()
  if (any(!is.finite(canvas))) return(NULL)
  text_style <- spec$theme$legend.text
  size <- if (inherits(text_style,"element_text") && is.numeric(text_style@size)) text_style@size else 8 * spec$size_scale
  metrics <- if(inherits(text_style,"element_blank")) matrix(c(0,0),nrow=1,dimnames=list(NULL,c("width","height"))) else
    do.call(rbind,lapply(as.character(labels),function(label) {
      g <- textGrob(label,gp=gpar(fontsize=size,
        fontfamily=if(inherits(text_style,"element_text"))text_style@family %||% "" else "",
        fontface=if(inherits(text_style,"element_text"))text_style@face %||% "plain" else "plain"))
      c(width=convertWidth(grobWidth(g),"mm",TRUE),
        height=convertHeight(grobHeight(g),"mm",TRUE))
    }))
  dimension <- function(field,fallback) {
    value<-spec$theme[[field]]
    if(is.unit(value))as.numeric(convertUnit(value,"mm",valueOnly=TRUE))[1] else fallback
  }
  key_height<-dimension("legend.key.height",4*spec$size_scale)
  key_width<-dimension("legend.key.width",5.5*spec$size_scale)
  pitch <- max(key_height,max(metrics[,"height"])+2)
  if (spec$position %in% c("top","bottom")) {
    capacity <- max(1,floor(canvas[1]*.85/(max(metrics[,"width"])+key_width+2)))
    return(as.integer(min(length(labels),capacity)))
  }
  as.integer(max(1,ceiling(length(labels)*pitch/(spec$height_budget %||% (canvas[2]*.72)))))
}

# Share the side budget between the actual visible roles; the final gtable
# measurement remains authoritative and can report infeasible explicit styles.
ggchord_guide_height_budget <- function(plot,position) {
  canvas <- ggchord_render_dimensions()
  if(any(!is.finite(canvas)))return(NULL)
  roles <- character()
  for(layer in plot$layers) {
    if(identical(layer$show.legend,FALSE))next
    for(aes in names(ggchord_effective_mapping(layer))) {
      if (isFALSE(layer$show.legend[aes])) next
      role <- if(grepl("^track_(fill|colour)_",aes)) aes else switch(aes,
        seq_colour="seq",ribbon_fill="ribbon",ribbon_colour="ribbon",
        gene_fill="gene",feature_fill="feature",region_fill="region",NULL)
      if(is.null(role))next
      scale <- plot$scales$get_scales(aes)
      guide <- plot$guides$guides[[aes]] %||% scale$guide
      if (identical(guide, "none") || isFALSE(guide) || inherits(guide, "GuideNone")) next
      where <- if (inherits(guide, "Guide")) guide$params$position else NULL
      where <- where %||% plot$theme$legend.position %||% if (role == "ribbon") "left" else "right"
      if (identical(where, position)) roles <- c(roles, role)
    }
  }
  max(12,canvas[2]*.72/max(1,length(unique(roles)))-5)
}
