#' Add a circular-sequence centre label
#'
#' @param mapping,data Standard layer inputs.
#' @param name Show the sequence name. A character value overrides the name.
#' @param length Show the formatted sequence length.
#' @param separator Separator between name and length.
#' @param name_style,length_style Named lists of fixed text-style overrides
#'   for the name and length grobs. Supported fields are `colour`, `size`,
#'   `family`, `fontface`, and `lineheight`.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to `geom_text()`.
#' @return A ggplot2 layer.
#' @export
geom_seq_center_label <- function(
    mapping = NULL, data = NULL,
    name = TRUE, length = TRUE, separator = "\n",
    name_style = list(), length_style = list(),
    show.legend = FALSE, inherit.aes = FALSE, ...) {
  if (!(is.logical(name) && base::length(name) == 1L && !is.na(name)) &&
      !(is.character(name) && base::length(name) == 1L && !is.na(name))) {
    ggchord_stop("geom_seq_center_label(): name must be TRUE/FALSE or one character value")
  }
  if (!is.logical(length) || base::length(length) != 1L || is.na(length)) {
    ggchord_stop("geom_seq_center_label(): length must be TRUE or FALSE")
  }
  if (!is.character(separator) || base::length(separator) != 1L || is.na(separator)) {
    ggchord_stop("geom_seq_center_label(): separator must be one character value")
  }
  allowed_style <- c("colour", "color", "size", "family", "fontface",
    "lineheight")
  for (item in list(name_style = name_style, length_style = length_style)) {
    if (!is.list(item) || (length(item) &&
        (is.null(names(item)) || any(!names(item) %in% allowed_style)))) {
      ggchord_stop(
        "geom_seq_center_label(): style overrides must be named lists using ",
        paste(allowed_style, collapse = ", ")
      )
    }
  }
  lyr <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), label = character()),
    mapping = ggplot2::aes(
      x = x, y = y, label = label, .component = I(.component)
    ),
    stat = "identity", geom = GeomSeqCenterLabel, position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes,
    check.aes = FALSE, check.param = FALSE,
    params = c(list(na.rm = FALSE, center_name_params = name_style,
      center_length_params = length_style), list(...))
  )
  lyr$ggchord_type <- "seq_center_label"
  lyr$ggchord_theme_components <- c(
    center_name_params = "ggchord.seq.center.name",
    center_length_params = "ggchord.seq.center.length"
  )
  lyr$ggchord_params <- list(
    type = "seq_center_label", name = name,
    length = length, separator = separator
  )
  ggchord_capture_layer_input(
    lyr, data, mapping, c("accver", "length", "label")
  )
}

# Resolve the same fixed component styles used by the centre-label renderer
# before measuring the title block. A theme value only fills an unspecified
# component field; explicit name_style/length_style values stay authoritative.
ggchord_seq_center_label_styles <- function(layer, plot, seq_data) {
  defaults <- lapply(GeomSeqCenterLabel$default_aes, rlang::eval_tidy)
  mapping <- ggchord_effective_mapping(layer)
  style_fields <- c("size", "family", "fontface", "lineheight", "angle",
    "hjust", "vjust", "colour", "alpha")
  mapped <- NULL
  if (any(style_fields %in% names(mapping))) {
    # Resolve visual mappings at ggplot2's own stages, using the very input
    # rows that will be attached to these two rendered components. This also
    # respects identity aesthetics and user scales without hand-evaluating an
    # after_stat()/after_scale() expression as if it were a start expression.
    geometry <- ggchord_seq_center_label_geometry(seq_data,
      layer$ggchord_params)
    geometry <- ggchord_attach_input_columns(geometry, seq_data)
    style_plot <- ggplot2::ggplot()
    style_plot$layers <- list(reconstruct_layer(layer, geometry))
    style_plot$scales <- plot$scales$clone()
    style_plot$theme <- plot$theme
    mapped <- ggplot2::ggplot_build(style_plot)$data[[1L]]
  }
  styles <- list()
  for (component in c("name", "length")) {
    values <- utils::modifyList(defaults, layer$aes_params)
    if (is.data.frame(mapped) && ".component" %in% names(mapped)) {
      row <- which(mapped$.component == component)[1L]
      if (length(row) && is.finite(row)) {
        for (field in intersect(style_fields, names(mapped))) {
          values[[field]] <- mapped[[field]][row]
        }
      }
    }
    override <- layer$geom_params[[paste0("center_", component,
      "_params")]] %||% list()
    if (!is.null(override$color)) {
      override$colour <- override$color
      override$color <- NULL
    }
    element <- ggchord_theme_element(plot,
      paste0("ggchord.seq.center.", component))
    if (inherits(element, "element_blank")) {
      override$alpha <- 0
    } else if (inherits(element, "element_text")) {
      theme_values <- list(colour = element@colour,
        size = element@size / ggplot2::.pt, family = element@family,
        fontface = element@face, lineheight = element@lineheight)
      for (field in names(theme_values)) {
        if (is.null(override[[field]]) && is.null(mapping[[field]]) &&
            !is.null(theme_values[[field]])) {
          override[[field]] <- theme_values[[field]]
        }
      }
    }
    styles[[component]] <- utils::modifyList(values, override)
  }
  styles
}

ggchord_seq_center_label_geometry <- function(seq_data, params,
                                               units_per_inch = .30,
                                               styles = list()) {
  if (!is.data.frame(seq_data) || nrow(seq_data) != 1L) {
    ggchord_stop("geom_seq_center_label() requires exactly one sequence")
  }
  requested_name <- params$name
  name <- if (is.character(requested_name)) {
    requested_name
  } else if (isTRUE(requested_name)) {
    if ("label" %in% names(seq_data) && !is.na(seq_data$label[1L]) &&
        nzchar(as.character(seq_data$label[1L]))) {
      as.character(seq_data$label[1L])
    } else as.character(seq_data$accver[1L])
  } else character()
  length_text <- if (isTRUE(params$length)) {
    paste0(format(round(seq_data$length[1L]), big.mark = ",", scientific = FALSE), " bp")
  } else character()
  text <- paste(c(name, length_text), collapse = params$separator %||% "\n")
  labels <- c(name = name, length = length_text)
  labels <- labels[nzchar(labels)]
  geometry <- data.frame(
    x = rep(0, length(labels)), y = rep(0, length(labels)),
    label = unname(labels),
    .component = as.character(names(labels)),
    combined_label = rep(text, length(labels)),
    source_row = rep(1L, length(labels)), stringsAsFactors = FALSE
  )
  for (row in seq_len(nrow(geometry))) {
    style <- styles[[geometry$.component[row]]] %||%
      lapply(GeomSeqCenterLabel$default_aes, rlang::eval_tidy)
    for (field in intersect(names(style), c("size", "family", "fontface",
        "lineheight", "angle", "hjust", "vjust", "colour", "alpha"))) {
      geometry[[field]][row] <- style[[field]]
    }
  }
  if (nrow(geometry) == 2L) {
    measured <- geometry
    measured$text <- measured$label
    measured$text_x <- measured$x
    measured$text_y <- measured$y
    measured$text_angle <- measured$angle %||% 0
    boxes <- ggchord_text_boxes(measured,
      units_per_inch = units_per_inch, box_padding = 0)
    # Reserve an absolute reading gap and centre the complete measured block.
    # Data-unit offsets would collapse when an automatic preview shrinks the
    # physical sequence circle while retaining readable point-sized text.
    gap <- .03 * units_per_inch
    geometry$y <- c(gap / 2 - boxes$ymin[1L],
      -gap / 2 - boxes$ymax[2L])
    extent <- c(boxes$ymin + geometry$y, boxes$ymax + geometry$y)
    geometry$y <- geometry$y - mean(range(extent))
  }
  geometry
}

GeomSeqCenterLabel <- ggplot2::ggproto(
  "GeomSeqCenterLabel", ggplot2::Geom,
  required_aes = c("x", "y", "label", ".component"),
  default_aes = ggplot2::aes(
    colour = "#202020", size = 3.5, angle = 0, hjust = .5, vjust = .5,
    alpha = 1, family = "", fontface = 1, lineheight = 1.2
  ),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        center_name_params = list(),
                        center_length_params = list()) {
    grobs <- list()
    stacked <- all(c("name", "length") %in% data$.component)
    for (component in c("name", "length")) {
      part <- data[data$.component == component, , drop = FALSE]
      if (!nrow(part)) next
      params <- if (component == "name") {
        center_name_params
      } else {
        center_length_params
      }
      for (nm in names(params)) {
        target <- if (nm == "color") "colour" else nm
        if (target %in% names(part)) part[[target]] <- params[[nm]]
      }
      if (stacked) part$y <- 0
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
        part, panel_params, coord, parse = FALSE, check_overlap = FALSE,
        na.rm = na.rm
      )
    }
    if (stacked && length(grobs) == 2L) {
      # Shared scales can change a mapped font after the layout's preliminary
      # measurement. Use the final styled grobs and absolute grid units here
      # so even that later change cannot collapse the two printed lines.
      gap <- .03
      heights <- vapply(grobs, function(grob) grid::convertHeight(
        grid::grobHeight(grob), "inches", valueOnly = TRUE), numeric(1L))
      shifts <- c((heights[2L] + gap) / 2,
        -(heights[1L] + gap) / 2)
      for (i in seq_along(grobs)) {
        # grobY() is a location, not a height. Convert each location before
        # subtracting: device-native y axes may run in the opposite direction,
        # so convertHeight(grobY(...) - y) can add a whole device height.
        bounds <- vapply(c(90, 270), function(angle) grid::convertY(
          grid::grobY(grobs[[i]], angle), "inches", valueOnly = TRUE
        ), numeric(1L))
        anchor <- grid::convertY(grobs[[i]]$y,
          "inches", valueOnly = TRUE)
        offset <- mean(bounds) - anchor
        grobs[[i]]$y <- grobs[[i]]$y +
          grid::unit(shifts[i] - offset, "inches")
      }
    }
    do.call(grid::grobTree, grobs)
  }
)
