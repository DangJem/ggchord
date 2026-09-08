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

ggchord_seq_center_label_geometry <- function(seq_data, params) {
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
  offsets <- if (length(labels) == 2L) c(.026, -.026) else 0
  data.frame(
    x = 0, y = offsets, label = unname(labels),
    .component = names(labels), combined_label = text,
    source_row = 1L, stringsAsFactors = FALSE
  )
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
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
        part, panel_params, coord, parse = FALSE, check_overlap = FALSE,
        na.rm = na.rm
      )
    }
    do.call(grid::grobTree, grobs)
  }
)
