#' Add a circular-sequence centre label
#'
#' @param mapping,data Standard layer inputs.
#' @param name Show the sequence name. A character value overrides the name.
#' @param length Show the formatted sequence length.
#' @param separator Separator between name and length.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to `geom_text()`.
#' @return A ggplot2 layer.
#' @export
geom_seq_center_label <- function(
    mapping = NULL, data = NULL,
    name = TRUE, length = TRUE, separator = "\n",
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
  lyr <- ggplot2::geom_text(
    data = data.frame(x = numeric(), y = numeric(), label = character()),
    mapping = ggplot2::aes(x = x, y = y, label = label),
    position = "identity", show.legend = show.legend,
    inherit.aes = inherit.aes, ...
  )
  lyr$ggchord_type <- "seq_center_label"
  lyr$ggchord_theme_element <- "ggchord.seq.center.label"
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
  data.frame(x = 0, y = 0, label = text, source_row = 1L,
             stringsAsFactors = FALSE)
}
