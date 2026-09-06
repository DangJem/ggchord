#' Scales for point links
#'
#' Point-link aesthetics are independent of sequence and ribbon aesthetics.
#' @param ... Arguments passed to the corresponding ggplot2 scale.
#' @param values Manual aesthetic values.
#' @return A ggplot2 scale.
#' @name scale_link
NULL

#' @rdname scale_link
#' @export
scale_link_colour_manual <- function(..., values) ggplot2::scale_colour_manual(..., values = values, aesthetics = "link_colour")
#' @rdname scale_link
#' @export
scale_link_color_manual <- scale_link_colour_manual
#' @rdname scale_link
#' @export
scale_link_alpha_manual <- function(..., values) ggplot2::scale_alpha_manual(..., values = values, aesthetics = "link_alpha")
#' @rdname scale_link
#' @export
scale_link_linewidth_manual <- function(..., values) ggplot2::scale_linewidth_manual(..., values = values, aesthetics = "link_linewidth")
#' @rdname scale_link
#' @export
scale_link_linetype_manual <- function(..., values) ggplot2::scale_linetype_manual(..., values = values, aesthetics = "link_linetype")
#' @rdname scale_link
#' @export
scale_link_colour_gradientn <- function(...) ggplot2::scale_colour_gradientn(..., aesthetics = "link_colour")
#' @rdname scale_link
#' @export
scale_link_color_gradientn <- scale_link_colour_gradientn
#' @rdname scale_link
#' @export
scale_link_alpha_continuous <- function(...) ggplot2::scale_alpha_continuous(..., aesthetics = "link_alpha")
#' @rdname scale_link
#' @export
scale_link_linewidth_continuous <- function(...) ggplot2::scale_linewidth_continuous(..., aesthetics = "link_linewidth")

ggchord_add_link_scales <- function(plot) {
  for (aesthetic in c("link_colour", "link_alpha", "link_linewidth", "link_linetype")) {
    if (plot$scales$has_scale(aesthetic)) next
    info <- ggchord_mapping_info(plot, aesthetic)
    if (is.null(info)) next
    guide <- ggchord_role_guide(plot, "link", colourbar = FALSE, order = 2)
    if (info$kind == "continuous") {
      scale <- switch(aesthetic,
        link_colour = scale_link_colour_gradientn(name = info$label, colours = c("#34457E", "#2FA96B"), guide = guide),
        link_alpha = scale_link_alpha_continuous(name = info$label, guide = guide),
        link_linewidth = scale_link_linewidth_continuous(name = info$label, guide = guide),
        ggchord_stop("link_linetype requires discrete values"))
    } else {
      levels <- unique(as.character(info$values)); levels <- levels[!is.na(levels)]
      n <- length(levels)
      values <- switch(aesthetic, link_colour = chord_default_palette(n),
        link_alpha = seq(.35, .9, length.out = n), link_linewidth = seq(.3, 1.2, length.out = n),
        link_linetype = rep(1:6, length.out = n))
      values <- stats::setNames(values, levels)
      fun <- switch(aesthetic, link_colour = scale_link_colour_manual,
        link_alpha = scale_link_alpha_manual, link_linewidth = scale_link_linewidth_manual,
        link_linetype = scale_link_linetype_manual)
      scale <- fun(name = info$label, values = values, guide = guide)
    }
    scale$ggchord_managed <- TRUE
    plot$scales$add(scale)
  }
  plot
}
