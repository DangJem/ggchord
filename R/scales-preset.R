#' Apply the traditional ggchord data palettes
#'
#' Choose palettes for mapped sequence, strand, feature type, alignment identity
#' and track data. Explicit ordinary or scoped scales take precedence. Adding
#' the preset a second time replaces the previous declaration.
#'
#' @param roles `"auto"` or a vector of roles selected from `"sequence"`,
#'   `"strand"`, `"feature_type"`, `"pident"`, `"link_pident"` and `"signal"`.
#' @return A palette declaration to add to a ggchord plot.
#' @export
#' @examples
#' library(ggplot2)
#' ggchord(seq_data_example) +
#'   geom_seq(aes(colour = accver)) + scales_ggchord()
scales_ggchord <- function(roles = "auto") {
  ggchord_local_console_errors()
  ggchord_id_vector(roles, "roles")
  known <- c("sequence", "strand", "feature_type", "pident", "link_pident", "signal")
  if (!identical(roles, "auto") && any(!roles %in% known))
    ggchord_stop("scales_ggchord(): roles must be 'auto' or chosen from ", paste(known, collapse = ", "))
  structure(list(roles = roles), class = "ggchord_palette_preset")
}

#' Sequence identifier palette
#' @param ids Sequence identifiers in reference order.
#' @return A named colour vector for distinct identifiers.
#' @export
pal_sequence <- function(ids) {
  ids <- unique(as.character(ids[!is.na(ids)]))
  setNames(pal_ggchord(length(ids)), ids)
}

#' Strand palette
#' @return A named colour vector for forward, reverse, bidirectional and unknown strands.
#' @export
pal_strand <- function() {
  c(process_strand_colors(NULL), "+/-" = "#795788", "." = "#767B80")
}

#' Traditional continuous palettes
#' @return Colour vectors for ribbons, point links and quantitative signals.
#' @export
pal_pident <- function() {
  c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
    "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725")
}

#' @rdname pal_pident
#' @export
pal_link_pident <- function() c("#34457E", "#2FA96B")

#' @rdname pal_pident
#' @export
pal_signal <- function() c("#31507D", "#65B4A8", "#F0D56A")

#' Apply a traditional ggchord palette with a native scale
#' @param ... Arguments for the corresponding ggplot2 scale.
#' @param aesthetics Aesthetics trained by the discrete scale.
#' @return A ggplot2 scale.
#' @export
scale_colour_ggchord <- function(..., aesthetics = "colour") {
  discrete_scale(aesthetics, palette = pal_ggchord, ...)
}

#' @rdname scale_colour_ggchord
#' @export
scale_fill_ggchord <- function(..., aesthetics = "fill") {
  scale_colour_ggchord(..., aesthetics = aesthetics)
}

#' @rdname scale_colour_ggchord
#' @export
scale_fill_plasmid <- function(...) scale_fill_feature_type(...)

#' @rdname scale_colour_ggchord
#' @export
scale_colour_link_pident <- function(...) {
  defaults <- list(name = "Identity (%)", colours = pal_link_pident(),
    limits = c(0, 100), breaks = c(0, 50, 80, 90, 95, 100))
  do.call(ggplot2::scale_colour_gradientn, utils::modifyList(defaults, list(...), keep.null = TRUE))
}

#' @rdname scale_colour_ggchord
#' @export
scale_fill_signal <- function(...) ggplot2::scale_fill_gradientn(colours = pal_signal(), ...)

#' @rdname scale_colour_ggchord
#' @export
scale_colour_signal <- function(...) ggplot2::scale_colour_gradientn(colours = pal_signal(), ...)
