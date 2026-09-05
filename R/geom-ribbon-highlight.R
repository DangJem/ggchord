# geom-ribbon-highlight.R - ribbon highlight layer (v0.9.0)

#' Highlight selected alignment ribbons
#'
#' Draws a second polygon layer on top of selected ribbons so they can be
#' emphasized without changing the underlying Identity(%) legend. Selection is
#' done with safe, explicit filters (row numbers, query/subject IDs, pident and
#' length ranges) or a predicate function.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param ribbon_ids Optional integer vector of original ribbon row numbers to
#'   highlight.
#' @param qaccver Optional character vector; only ribbons whose query ID is in
#'   this set are highlighted.
#' @param saccver Optional character vector; only ribbons whose subject ID is in
#'   this set are highlighted.
#' @param min_pident Optional numeric. Minimum percent identity.
#' @param max_pident Optional numeric. Maximum percent identity.
#' @param min_length Optional numeric. Minimum alignment length.
#' @param max_length Optional numeric. Maximum alignment length.
#' @param predicate Optional function taking the ribbon data.frame and returning
#'   a logical vector with one element per row. Evaluated safely (no string
#'   parsing).
#' @param fill,alpha,colour,linewidth,linetype Standard fixed polygon styles.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A ggplot2 layer.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq() + geom_link_ribbon() + geom_ribbon_highlight(ribbon_ids = 1)
#' p
geom_ribbon_highlight <- function(mapping = NULL, data = NULL,
                                  ribbon_ids = NULL,
                                  qaccver = NULL,
                                  saccver = NULL,
                                  min_pident = NULL,
                                  max_pident = NULL,
                                  min_length = NULL,
                                  max_length = NULL,
                                  predicate = NULL,
                                  fill = "#C51B7D",
                                  alpha = 0.75,
                                  colour = NA,
                                  linewidth = 0.3,
                                  linetype = 1,
                                  position = "identity",
                                  show.legend = FALSE,
                                  inherit.aes = FALSE,
                                  ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  alias <- ggchord_colour_alias(
    colour, dots, "geom_ribbon_highlight()", sys.call()
  )
  colour <- alias$colour
  dots <- alias$dots
  ggchord_reject_retired(dots, "geom_ribbon_highlight()", c(
    highlight_color = "fill",
    highlight_alpha = "alpha",
    highlight_outline_color = "colour",
    highlight_outline_width = "linewidth",
    show_legend = "show.legend"
  ))

  if (!is.null(predicate) && !is.function(predicate)) {
    ggchord_stop("predicate must be a function taking ribbon_data and returning a logical vector")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    ggchord_stop("alpha must be in [0, 1]")
  }
  validate_optional_numeric <- function(x, name, integer = FALSE) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.numeric(x) || length(x) == 0 || any(!is.finite(x)) ||
        (integer && any(x != as.integer(x)))) {
      suffix <- if (integer) "finite integer values" else "finite numeric values"
      ggchord_stop("geom_ribbon_highlight(): ", name, " must contain ", suffix)
    }
  }
  validate_optional_numeric(ribbon_ids, "ribbon_ids", integer = TRUE)
  for (nm in c("min_pident", "max_pident", "min_length", "max_length")) {
    value <- get(nm)
    validate_optional_numeric(value, nm)
    if (!is.null(value) && length(value) != 1) {
      ggchord_stop("geom_ribbon_highlight(): ", nm, " must be a single value")
    }
  }
  if (!is.null(ribbon_ids) && any(ribbon_ids < 1)) {
    ggchord_stop("geom_ribbon_highlight(): ribbon_ids must be positive row numbers")
  }
  if (!is.null(min_pident) && !is.null(max_pident) && min_pident > max_pident) {
    ggchord_stop("geom_ribbon_highlight(): min_pident cannot exceed max_pident")
  }
  if (!is.null(min_length) && !is.null(max_length) && min_length > max_length) {
    ggchord_stop("geom_ribbon_highlight(): min_length cannot exceed max_length")
  }
  if (!is.null(qaccver) && (!is.character(qaccver) || anyNA(qaccver))) {
    ggchord_stop("geom_ribbon_highlight(): qaccver must be a character vector without NA")
  }
  if (!is.null(saccver) && (!is.character(saccver) || anyNA(saccver))) {
    ggchord_stop("geom_ribbon_highlight(): saccver must be a character vector without NA")
  }

  empty_polys <- data.frame(
    x = numeric(0), y = numeric(0), group = integer(0),
    stringsAsFactors = FALSE
  )
  lyr <- ggplot2::layer(
    data        = empty_polys,
    mapping     = ggplot2::aes(x = x, y = y, group = group),
    stat        = "identity",
    geom        = ggplot2::GeomPolygon,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    check.aes   = FALSE,
    check.param = FALSE,
    params      = c(dots, list(
      fill = fill, alpha = alpha, colour = colour,
      linewidth = linewidth, linetype = linetype
    ))
  )
  lyr$ggchord_type <- "ribbon_highlight"
  lyr$ggchord_params <- list(
    type                     = "ribbon_highlight",
    ribbon_ids               = ribbon_ids,
    qaccver                  = qaccver,
    saccver                  = saccver,
    min_pident               = min_pident,
    max_pident               = max_pident,
    min_length               = min_length,
    max_length               = max_length,
    predicate                = predicate,
    highlight_color          = fill,
    highlight_alpha          = alpha,
    highlight_outline_color  = colour,
    highlight_outline_width  = linewidth
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping,
    c("qaccver", "saccver", "length", "pident", "qstart", "qend",
      "sstart", "send")
  )
  lyr
}
