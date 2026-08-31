# Role-specific scales for ggchord. Keeping roles separate lets one plot use
# independent sequence, ribbon, gene, feature and region scales without the
# usual fill/colour collisions.

#' Sequence colour scales
#'
#' @param ... Arguments passed to [ggplot2::scale_colour_manual()].
#' @param values A set of aesthetic values to map data values to.
#' @return A ggplot2 scale.
#' @export
scale_seq_colour_manual <- function(..., values) {
  ggplot2::scale_colour_manual(..., values = values, aesthetics = "seq_colour")
}

#' @rdname scale_seq_colour_manual
#' @export
scale_seq_color_manual <- scale_seq_colour_manual

#' Ribbon fill scales
#'
#' @param ... Arguments passed to the corresponding ggplot2 scale.
#' @param colours,colors Gradient colours.
#' @param values Manual values, or positions for gradient colours where
#'   supported by the corresponding ggplot2 scale.
#' @param guide A guide object; defaults accept the `ribbon_fill` aesthetic.
#' @return A ggplot2 scale.
#' @export
scale_ribbon_fill_stepsn <- function(
    ..., colours, values = NULL, colors,
    guide = ggplot2::guide_coloursteps(available_aes = "ribbon_fill")) {
  if (missing(colours)) colours <- colors
  ggplot2::scale_fill_stepsn(
    ..., colours = colours, values = values, aesthetics = "ribbon_fill",
    guide = guide
  )
}

#' @rdname scale_ribbon_fill_stepsn
#' @export
scale_ribbon_fill_gradientn <- function(
    ..., colours, values = NULL, colors,
    guide = ggplot2::guide_colourbar(available_aes = "ribbon_fill")) {
  if (missing(colours)) colours <- colors
  ggplot2::scale_fill_gradientn(
    ..., colours = colours, values = values, aesthetics = "ribbon_fill",
    guide = guide
  )
}

#' @rdname scale_ribbon_fill_stepsn
#' @export
scale_ribbon_fill_manual <- function(..., values) {
  ggplot2::scale_fill_manual(..., values = values, aesthetics = "ribbon_fill")
}

#' @rdname scale_ribbon_fill_stepsn
#' @export
scale_ribbon_fill_identity <- function(...) {
  ggplot2::scale_fill_identity(..., aesthetics = "ribbon_fill")
}

#' Ribbon alpha, outline and linetype scales
#'
#' @param ... Arguments passed to the corresponding ggplot2 scale.
#' @param range Output alpha range for a continuous scale.
#' @param values Manual aesthetic values.
#' @return A ggplot2 scale.
#' @export
scale_ribbon_alpha_continuous <- function(..., range = NULL) {
  ggplot2::scale_alpha_continuous(
    ..., range = range, aesthetics = "ribbon_alpha"
  )
}

#' @rdname scale_ribbon_alpha_continuous
#' @export
scale_ribbon_alpha_manual <- function(..., values) {
  ggplot2::scale_alpha_manual(..., values = values, aesthetics = "ribbon_alpha")
}

#' @rdname scale_ribbon_alpha_continuous
#' @export
scale_ribbon_colour_manual <- function(..., values) {
  ggplot2::scale_colour_manual(
    ..., values = values, aesthetics = "ribbon_colour"
  )
}

#' @rdname scale_ribbon_alpha_continuous
#' @export
scale_ribbon_color_manual <- scale_ribbon_colour_manual

#' @rdname scale_ribbon_alpha_continuous
#' @export
scale_ribbon_linetype_manual <- function(..., values) {
  ggplot2::scale_linetype_manual(
    ..., values = values, aesthetics = "ribbon_linetype"
  )
}

#' Gene, feature and region fill scales
#'
#' @param ... Arguments passed to [ggplot2::scale_fill_manual()].
#' @param values A set of fill values to map data values to.
#' @return A ggplot2 scale.
#' @export
scale_gene_fill_manual <- function(..., values) {
  ggplot2::scale_fill_manual(..., values = values, aesthetics = "gene_fill")
}

#' @rdname scale_gene_fill_manual
#' @export
scale_feature_fill_manual <- function(..., values) {
  ggplot2::scale_fill_manual(..., values = values, aesthetics = "feature_fill")
}

#' Feature shape scale
#'
#' Maps feature categories to the four geometry types understood by
#' [geom_feature()]: `"arrow"`, `"block"`, `"chevron"`, and `"lollipop"`.
#' Shape values affect the actual feature geometry as well as its legend key.
#'
#' @param ... Arguments passed to [ggplot2::discrete_scale()].
#' @param values A named or unnamed character vector of feature geometry names.
#' @param name Scale and guide title, default `"Feature"`.
#' @param limits Optional category order. Named `values` use their name order
#'   by default so shape and fill guides can merge cleanly.
#' @param guide Guide specification. The default is `"none"` because
#'   [geom_feature()] combines shape silhouettes into its fill legend when
#'   both aesthetics describe the same feature categories. Supply
#'   [guide_ggchord_legend()] for a separate shape guide.
#' @return A ggplot2 discrete scale for the `feature_shape` aesthetic.
#' @export
scale_feature_shape_manual <- function(..., values, name = "Feature",
                                       limits = NULL, guide = "none") {
  allowed <- c("arrow", "block", "chevron", "lollipop")
  if (!is.character(values) || length(values) == 0L || anyNA(values) ||
      any(!values %in% allowed)) {
    ggchord_stop(
      "scale_feature_shape_manual(): values must use 'arrow', 'block', ",
      "'chevron', or 'lollipop'"
    )
  }
  if (is.null(limits) && !is.null(names(values))) limits <- names(values)
  ggplot2::discrete_scale(
    aesthetics = "feature_shape",
    palette = scales::manual_pal(values),
    ..., name = name, limits = limits, na.value = "arrow", guide = guide
  )
}

#' @rdname scale_gene_fill_manual
#' @export
scale_region_fill_manual <- function(..., values) {
  ggplot2::scale_fill_manual(..., values = values, aesthetics = "region_fill")
}

#' Genomic sequence position scale
#'
#' Controls major/minor genomic position breaks and labels independently for
#' every sequence. The scale is trained against each sequence's `[0, length]`
#' range during chord layout.
#'
#' @param name Scale name; position guides are disabled by default.
#' @param breaks,minor_breaks,labels,limits,expand,oob,transform Standard
#'   continuous-scale controls.
#' @param ... Additional arguments passed to [ggplot2::continuous_scale()].
#' @return A ggplot2 continuous scale for the `seq_position` role.
#' @export
scale_seq_position_continuous <- function(
    name = ggplot2::waiver(), breaks = ggplot2::waiver(),
    minor_breaks = ggplot2::waiver(), labels = ggplot2::waiver(),
    limits = NULL, expand = ggplot2::waiver(), oob = scales::censor,
    transform = "identity", ...) {
  ggplot2::continuous_scale(
    aesthetics = "seq_position",
    palette = function(x) x,
    name = name, breaks = breaks, minor_breaks = minor_breaks,
    labels = labels, limits = limits, expand = expand, oob = oob,
    transform = transform, guide = "none", ...
  )
}
