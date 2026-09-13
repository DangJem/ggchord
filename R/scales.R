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

#' Explicit sequence ring scale
#'
#' Maps user-defined ring identifiers to positive sequence radii. Ring
#' membership is supplied with `geom_seq(aes(seq_ring = ...))`; ggchord never
#' guesses the number, order, or spacing of rings.
#'
#' @param ... Arguments passed to [ggplot2::discrete_scale()].
#' @param values A named or unnamed numeric vector of positive radii.
#' @param name Scale title. The guide is hidden by default because the scale
#'   controls geometry rather than a visible legend aesthetic.
#' @param limits Optional ring order.
#' @param guide Guide specification, default `"none"`.
#' @return A ggplot2 discrete scale for the `seq_ring` role.
#' @export
scale_seq_ring_manual <- function(..., values, name = "Ring",
                                  limits = NULL, guide = "none") {
  if (!is.numeric(values) || length(values) == 0L || anyNA(values) ||
      any(!is.finite(values)) || any(values <= 0)) {
    ggchord_stop(
      "scale_seq_ring_manual(): values must be finite positive radii"
    )
  }
  if (is.null(limits) && !is.null(names(values))) limits <- names(values)
  ggplot2::discrete_scale(
    aesthetics = "seq_ring", palette = scales::manual_pal(values),
    ..., name = name, limits = limits, na.value = NA_real_, guide = guide
  )
}

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
  has_colours <- !missing(colours)
  has_colors <- !missing(colors)
  if (has_colours && has_colors)
    ggchord_stop("scale_ribbon_fill_stepsn(): use only one of colours and colors")
  if (!has_colours && !has_colors)
    ggchord_stop("scale_ribbon_fill_stepsn(): colours (or colors) is required")
  if (!has_colours) colours <- colors
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
  has_colours <- !missing(colours)
  has_colors <- !missing(colors)
  if (has_colours && has_colors)
    ggchord_stop("scale_ribbon_fill_gradientn(): use only one of colours and colors")
  if (!has_colours && !has_colors)
    ggchord_stop("scale_ribbon_fill_gradientn(): colours (or colors) is required")
  if (!has_colours) colours <- colors
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
#' @param limits Optional display order. Named feature values use their name
#'   order by default.
#' @return A ggplot2 scale.
#' @export
scale_gene_fill_manual <- function(..., values) {
  ggplot2::scale_fill_manual(..., values = values, aesthetics = "gene_fill")
}

#' @rdname scale_gene_fill_manual
#' @export
scale_feature_fill_manual <- function(..., values, limits = NULL) {
  if (is.null(limits) && !is.null(names(values))) limits <- names(values)
  out <- ggplot2::scale_fill_manual(..., values = values, limits = limits,
    aesthetics = "feature_fill")
  attr(out, "ggchord_scale_priority") <- 2L
  out
}

#' Feature shape scale
#'
#' Maps feature categories to the geometry types understood by
#' [geom_feature()], including full/compact/promoter/primer arrows, markers,
#' blocks, chevrons, and lollipops.
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
  allowed <- c(
    "arrow", "compact_arrow", "promoter_arrow", "primer_arrow", "marker",
    "block", "chevron", "lollipop"
  )
  if (!is.character(values) || length(values) == 0L || anyNA(values) ||
      any(!values %in% allowed)) {
    ggchord_stop(
      "scale_feature_shape_manual(): values contain an unknown geometry"
    )
  }
  if (is.null(limits) && !is.null(names(values))) limits <- names(values)
  out <- ggplot2::discrete_scale(
    aesthetics = "feature_shape",
    palette = scales::manual_pal(values),
    ..., name = name, limits = limits, na.value = "arrow", guide = guide
  )
  attr(out, "ggchord_scale_priority") <- 2L
  out
}

ggchord_plasmid_feature_colours <- function() {
  c(
    CDS="#993366", gene="#993366", resistance_gene="#CCFFCC",
    selection_marker="#CCFFCC", reporter="#05FD14",
    peptide="#CC99B2", tag="#CC99B2",
    promoter="#FFFFFF", rep_origin="#FFFF00", replication_origin="#FFFF00",
    ori="#FFFF00", origin="#FFFF00", primer_bind="#A020F0",
    primer="#A020F0", terminator="#FFFFFF", enhancer="#FFFFFF",
    protein_bind="#31849B", binding_site="#31849B", operator="#31849B",
    RBS="#A6ACB3", polyA_signal="#A6ACB3", poly_a_signal="#A6ACB3",
    regulatory="#A6ACB3", LTR="#FFE4C4", repeat_region="#FFE4C4",
    misc_RNA="#00CCFF", MCS="#99CCFF", misc_feature="#A6ACB3"
  )
}

#' Lighten a feature colour for an external callout
#' @noRd
ggchord_feature_callout_fill <- function(fill, amount = .78) {
  vapply(as.character(fill), function(value) {
    rgb <- tryCatch(
      grDevices::col2rgb(value, alpha = TRUE)[, 1L] / 255,
      error = function(e) c(184, 189, 195, 255) / 255
    )
    mixed <- rgb[1:3] + (1 - rgb[1:3]) * amount
    grDevices::rgb(mixed[1L], mixed[2L], mixed[3L], alpha = rgb[4L])
  }, character(1L), USE.NAMES = FALSE)
}

ggchord_contrast_colour <- function(fill) {
  vapply(as.character(fill), function(value) {
    rgb <- tryCatch(grDevices::col2rgb(value)[, 1L] / 255,
      error = function(e) c(1, 1, 1))
    linear <- ifelse(rgb <= .04045, rgb / 12.92,
      ((rgb + .055) / 1.055)^2.4)
    if (sum(linear * c(.2126, .7152, .0722)) < .36) "#FFFFFF" else "#202020"
  }, character(1), USE.NAMES = FALSE)
}

#' Plasmid-map feature presets
#' @param ... Additional scale arguments.
#' @param limits Optional feature-type order.
#' @param guide Guide used by the shape preset.
#' @return A feature fill or shape scale.
#' @export
scale_feature_fill_plasmid <- function(..., limits = NULL) {
  values <- ggchord_plasmid_feature_colours()
  if (is.null(limits)) limits <- names(values)
  out <- ggplot2::scale_fill_manual(
    ..., values=values, limits=limits, aesthetics="feature_fill", na.value="#B8BDC3"
  )
  attr(out,"ggchord_scale_priority") <- 1L
  out
}

#' @rdname scale_feature_fill_plasmid
#' @export
scale_feature_shape_plasmid <- function(..., limits = NULL, guide = "none") {
  values <- c(
    CDS="arrow", gene="arrow", resistance_gene="arrow",
    selection_marker="arrow", reporter="arrow", peptide="compact_arrow",
    tag="compact_arrow",
    rep_origin="arrow", replication_origin="arrow", ori="arrow",
    promoter="promoter_arrow", primer_bind="primer_arrow",
    primer="primer_arrow", protein_bind="block", binding_site="block",
    operator="block", terminator="block", enhancer="block",
    regulatory="compact_arrow", RBS="block", polyA_signal="block",
    poly_a_signal="block", LTR="block", repeat_region="block",
    misc_RNA="block", MCS="block", misc_feature="block"
  )
  if (is.null(limits)) limits <- names(values)
  out <- ggplot2::discrete_scale(
    aesthetics="feature_shape",palette=scales::manual_pal(values),
    ...,limits=limits,na.value="block",guide=guide
  )
  attr(out,"ggchord_scale_priority") <- 1L
  out
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
