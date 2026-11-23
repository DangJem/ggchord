#' Colour sequence identifiers
#'
#' Use the ggchord sequence palette on a discrete colour aesthetic. In a
#' ggchord plot each registered reference keeps its colour across layers,
#' including layers displaying a subset of references. On a plain ggplot
#' colours follow the scale limits.
#'
#' Use `aesthetics = c("colour", "fill")` to share one scale between a
#' backbone mapped with `aes(colour = accver)` and features mapped with
#' `aes(fill = accver)`. Their glyphs can then form one composition per category.
#' @param name Legend title.
#' @param ... Additional arguments for the underlying ggplot2 discrete scale,
#'   such as `limits`, `breaks`, `labels`, `na.value`, or `guide`.
#' @param aesthetics The colour or fill aesthetics to control.
#' @return A ggplot2 scale that can also be used inside [scale_for()].
#' @export
#' @examples
#' library(ggplot2)
#' ggchord(seq_data_example) +
#'   geom_seq(aes(colour = accver)) +
#'   scale_colour_sequence(name = "Sequence")
#'
#' ggchord(seq_data_example) +
#'   geom_seq(aes(colour = accver)) +
#'   geom_gene(aes(fill = accver), data = gene_data_example) +
#'   scale_colour_sequence(name = "Sequence", aesthetics = c("colour", "fill"))
scale_colour_sequence <- function(name = waiver(), ..., aesthetics = "colour") {
  ggchord_local_console_errors()
  scale <- discrete_scale(aesthetics, palette = chord_default_palette, name = name, ...)
  scale$ggchord_palette <- "sequence"
  scale
}

#' Colour features by biological strand
#'
#' Forward and reverse strands use red and blue. Bidirectional features use
#' purple and undirected features use grey.
#' @inheritParams scale_colour_sequence
#' @inherit scale_colour_sequence return
#' @export
#' @examples
#' library(ggplot2)
#' ggchord(seq_data_example) +
#'   geom_gene(aes(fill = strand), data = gene_data_example) +
#'   scale_fill_strand(name = "Strand")
scale_fill_strand <- function(name = waiver(), ..., aesthetics = "fill") {
  ggchord_local_console_errors()
  palette <- c(process_strand_colors(NULL), "+/-" = "#795788", "." = "#767B80")
  scale_fill_manual(name = name, values = palette, aesthetics = aesthetics, ...)
}

#' Colour genomic feature types
#'
#' Known feature types use the plasmid palette. Other type names receive a
#' colour derived from their spelling, so their colours remain stable when
#' filtering data or adding types. Missing types use the scale's `na.value`.
#' @inheritParams scale_colour_sequence
#' @inherit scale_colour_sequence return
#' @export
#' @examples
#' library(ggplot2)
#' features <- transform(gene_data_example, type = "CDS")
#' ggchord(seq_data_example) +
#'   geom_gene(aes(fill = type), data = features) +
#'   scale_fill_feature_type()
scale_fill_feature_type <- function(name = waiver(), ..., aesthetics = "fill") {
  ggchord_local_console_errors()
  scale <- discrete_scale(aesthetics, palette = ggchord_feature_type_palette,
    name = name, ...)
  scale
}

ggchord_feature_type_palette <- function(self, n) {
  keys <- as.character(self$get_limits())
  keys <- keys[!is.na(keys)]
  palette <- ggchord_plasmid_feature_colours()
  extra <- setdiff(keys, names(palette))
  if (length(extra)) {
    hue <- vapply(enc2utf8(extra), function(key) {
      hash <- 0
      for (code in utf8ToInt(key)) hash <- (hash * 31 + code) %% 360
      hash
    }, numeric(1))
    palette[extra] <- grDevices::hcl(h = hue, c = 52, l = 62)
  }
  palette[keys]
}

#' Colour ribbons by percentage identity
#'
#' Map percentage identity to a stepped purple-to-yellow fill. The default
#' guide displays intervals in proportion to their widths.
#' @param name Legend title.
#' @param ... Additional arguments for [ggplot2::scale_fill_stepsn()].
#' @param limits Numeric limits of the percentage domain.
#' @param breaks Numeric boundaries of the identity intervals.
#' @param guide A fill guide. The default uses a proportional stepped colour bar.
#' @inherit scale_colour_sequence return
#' @export
#' @examples
#' library(ggplot2)
#' ggchord(seq_data_example) +
#'   geom_link_ribbon(aes(fill = pident), data = ribbon_data_example) +
#'   scale_fill_pident()
scale_fill_pident <- function(name = "Identity (%)", ..., limits = c(0, 100),
    breaks = c(0, 50, 80, 90, 95, 100), guide = waiver()) {
  ggchord_local_console_errors()
  automatic <- is_waiver(guide)
  if (automatic) guide <- ggchord_identity_guide(aesthetic = "fill")
  scale <- ggplot2::scale_fill_stepsn(name = name,
    colours = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
      "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"),
    limits = limits, breaks = breaks, guide = guide, ...)
  if (automatic) scale$ggchord_palette <- "pident"
  scale
}

# Resolve reference order and legend dimensions once the plot exists.
ggchord_prepare_domain_scale <- function(scale, plot) {
  if (is.null(scale$ggchord_palette)) return(scale)
  copy <- scale$clone()
  if (identical(copy$ggchord_palette, "sequence")) {
    layout <- plot$ggchord$ref$layout
    palette <- layout$legend_seq_colors %||% layout$seq_colors
    if (length(palette)) {
      copy$palette <- function(self, n) {
        extra <- setdiff(as.character(stats::na.omit(self$get_limits())), names(palette))
        if (!length(extra)) return(palette)
        colours <- chord_default_palette(length(palette) + length(extra))
        c(palette, setNames(tail(colours, length(extra)), extra))
      }
    }
  } else if (identical(copy$ggchord_palette, "pident")) {
    copy$guide <- ggchord_identity_guide(plot, "fill")
  }
  copy
}
