# geom-gene.R - gene arrow layer and gene label layer
# The gene shortcut shares feature construction; layer-owned inputs and
# parameters are resolved by the plot-owned layout calculation at build time.

#' Draw directional sequence features
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param shape Default shape. An explicitly supplied fixed shape takes priority over its mapping.
#' @param thickness Physical thickness in mm or absolute grid units. Named values select references.
#' @param stat A compatible ggplot2 Stat or its name. See [geom_feature()] for
#'   source grouping and reuse during rendering.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical controlling legend participation
#'   as in ggplot2. The default `NA` uses the layer's mappings, and `FALSE`
#'   omits the layer from guides. `TRUE` can include a fixed-style gene key in
#'   an existing discrete guide trained by another layer. A fixed-style arrow
#'   can float above the backbone when it joins one compatible legend group.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details Provide one-based closed intervals with `accver`, `start` and `end`. A strand or disjoint segments can be supplied where appropriate. The gene shortcut uses strand placement and a shoulderless arrow by default.
#' Set `head_style = "shouldered"` through `...` to use a wider arrowhead.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' gene_colours <- c(CDS = "#187B77", tRNA = "#CC7D3A")
#'
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_gene(aes(fill = type), data = genes, id = "genes") +
#'   scale_fill_manual(values = gene_colours) +
#'   coord_collinear()
geom_gene <- function(mapping = NULL, data = NULL, shape = "arrow",
    thickness = 2, stat = "identity", position = position_strand(), ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_annotation_layer(mapping, data, "interval", shape, !missing(shape), thickness,
    stat, position, list(...), na.rm, show.legend, inherit.aes, id,
    default_head_style = "flush", feature_role = "gene")
}

#' Label gene intervals
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional label table or function. With a bound anchor, the default data consists of the target entities.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param anchor_layer Use [auto_anchor()] to follow the unique compatible layer, a layer ID to select a target, or `NULL` for independent labels.
#' @param anchor Position on the entity or reference path that anchors the text.
#' @param position Position adjustment. For bound labels, leave this unset to follow the target layer.
#' @param label_layout A fixed or repelled label layout declaration.
#' @param leader How to route the optional label leader: `"auto"`, `"none"`,
#'   `"straight"` or `"elbow"`. Set fixed `segment_colour`, `segment_linewidth`,
#'   `segment_linetype` or `segment_alpha` through `...`.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details The source entities have `accver`, `start`, `end` and a label or label mapping. A bound label follows the target layer by identity. Text size, collision handling and leader placement are independent of the underlying feature geometry.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' ggchord(seqs) +
#'   geom_gene(data = genes, id = "genes") +
#'   geom_gene_label(aes(label = label), anchor_layer = "genes") +
#'   coord_collinear()
geom_gene_label <- function(mapping = NULL, data = NULL, stat = "identity",
    anchor_layer = auto_anchor(), anchor = "middle", position = NULL,
    label_layout = layout_label_fixed(orientation = "horizontal"), leader = "auto",
    ..., na.rm = FALSE,
    show.legend = FALSE, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  layer <- geom_feature_label(mapping, data, stat, anchor_layer, anchor, position, label_layout, leader,
    ..., na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes, id = id)
  layer$ggchord_params$anchor_role <- "gene"
  layer
}

#' Label primer intervals
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional label table or function. With a bound anchor, the default data consists of the target entities.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param anchor_layer Use [auto_anchor()] to follow the unique compatible layer, a layer ID to select a target, or `NULL` for independent labels.
#' @param anchor Position on the entity or reference path that anchors the text.
#' @param position Position adjustment. For bound labels, leave this unset to follow the target layer.
#' @param label_layout A fixed or repelled label layout declaration.
#' @param leader How to route the optional label leader: `"auto"`, `"none"`,
#'   `"straight"` or `"elbow"`. Set fixed `segment_colour`, `segment_linewidth`,
#'   `segment_linetype` or `segment_alpha` through `...`.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details The source entities have `accver`, `start`, `end` and a label or label mapping. A bound label follows the target layer by identity. Text size, collision handling and leader placement are independent of the underlying feature geometry.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' primers <- transform(genes, label = paste0("primer_", feature_id))
#' ggchord(seqs) +
#'   geom_primer(data = primers, id = "primers") +
#'   geom_primer_label(aes(label = label), anchor_layer = "primers") +
#'   coord_collinear()
geom_primer_label <- function(mapping = NULL, data = NULL, stat = "identity",
    anchor_layer = auto_anchor(), anchor = "head", position = NULL,
    label_layout = layout_label_fixed(), leader = "auto", ..., na.rm = FALSE,
    show.legend = FALSE, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  layer <- geom_feature_label(mapping, data, stat, anchor_layer, anchor, position, label_layout, leader,
    ..., na.rm = na.rm, show.legend = show.legend, inherit.aes = inherit.aes, id = id)
  layer$ggchord_params$anchor_role <- "primer"
  layer
}
