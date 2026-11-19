#' Plotting, data and output contracts
#'
#' Use this guide to choose where a plotting parameter belongs, how to specify
#' source coordinates and where to retrieve layout measurements.
#'
#' @section Plotting responsibilities:
#' [ggchord()] supplies reference sequences and optional default input tables.
#' Each layer retains its own data, mapping, parameters and identity. `aes()`
#' selects data roles, scales control mapping rules/guides, geoms draw entities,
#' Positions place them locally, coords select the canvas mode, and themes
#' control non-data appearance. Explicit fixed styles take precedence over
#' matching mappings; omit a fixed argument to retain its mapping. User scales
#' take precedence over automatically inferred defaults.
#'
#' `geom_feature(preset = "gene")` draws strictly directed genes. Use the
#' generic [geom_feature()] constructor for unknown-strand features. The plasmid
#' preset retains its compact placement and styling defaults. Fixed feature shape is `shape`; mapped shape is `feature_shape`.
#' Gene/feature widths use `width`, label controls use `label_*`, and point/area
#' links use `link_gap` and `link_ctrl_point`. These names do not merge the
#' distinct geometries or their defaults.
#'
#' Chord, collinear, radial synteny, Sankey, circular and Dotplot are six
#' coordinate modes. Fusion composes independently registered plots and is not
#' a coordinate mode. Dotplot requires its dedicated alignment segment layer
#' and is not a registered Fusion endpoint. Circular
#' maps reject ribbon and point-link layers; collinear maps reject
#' restriction/primer/centre-label layers. `geom_seq_scale()` is a collinear
#' sequence-unit reference bar; [scale_seq_position_continuous()] configures genomic axes.
#' A matching name never implies that every geom supports every coord.
#'
#' @section Names and migration errors:
#' The `colour` spelling is preferred in documentation; existing `color`
#' spelling aliases remain supported. Selected ggplot2/grid reexports retain
#' the original upstream object and follow that dependency's public contract.
#' Their internals and future upstream signatures are not frozen by ggchord.
#'
#' Retired arguments remain errors with replacement guidance; they do not
#' execute compatibility code. Fixed `feature_shape`, old width,
#' label and ribbon endpoint names were removed in favour of the names above.
#' `gene_offset`/`feature_offset` were also removed; use
#' `position = position_strand(offset = ...)`. Explicit strand-named offsets
#' are signed, so negate the former minus-strand magnitude when migrating it.
#'
#' Automatic gene/feature labels own their rotation, offsets, orientation and
#' overlap policy. Manual `label_rotation`, `label_radial_offset`,
#' `label_circum_offset`, `label_circum_limit`, `label_orientation` and
#' `label_overlap` supplied through their `...` now fail rather than being
#' ignored. Use [geom_gene_label()] or [geom_feature_label()] for manual
#' placement, or choose automatic geometry with `label_layout`.
#' The old `"aligned"` label request becomes `label_layout = "auto"`;
#' restriction `leader = "trunk"` becomes `leader = "radial"`.
#'
#' The constructor rejects unknown arguments. ggchord geoms check
#' extra names against the selected layer, geom and stat before construction.
#' For example, `geom_seq(linewditth = 2)` reports the spelling error and
#' suggests `linewidth`. Scale, guide and theme arguments follow the contract
#' of the function that receives them. Geometric roles
#' must be evaluated before layout, whereas supported visual mappings may use
#' `after_stat()`/`after_scale()`. Unknown legacy label prefixes never imply
#' that a similarly named replacement parameter exists.
#'
#' @section Data, units and results:
#' Input `accver` and query/subject identifiers are preserved without adding or
#' removing versions. Biological strand, display orientation and link direction
#' remain distinct. Coordinates are one-based inclusive in the declared native
#' bp or aa units after explicit
#' external-format adaptation. Missing values, defaults and `NULL` are defined
#' per responsibility, not by one package-wide coercion rule.
#' Standalone cleaning defaults to sorting reversed endpoints; preparation
#' retains direction by default. Readers retain `file`/`files` with shared
#' path expansion. See [ggchord_data_contract] for source rows, optional
#' dependencies, transformations, aggregation estimates and search limits.
#'
#' Native layout distances, sequence units, angles, numerical track domains and physical
#' mm are separate quantities. Track IDs share numerical domains and consistent
#' band declarations; their physical placement is not a numeric value scale.
#' Existing documented return members remain supported. Prefer
#' [export_ggchord_layout()] for layer-aware tables; [get_chord_layout()] remains
#' available for the complete layout. Preview paths and target-size diagnostic
#' reports retain separate purposes. See [ggchord_layout_contract] for stable
#' fields, provenance, measurement status and plot/render cache ownership.
#'
#' @name ggchord_api_contract
#' @md
#' @seealso [ggchord_data_contract], [ggchord_layout_contract]
#' @examples
#' library(dplyr)
#' data(plasmid_example_pUC19)
#' seqs <- plasmid_example_pUC19 |>
#'   select(accver, length)
#' feature_source <- find_common_features(plasmid_example_pUC19)
#' features <- feature_source |>
#'   select(accver, start, end, strand, type, anno)
#' # ASCII display copy for the default PDF help-example device.
#' features <- features |>
#'   mutate(anno = gsub("\u03b1", "alpha", anno, fixed = TRUE))
#' p <- ggchord(data = seqs) +
#'   geom_seq() +
#'   geom_feature(data = features) +
#'   coord_circular()
#' export_ggchord_layout(p, include = "feature")$metadata$coordinate
NULL
