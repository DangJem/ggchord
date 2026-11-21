# ggchord.R - constructor
# Entry point of the layered ggplot2 API for chord diagrams
# v0.6.0: plot objects are self-contained. Data and parameters are stored on
# the plot object itself, and the layout is computed at build time so that
# print(), ggsave(), ggplot_build() and other ggplot2 workflows all work.

# Global variable declarations (to avoid R CMD check NOTEs)
globalVariables(c(
  ".axis_plotmath",
  ".flow_from", ".flow_tip", ".flow_to", ".fusion_row", ".fusion_x", ".fusion_y",
  ".q_n", ".s_n", ".track_group", ".track_input_index", ".track_width", ".track_x",
  ".coverage_start", ".coverage_end", "coverage", "pos", "track_value", "key",
  "x", "y", "group", "pident", "fill", "colour", "strand", "anno", "accver",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "xend", "yend", ".component",
  "label", "label_x", "label_y", "size", "angle",
  "fill_col", "alpha", "label_hjust", "label_vjust", "label_angle",
  "linetype", "zcolour", "zregionfill", "zoutline", "zlinetype",
  "outline_col", "linetype_val", "value", "source_row", "direction",
  "seq_colour", "ribbon_fill", "ribbon_alpha",
  "ribbon_colour", "ribbon_linetype", ".ggchord_ribbon_linewidth", "gene_fill", "feature_fill",
  "feature_shape", "seq_ring", "bundle_n", "bundle_weight", "density",
  "region_fill", ".bundle_n", ".bundle_weight", ".bundle_density",
  ".feature_shape_raw", ".feature_boundaries", ".feature_internal_exit", "feature_label_colour",
  "feature_label_fill", "feature_label_mode",
  "feature_label_orientation", "plotmath_label", "combined_label",
  "enzyme_label", "coordinate_label", "label_order", "enzyme_fontface",
  ".draw_as_arc", "annotation_class", "primer_name", "primer_start",
  "primer_end", "primer_show_location"
))

#' Layered sequence comparisons, circular maps and dot plots
#'
#' Compare multiple sequences with [coord_chord()], [coord_collinear()],
#' [coord_radial_synteny()] or [coord_sankey()]. Use [coord_circular()] for
#' one sequence, [coord_dotplot()] for alignment segments and [ggchord_fusion()]
#' to compose independent panels. Add geometry with `geom_*()`, and use
#' mappings/scales for data-driven styles and themes for non-data appearance.
#' The layout is computed lazily when the plot is built (e.g. via \code{print()},
#' \code{ggsave()}, or \code{ggplot_build()}).
#'
#' @param data A complete reference table with `accver` and `length`.
#'   Omitted units default to bp. An explicit reference table supplied to
#'   [geom_seq()] can also define the reference when the plot has no data.
#' @param mapping Default aesthetic mapping, constructed with [`aes()`][ggplot2::aes] from ggplot2.
#'   Evaluate inherited mappings on the data of each layer.
#' @param ... Must be empty.
#'
#' @return A ggchord object (inherits from ggplot) to which geom_* layers can be added with +
#' @seealso [ggchord_api_contract], [ggchord_data_contract],
#'   [ggchord_layout_contract], [view_ggchord()]
#' @importFrom utils head
#' @importFrom grDevices as.raster
#' @export
#'
#' @importFrom ggplot2 ggplot_build
#' @examples
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
#' ggchord(data = seqs) +
#'   geom_seq() + geom_gene(data = genes) + coord_collinear()
#'
#' ggchord() +
#'   geom_seq(data = seqs) + geom_gene(data = genes) + coord_chord()
#'
#' ggchord(data = seqs) + geom_gene(data = genes) + coord_collinear()
ggchord <- function(data = NULL, mapping = aes(), ...) {
  ggchord_check_argument_names(sys.call(), sys.function())
  if (length(list(...))) ggchord_stop("ggchord(): `...` must be empty")
  if (!ggplot2::is_mapping(mapping)) ggchord_stop("ggchord(): mapping must be aes()")
  if (!is.null(data)) data <- ggchord_reference_table(data, "ggchord()")

  p <- ggplot(data = data, mapping = ggchord_normalize_mapping(mapping)) +
    coord_chord() + theme_ggchord()
  p$coordinates$default <- TRUE
  p$ggchord <- list(
    data = list(seq_data = data),
    global = list(rotation = 0, debug = FALSE),
    ref = new.env(parent = emptyenv())
  )
  class(p) <- c("ggchord", class(p))
  p
}

#' Assign an identifier to a ggplot2 layer
#'
#' @param layer An existing ggplot2 layer.
#' @param id A nonempty identifier, unique within the plot.
#' @return A layer with the specified identity.
#' @export
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' ggchord(seqs) + geom_seq() +
#'   layer_id(geom_point(data = data.frame(x = 0, y = 0),
#'     aes(x, y), inherit.aes = FALSE), "canvas_points")
layer_id <- function(layer, id) {
  if (!inherits(layer, "Layer")) ggchord_stop("layer_id(): layer must be a ggplot2 Layer")
  ggchord_check_layer_id(id)
  copy <- ggproto(NULL, layer)
  copy$ggchord_layer_id <- id
  if (!is.null(layer$ggchord_type)) copy$ggchord_params$id <- id
  copy
}

# ====================================================================
# +.ggchord method
# ====================================================================

#' Combine a ggchord plot with ggplot2 objects
#'
#' Uses ggplot2's standard composition semantics for layers, scales, themes,
#' coordinates and annotations.
#'
#' @param e1 A ggchord object
#' @param e2 A ggplot2 component.
#' @return A ggchord object
#' @export
`+.ggchord` <- function(e1, e2) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  # Compose ordinary lists through the same per-component bookkeeping as +.
  if (identical(class(e2), "list")) return(Reduce(`+`, e2, init = e1))

  if (inherits(e2, "ggchord_scoped_scale")) {
    e1$ggchord$scoped_scales <- c(e1$ggchord$scoped_scales, list(e2))
    e1$ggchord$ref <- new.env(parent = emptyenv())
    e1$ggchord$layout <- NULL
    return(e1)
  }

  old_settings <- attr(e1$theme, "ggchord.settings")
  new_settings <- if (inherits(e2, "theme")) {
    attr(e2, "ggchord.settings")
  } else NULL
  if (inherits(e2, "theme") && is.null(new_settings) &&
      !is.null(old_settings)) {
    # Preserve chord settings across an ordinary theme(), while recording
    # explicit generic legend fields so role guides inherit those values
    # instead of mistaking ggplot2's built-in line-unit defaults for user
    # choices.
    new_settings <- old_settings
    for (field in names(new_settings$legend)) {
      value <- e2[[paste0("legend.", field)]]
      if (!is.null(value)) new_settings$legend[[field]] <- value
    }
  }
  if (inherits(e2, 'theme') && is.null(attr(e2, 'ggchord.settings')) &&
      !is.null(e2[['ggchord.axis.text.offset']])) {
    new_settings <- new_settings %||% ggchord_plot_settings(e1)
    new_settings$axis$text.offset.explicit <- TRUE
  }
  if (inherits(e2, "Scale")) {
    aesthetic <- e2$aesthetics[1L]
    existing <- if (!is.null(aesthetic)) e1$scales$get_scales(aesthetic) else NULL
    if (!is.null(e2$ggchord_track_id) && !is.null(existing$ggchord_track_id) &&
        (!identical(existing$ggchord_domain_policy, e2$ggchord_domain_policy) ||
          !identical(existing$ggchord_facet_domain, e2$ggchord_facet_domain) ||
          !identical(existing$ggchord_track_unit, e2$ggchord_track_unit)))
      ggchord_stop("Conflicting value scales for track '", e2$ggchord_track_id, "'")
    new_priority <- attr(e2, "ggchord_scale_priority") %||% 2L
    old_priority <- if (is.null(existing)) -Inf else
      attr(existing, "ggchord_scale_priority") %||% 2L
    if (!is.null(existing) && new_priority < old_priority) return(e1)
    # A user-supplied scale intentionally replaces the ggchord-managed default
    # scale of the same aesthetic; muffle ggplot2's "already present" message.
    p <- withCallingHandlers(
      NextMethod(),
      message = function(m) {
        if (grepl("already present", conditionMessage(m)) &&
            !is.null(findRestart("muffleMessage"))) {
          invokeRestart("muffleMessage")
        }
      }
    )
  } else {
    p <- NextMethod()
  }
  if (length(p$layers) > length(e1$layers))
    p <- ggchord_register_annotation_sources(p, length(e1$layers))
  if (inherits(e2, "theme")) {
    incoming <- attr(e2, "ggchord.explicit_theme") %||% as.list(e2)
    p$ggchord$theme_explicit <- ggchord_merge_theme_constraints(
      e1$ggchord$theme_explicit %||% list(), incoming, isTRUE(attr(e2, "complete")))
    if (is.null(attr(e2, "ggchord.settings"))) {
      p$ggchord$native_theme <- if (is.null(e1$ggchord$native_theme)) e2 else e1$ggchord$native_theme + e2
    } else p$ggchord$native_theme <- NULL
  }
  # ggplot2 copies plot metadata, but environments remain shared references.
  # A new composition owns a new lazy cache; invalidating e1's environment
  # would clear the base plot and let sibling plots reuse each other's layout.
  p$ggchord$ref <- new.env(parent = emptyenv())
  p$ggchord$layout <- NULL
  attr(p$theme, "ggchord.settings") <- new_settings %||% old_settings
  class(p) <- unique(c("ggchord", class(p)))
  p
}
