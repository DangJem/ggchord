# ggchord.R - constructor
# Entry point of the layered ggplot2 API for chord diagrams
# v0.6.0: plot objects are self-contained. Data and parameters are stored on
# the plot object itself, and the layout is computed at build time so that
# print(), ggsave(), ggplot_build() and other ggplot2 workflows all work.

# Global variable declarations (to avoid R CMD check NOTEs)
globalVariables(c(
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
#' @param data A reference data.frame/tibble with accver and length, or a `GRanges`,
#'   `Seqinfo` or `GInteractions` object with known full sequence lengths,
#'   or a named nucleotide or protein Biostrings `XStringSet`. Unnamed sequence
#'   objects must first use [as_ggchord_seq()] with explicit `ids`.
#'   Optional S4 inputs use [as_ggchord_seq()]. The `sequence_unit` column is
#'   `"bp"` (default when omitted) or `"aa"`. All rows use the same unit.
#'   Protein lengths come from complete proteins. Native aa plots accept
#'   links and generic features; gene, primer and restriction layers require bp.
#'   Arbitrary character widths are not biological sequence lengths.
#'   A result from [prepare_ggchord_data()] supplies all constructor tables.
#'   A single eligible table supplies the layer default. When several qualify,
#'   select one explicitly, e.g. `geom_feature(data = prepared$features$cds)`.
#' @param mapping Default aesthetic mapping, made with [ggplot2::aes()].
#'   Each inheriting layer uses its own source table for domain roles and
#'   starting visual mappings. With `aes(accver = ref, fill = accver)`,
#'   `ref` selects the reference and the source `accver` supplies the fill.
#'   Generated fields such as `.component` are available after projection.
#'   [ggplot2::after_stat()] and [ggplot2::after_scale()] use their respective
#'   ggplot2 stages. A layer mapping overrides the same plot aesthetic.
#'   Native ggplot2 layers inherit the reference table as their plot data.
#' @param debug Logical. Whether to output debug information, default FALSE
#' @param validate Character, default \code{"warn"}. How to run the structured
#'   input-data validation (see \code{\link{validate_ggchord_data}}):
#'   \code{"warn"} emits a single summary warning when the data has problems
#'   and caches the full report on the plot object
#'   (\code{p$ggchord$validation}); \code{"error"} stops on severe problems;
#'   \code{"none"} skips the diagnostic validation (the cheap structural
#'   checks that prevent crashes are still performed).
#' @param ... Must be empty. Add [labs()], [coord_chord()] or [theme()] to
#'   the plot to set labels, select coordinates or style its appearance.
#'
#' @return A ggchord object (inherits from ggplot) to which geom_* layers can be added with +
#' @seealso [ggchord_api_contract], [ggchord_data_contract],
#'   [ggchord_layout_contract], [view_ggchord()]
#' @importFrom utils head
#' @importFrom grDevices as.raster
#' @export
#'
#' @examples
#' library(dplyr)
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' data(gene_data_example)
#'
#' seqs <- seq_data_example |>
#'   select(accver, length)
#' links <- ribbon_data_example |>
#'   select(qaccver, saccver, qstart, qend, sstart, send)
#' genes <- gene_data_example |>
#'   transmute(accver, start, end, strand = source_strand, type = "CDS")
#' p <- ggchord(seqs) +
#'   geom_seq() +
#'   geom_link_ribbon(data = links) +
#'   geom_feature(data = genes)
#' print(p)
#' @importFrom ggplot2 ggplot_build
ggchord <- function(data, mapping = NULL, ..., debug = FALSE,
                    validate = c("warn", "error", "none")) {
  ggchord_check_argument_names(sys.call(), sys.function())
  if (missing(data)) ggchord_stop("ggchord() needs data = a sequence reference or prepared tables")
  if (!is.null(mapping) && !ggplot2::is_mapping(mapping)) ggchord_stop(
    "ggchord(): mapping must be aes(). Put extra tables in layer data or prepare_ggchord_data().")
  mapping <- ggchord_normalize_mapping(mapping)
  prepared <- if (inherits(data, "ggchord_prepared")) data else NULL
  seq_data <- if (is.null(prepared)) data else prepared$seq_data
  ribbon_data <- prepared$ribbon_data
  if (!is.null(prepared[["gene_data"]])) ggchord_stop(
    "ggchord(): this saved prepared object contains gene_data. Rebuild it with ",
    "feature_data for classified annotations, or pass its gene_data explicitly ",
    "to geom_feature(data = ..., preset = 'gene').")
  feature_data <- prepared$feature_data
  if (inherits(seq_data, "GRanges") || inherits(seq_data, "Seqinfo") ||
      inherits(seq_data, "GInteractions") || inherits(seq_data, "XStringSet") ||
      inherits(seq_data, "XString")) seq_data <- as_ggchord_seq(seq_data)
  if (inherits(ribbon_data, "GInteractions")) ribbon_data <- as_ggchord_link(ribbon_data)
  seq_data <- ggchord_normalize_accver(seq_data)
  feature_data <- ggchord_normalize_feature_input(feature_data)
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  if ("gene_data" %in% names(dots)) ggchord_stop(
    "ggchord(): use geom_feature(data = genes, preset = 'gene') or prepare_ggchord_data(feature_data = ...) with classified annotations.")
  retired <- intersect(names(dots), c("seq_data", "feature_data", "ribbon_data"))
  if (length(retired)) ggchord_stop("ggchord(): move ", paste(retired, collapse = ", "),
    " into prepare_ggchord_data() or an explicit layer data argument.")
  ggchord_reject_retired(dots, "ggchord()", c(
    title = "labs(title = ...)",
    rotation = "coord_chord(rotation = ...)",
    panel_margin = "theme(plot.margin = margin(...))",
    show_legend = "theme(legend.position = 'none')"
  ))
  if (length(dots)) {
    ggchord_stop("ggchord(): unused argument(s): ",
                 paste(names(dots), collapse = ", "))
  }

  validate <- match.arg(validate)
  # ====================================================================
  # 1. Validate data
  # ====================================================================
  if (!is.logical(debug) || length(debug) != 1 || is.na(debug)) {
    ggchord_stop("debug must be TRUE or FALSE")
  }
  # One structured validation pass supplies both the always-on crash-prevention
  # checks and the optional diagnostics. `validate = "none"` skips the costly
  # coordinate, duplicate and self-link checks but still rejects unsafe input.
  validation_result <- validate_ggchord_data(
    seq_data, ribbon_data, feature_data, strict = FALSE,
    check_coordinates = validate != "none",
    check_duplicates = validate != "none",
    check_self_links = validate != "none"
  )
  structural_errors <- ggchord_structural_validation_errors(validation_result)
  if (nrow(structural_errors) > 0) {
    ggchord_stop("ggchord(): ", structural_errors$message[1])
  }

  if (!is.null(ribbon_data)) {
    if (nrow(ribbon_data) == 0) warning("No valid alignment data in ribbon_data")
    if (debug) cat("Number of alignment data rows: ", nrow(ribbon_data), "\n")
  }
  if (debug && !is.null(feature_data)) cat("Number of feature annotation rows: ", nrow(feature_data), "\n")

  validation <- if (validate == "none") NULL else validation_result
  if (!is.null(validation)) {
    if (validate == "error" && !validation$valid) {
      ggchord_stop(sprintf(
        "ggchord(): input data failed validation (%d severe error(s); first: %s). Run validate_ggchord_data(..., strict = FALSE) for the full report.",
        nrow(validation$errors), validation$errors$message[1]),
        call. = FALSE)
    }
    if (validate == "warn") {
      n_err <- nrow(validation$errors)
      n_warn <- nrow(validation$warnings)
      if (n_err > 0) {
        warning(sprintf(
          "ggchord(): input data has %d severe validation error(s) (e.g. \"%s\"). The plot may be misleading; run validate_ggchord_data(..., strict = FALSE) for details.",
          n_err, validation$errors$message[1]), call. = FALSE)
      } else if (n_warn > 0) {
        warning(sprintf(
          "ggchord(): input data has %d validation warning(s) (e.g. \"%s\"). Run validate_ggchord_data(...) for details.",
          n_warn, validation$warnings$message[1]), call. = FALSE)
      }
    }
  }

  # ====================================================================
  # 2. Build the base ggplot object and store data + global parameters
  #    on the plot itself so the object is fully self-contained.
  # ====================================================================
  p <- ggplot2::ggplot(data = seq_data, mapping = mapping %||% ggplot2::aes()) +
    coord_chord(rotation = 0) +
    theme_ggchord()

  p$ggchord <- list(
    data   = list(seq_data = seq_data, ribbon_data = ribbon_data,
                  feature_data = feature_data),
    collections = list(features = prepared$features, links = prepared$links),
    sources = list(sequence = seq_data, relationship = ribbon_data,
      feature = feature_data),
    global = list(rotation = 0, debug = debug, validate = validate),
    validation = validation,
    # Shared reference environment for plot-owned layout caching.
    ref    = new.env(),
    layout = NULL
  )
  for (kind in c("features", "links")) {
    named <- p$ggchord$collections[[kind]]
    primary <- if (kind == "features") "feature_data" else "ribbon_data"
    present <- primary[!vapply(p$ggchord$data[primary], is.null, logical(1))]
    p$ggchord$candidates[[kind]] <- c(present, if (length(named)) paste0(kind, "$", names(named)))
    if (!length(present) && length(named) == 1L) {
      p$ggchord$data[[primary[1L]]] <- named[[1L]]
      p$ggchord$sources[[if (kind == "features") "feature" else "relationship"]] <- named[[1L]]
    }
  }
  class(p) <- c("ggchord", class(p))

  p
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
  } else if (inherits(e2, "Coord")) {
    p <- withCallingHandlers(
      NextMethod(),
      message = function(m) {
        if (grepl("Coordinate system already present", conditionMessage(m)) &&
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
  ggchord_validate_coord_contract(p)
  # ggplot2 copies plot metadata, but environments remain shared references.
  # A new composition owns a new lazy cache; invalidating e1's environment
  # would clear the base plot and let sibling plots reuse each other's layout.
  p$ggchord$ref <- new.env(parent = emptyenv())
  p$ggchord$layout <- NULL
  attr(p$theme, "ggchord.settings") <- new_settings %||% old_settings
  class(p) <- unique(c("ggchord", class(p)))
  p
}
