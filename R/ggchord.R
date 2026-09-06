# ggchord.R - constructor
# Entry point of the layered ggplot2 API for chord diagrams
# v0.6.0: plot objects are self-contained. Data and parameters are stored on
# the plot object itself, and the layout is computed at build time so that
# print(), ggsave(), ggplot_build() and other ggplot2 workflows all work.

# Global variable declarations (to avoid R CMD check NOTEs)
globalVariables(c(
  "x", "y", "group", "pident", "fill", "colour", "strand", "anno", "accver",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "xend", "yend", ".component",
  "label", "label_x", "label_y", "size", "angle",
  "fill_col", "alpha", "label_hjust", "label_vjust", "label_angle",
  "linetype", "zcolour", "zregionfill", "zoutline", "zlinetype",
  "outline_col", "linetype_val", "value", "source_row", "direction",
  "seq_colour", "ribbon_fill", "ribbon_alpha",
  "ribbon_colour", "ribbon_linetype", "gene_fill", "feature_fill",
  "feature_shape", "seq_ring", "bundle_n", "bundle_weight", "density",
  "region_fill", ".bundle_n", ".bundle_weight", ".bundle_density"
))

#' ggchord: layered multi-sequence alignment chord diagrams for ggplot2
#'
#' ggchord visualizes multi-sequence alignment results using ggplot2's layered grammar.
#' The \code{ggchord()} constructor handles data validation and global settings;
#' the \code{geom_*} layers are stacked as needed, each responsible for its own layout parameters and visual rendering.
#' The layout is computed lazily when the plot is built (e.g. via \code{print()},
#' \code{ggsave()}, or \code{ggplot_build()}).
#'
#' @param seq_data data.frame/tibble with accver and length.
#' @param ribbon_data Optional alignment table with qaccver, saccver, qstart,
#'   qend, sstart and send. length and pident are optional unless requested
#'   by a mapping, filter or statistic.
#' @param gene_data Optional table with accver, start, end and strand.
#'   Annotation (anno) may be omitted.
#' @param debug Logical. Whether to output debug information, default FALSE
#' @param validate Character, default \code{"warn"}. How to run the structured
#'   input-data validation (see \code{\link{validate_ggchord_data}}):
#'   \code{"warn"} emits a single summary warning when the data has problems
#'   and caches the full report on the plot object
#'   (\code{p$ggchord$validation}); \code{"error"} stops on severe problems;
#'   \code{"none"} skips the diagnostic validation (the cheap structural
#'   checks that prevent crashes are still performed).
#' @param ... Reserved for clear migration errors from removed constructor
#'   arguments. Use [labs()], [coord_chord()] and [theme()] for presentation.
#'
#' @return A ggchord object (inherits from ggplot) to which geom_* layers can be added with +
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' data(gene_data_example)
#'
#' p <- ggchord(
#'   seq_data = seq_data_example,
#'   ribbon_data = ribbon_data_example,
#'   gene_data = gene_data_example
#' ) +
#'   geom_seq() +
#'   geom_link_ribbon() +
#'   geom_gene()
#' print(p)
#'
#' @importFrom ggplot2 ggplot_build
ggchord <- function(
    seq_data,
    ribbon_data = NULL,
    gene_data = NULL,
    debug = FALSE,
    validate = c("warn", "error", "none"),
    ...
) {
  seq_data <- ggchord_normalize_accver(seq_data)
  gene_data <- ggchord_normalize_accver(gene_data)

  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
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
    seq_data, ribbon_data, gene_data, strict = FALSE,
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
  if (!is.null(gene_data)) {
    if (nrow(gene_data) == 0) warning("No valid gene annotation data in gene_data")
    if (debug) cat("Number of gene annotation rows: ", nrow(gene_data), "\n")
  }

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
  p <- ggplot2::ggplot() +
    coord_chord(rotation = 45) +
    theme_ggchord()

  p$ggchord <- list(
    data   = list(seq_data = seq_data, ribbon_data = ribbon_data,
                  gene_data = gene_data),
    global = list(rotation = 45, debug = debug, validate = validate),
    validation = validation,
    # Shared reference environment for plot-owned layout caching.
    ref    = new.env(),
    layout = NULL
  )
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
  if (inherits(e2, "Scale")) {
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
  # Invalidate the plot-owned layout after adding any component.
  if (!is.null(e1$ggchord$ref)) e1$ggchord$ref$layout <- NULL
  attr(p$theme, "ggchord.settings") <- new_settings %||% old_settings
  class(p) <- unique(c("ggchord", class(p)))
  p
}
