# zzz.R - package environment and infrastructure
# The package keeps no global state that affects rendering: plot data and
# parameters are stored on the plot object itself. The package environment
# only holds a cache of the most recently computed layout so that the
# get_chord_layout() accessor can inspect it after rendering.

#' Package-level environment
#'
#' Internal environment that caches the most recently computed chord layout.
#'
#' @keywords internal
.chord_env <- new.env(parent = emptyenv())

# ====================================================================
# Layout cache (set at build time; used by the get_chord_layout() accessor)
# ====================================================================

#' Set the chord layout into the package environment
#' @keywords internal
set_chord_layout <- function(layout) {
  .chord_env$layout <- layout
}

#' Get the chord layout from the package environment
#' @keywords internal
get_chord_layout <- function() {
  layout <- .chord_env$layout
  if (is.null(layout)) {
    stop(
      "Chord layout data not found. Please render the plot first.",
      call. = FALSE
    )
  }
  layout
}

# ====================================================================
# Environment cleanup
# ====================================================================

#' Clear the package environment (used to reset state)
#' @keywords internal
clear_chord_env <- function() {
  rm(list = ls(.chord_env, all.names = TRUE), envir = .chord_env)
}

# ====================================================================
# Helper operators
# ====================================================================

#' NULL coalescing operator
#'
#' Returns \code{y} if \code{x} is NULL, otherwise returns \code{x}.
#'
#' @param x Any R object (may be NULL)
#' @param y Default value returned when \code{x} is NULL
#' @name null-coalescing-operator
#' @rdname null-coalescing-operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
