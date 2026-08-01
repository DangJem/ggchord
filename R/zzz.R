# zzz.R - package environment and infrastructure
# Package-level environment for passing data, parameters, and layout between ggchord() and the geom layers

#' Package-level environment
#'
#' Internal environment that stores the chord diagram's raw data, geom parameters, and the lazily computed layout result.
#'
#' @keywords internal
.chord_env <- new.env(parent = emptyenv())

# ====================================================================
# Data storage
# ====================================================================

#' Set raw data into the package environment
#' @keywords internal
set_chord_data <- function(seq_data, ribbon_data, gene_data) {
  .chord_env$seq_data    <- seq_data
  .chord_env$ribbon_data <- ribbon_data
  .chord_env$gene_data   <- gene_data
}

#' Get raw data
#' @keywords internal
get_chord_data <- function() {
  list(
    seq_data    = .chord_env$seq_data,
    ribbon_data = .chord_env$ribbon_data,
    gene_data   = .chord_env$gene_data
  )
}

# ====================================================================
# Parameter storage (each geom stores its parameters in this environment during composition)
# ====================================================================

# Global parameters
set_global_params <- function(params) {
  .chord_env$global_params <- params
}
get_global_params <- function() {
  .chord_env$global_params
}

# Sequence parameters (from geom_seq)
set_seq_params <- function(params) {
  .chord_env$seq_params <- params
}
get_seq_params <- function() {
  .chord_env$seq_params
}

# Ribbon parameters (from geom_ribbon)
set_ribbon_params <- function(params) {
  .chord_env$ribbon_params <- params
}
get_ribbon_params <- function() {
  .chord_env$ribbon_params
}

# Gene parameters (from geom_gene)
set_gene_params <- function(params) {
  .chord_env$gene_params <- params
}
get_gene_params <- function() {
  .chord_env$gene_params
}

# Axis parameters (from geom_axis)
set_axis_params <- function(params) {
  .chord_env$axis_params <- params
}
get_axis_params <- function() {
  .chord_env$axis_params
}

# ====================================================================
# Layout storage (lazily computed at print time)
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
      "Chord layout data not found. Please call ggchord() to create the layout first.",
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
