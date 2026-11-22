interval_recip_overlap <- function(a1, a2, b1, b2) {
  lo <- max(min(a1, a2), min(b1, b2))
  hi <- min(max(a1, a2), max(b1, b2))
  ovl <- max(0, hi - lo + 1)
  denom <- max(abs(a2 - a1) + 1, abs(b2 - b1) + 1)
  if (denom <= 0) return(0)
  ovl / denom
}

#' @export
print.ggchord_validation <- function(x, ...) {
  ggchord_local_console_errors()
  cat(sprintf("ggchord data validation: %s (%d errors, %d warnings)\n",
    if (x$valid) "valid" else "invalid", x$summary[["errors"]], x$summary[["warnings"]]))
  if (nrow(x$issues)) print(x$issues, row.names = FALSE)
  invisible(x)
}

#' @export
summary.ggchord_validation <- function(object, ...) {
  ggchord_local_console_errors()
  object$summary
}

#' @export
as.data.frame.ggchord_validation <- function(x, row.names = NULL, optional = FALSE, ...) {
  ggchord_local_console_errors()
  x$issues
}
