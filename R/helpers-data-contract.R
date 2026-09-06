# Input-only compatibility bridge. Never changes the user's object in place.
ggchord_normalize_accver <- function(data) {
  if (!is.data.frame(data) || !"seq_id" %in% names(data)) return(data)
  if ("accver" %in% names(data)) {
    ggchord_stop("Use only one of accver and seq_id; rename the legacy column explicitly")
  }
  warning("seq_id is deprecated; rename it to accver (removed in v0.13.0)", call. = FALSE)
  names(data)[names(data) == "seq_id"] <- "accver"
  data
}

ggchord_require_columns <- function(data, columns, caller) {
  if (!is.data.frame(data)) ggchord_stop(caller, ": data must be a data.frame")
  missing <- setdiff(columns, names(data))
  if (length(missing)) ggchord_stop(caller, ": missing required column(s): ", paste(missing, collapse = ", "))
  invisible(data)
}
