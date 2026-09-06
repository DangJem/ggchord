# Reject the v0.12 compatibility name at every public data boundary.
ggchord_normalize_accver <- function(data) {
  if (!is.data.frame(data) || !"seq_id" %in% names(data)) return(data)
  ggchord_stop("`seq_id` was removed in v0.13.0; rename the column to `accver`")
}

ggchord_require_columns <- function(data, columns, caller) {
  if (!is.data.frame(data)) ggchord_stop(caller, ": data must be a data.frame")
  missing <- setdiff(columns, names(data))
  if (length(missing)) ggchord_stop(caller, ": missing required column(s): ", paste(missing, collapse = ", "))
  invisible(data)
}
