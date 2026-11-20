# Shared primitives for interval relationship data.

ggchord_ribbon_required_columns <- function() {
  c("qaccver", "saccver", "qstart", "qend", "sstart", "send")
}

ggchord_typed_na <- function(x) {
  if (is.factor(x)) return(factor(NA_character_, levels = levels(x)))
  if (inherits(x, "Date")) return(as.Date(NA))
  if (inherits(x, "POSIXct")) return(as.POSIXct(NA_real_, origin = "1970-01-01"))
  x[NA_integer_][1]
}
