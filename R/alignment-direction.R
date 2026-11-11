# Shared direction validation for new HSP projections. Endpoints are already
# normalized by import; this helper never manufactures a biological direction.
ggchord_alignment_unknown <- function(input, caller) {
  unknown <- rep(FALSE, nrow(input))
  for (field in intersect(c("direction", "strand", "relative_strand"), names(input))) {
    value <- tolower(as.character(input[[field]]))
    unknown <- unknown | is.na(value) | value %in% c("unknown", ".", "?", "*", "0", "")
  }
  if ("orientation_known" %in% names(input))
    unknown <- unknown | is.na(input$orientation_known) | !input$orientation_known
  slope <- sign(input$qend - input$qstart) * sign(input$send - input$sstart)
  for (field in intersect(c("direction", "strand", "relative_strand"), names(input))) {
    value <- tolower(as.character(input[[field]]))
    if (any(!unknown & !value %in% c("same", "forward", "reverse", "+", "-", "1", "-1")))
      ggchord_stop(caller, ": unrecognized relative direction in ", field)
    reverse <- value %in% c("reverse", "-", "-1")
    same <- value %in% c("same", "forward", "+", "1")
    if (any(!unknown & ((reverse & slope > 0) | (same & slope < 0))))
      ggchord_stop(caller, ": direction conflicts with endpoints; normalize reverse endpoints first (read_paf() does this)")
  }
  unknown
}
