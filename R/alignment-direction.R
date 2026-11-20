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
      ggchord_stop(caller, ": direction conflicts with endpoints; use the original directional endpoints. When cleaning those records, set reversed_interval = 'keep'.")
  }
  unknown
}

# Relative direction is separate from the canonical ascending interval endpoints.
ggchord_reference_direction <- function(input, caller) {
  direction <- rep(NA_real_, nrow(input))
  merge_evidence <- function(evidence) {
    same <- !is.na(direction) & !is.na(evidence)
    if (any(direction[same] != evidence[same]))
      ggchord_stop(caller, ": conflicting link direction", ggchord_input_rows(same & direction != evidence))
    direction[!is.na(evidence)] <<- evidence[!is.na(evidence)]
  }
  signs <- c("+" = 1, "-" = -1)
  if (all(c("qstrand", "sstrand") %in% names(input))) {
    merge_evidence(unname(signs[as.character(input$qstrand)] * signs[as.character(input$sstrand)]))
  }
  codes <- c(same = 1, forward = 1, reverse = -1, opposite = -1, "+" = 1, "-" = -1, "1" = 1, "-1" = -1)
  for (key in intersect(c("direction", "relative_strand"), names(input))) {
    raw <- tolower(as.character(input[[key]]))
    unknown <- is.na(raw) | raw %in% c("unknown", ".", "?", "*", "0", "")
    if (any(!unknown & !raw %in% names(codes))) ggchord_stop(caller, ": invalid ", key)
    merge_evidence(unname(codes[raw]))
  }
  if (!is.null(input$orientation_known)) direction[!input$orientation_known %in% TRUE] <- NA_real_
  direction
}
