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
