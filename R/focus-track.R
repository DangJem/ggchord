
#' Crop a quantitative track to explicit windows
#'
#' Local references use the window IDs. Points are retained only when they
#' are present in the window. Signals keep their original values and
#' measurement units. Earliest source coordinates are available through
#' the source columns and [data_report()].
#'
#' @param data Track data frame.
#' @param sequences Complete reference dictionary.
#' @param loci Window table with accver, start, end and distinct locus_id.
#' @param type `"point"` uses pos. `"interval"` uses start and end.
#' @param boundary `"trim"` intersects intervals with the window, and
#'   `"drop"` keeps whole observations only.
#' @return A focused track table with source positions and a data report.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' signal <- data.frame(accver = rep(c("A", "B"), each = 3),
#'   pos = rep(c(100, 400, 700), 2), value = c(2, 5, 3, 1, 4, 2))
#'
#' loci <- data.frame(accver="A", start=50, end=500, locus_id="A_local")
#' local_signal <- focus_ggchord_track(signal, seqs, loci, type="point")
#' data_report(local_signal)$lineage
focus_ggchord_track <- function(data, sequences, loci, type = "point", boundary = "trim") {
  seq_data <- sequences
  type <- ggchord_match_arg(type, c("point", "interval"))
  boundary <- ggchord_match_arg(boundary, c("trim", "drop"))
  caller <- "focus_ggchord_track()"
  seq_data <- as_ggchord_seq(seq_data)
  if (!is.data.frame(data)) ggchord_stop(caller, ": data must be a data.frame")
  coords <- if (type == "point") "pos" else c("start", "end")
  ggchord_require_columns(data, c("accver", coords), caller)
  ids <- ggchord_conversion_ids(data$accver, caller)
  ggchord_conversion_coordinates(data, coords, caller)
  at <- match(ids, seq_data$accver)
  if (anyNA(at) || any(vapply(coords, function(name) any(data[[name]] > seq_data$length[at]), logical(1))))
    ggchord_stop(caller, ": coordinates must lie within known sequence lengths")
  if (type == "interval" && any(data$start > data$end)) ggchord_stop(caller,
    ": split origin-crossing intervals explicitly before linear window selection")
  windows <- ggchord_explicit_focus_loci(seq_data, loci, 0, caller)
  source <- ggchord_focus_source(seq_data)
  pieces <- list(); lineage <- integer()
  for (i in seq_len(nrow(data))) for (j in which(windows$.source_accver == ids[i])) {
    lo <- windows$.focus_start[j]; hi <- windows$.focus_end[j]
    a <- if (type == "point") data$pos[i] else data$start[i]
    b <- if (type == "point") a else data$end[i]
    keep <- if (boundary == "drop") a >= lo && b <= hi else b >= lo && a <= hi
    if (!keep) next
    row <- data[i, , drop = FALSE]
    row$accver <- windows$.output_accver[j]
    row$.source_accver <- ggchord_focus_original(row, ".source_accver", source$accver[at[i]])
    row$.source_row <- ggchord_focus_original(row, ".source_row", i)
    row$.source_origin <- source$origin[at[i]] + lo - 1
    if (type == "point") {
      row$.source_pos <- ggchord_focus_original(row, ".source_pos", a + source$origin[at[i]] - 1)
      row$pos <- a - lo + 1
    } else {
      row$.source_start <- ggchord_focus_original(row, ".source_start", a + source$origin[at[i]] - 1)
      row$.source_end <- ggchord_focus_original(row, ".source_end", b + source$origin[at[i]] - 1)
      row$start <- max(a, lo) - lo + 1; row$end <- min(b, hi) - lo + 1
      row$.truncated_left <- isTRUE(row[[".truncated_left"]]) || a < lo
      row$.truncated_right <- isTRUE(row[[".truncated_right"]]) || b > hi
    }
    pieces[[length(pieces) + 1L]] <- row
    lineage <- c(lineage, i)
  }
  out <- if (length(pieces)) ggchord_rbind_fill(pieces) else data[0, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "ggchord_track_focus") <- list(type = type, input_rows = nrow(data),
    output_input_rows = lineage, omitted_rows = setdiff(seq_len(nrow(data)), lineage),
    windows = windows, source_focus = attr(data, "ggchord_track_focus"))
  out$source_accver <- out$.source_accver %||% character()
  if (type == "point") out$source_position <- out$.source_pos %||% numeric()
  else {
    out$source_start <- out$.source_start %||% numeric()
    out$source_end <- out$.source_end %||% numeric()
  }
  ggchord_attach_data_report(out, "focus_ggchord_track", inputs = data_report(data),
    lineage = data.frame(output_row = seq_along(lineage), input_row = lineage),
    metadata = attr(out, "ggchord_track_focus"))
}
