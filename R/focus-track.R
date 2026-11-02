#' Focus quantitative intervals or samples without bridging omitted windows
#'
#' Applies the same explicit loci and output IDs as [focus_ggchord_data()] and
#' [segment_ggchord_data()]. It never averages, interpolates, rescales values,
#' or creates boundary samples. Missing values remain missing. Subsequent
#' track layers group by the returned window accession, so separate windows
#' cannot be joined by a line. Source coordinates and original row lineage
#' survive repeated focusing.
#' @param data Quantitative input table. Metadata and mapped values are retained.
#' @param seq_data Sequence reference table.
#' @param loci Explicit `accver`, `start`, `end` and optional `locus_id`, as in
#'   [focus_ggchord_data()].
#' @param type `"point"` requires `pos`; `"interval"` requires `start`/`end`.
#' @param boundary `"trim"` (default) or `"drop"` for partly retained intervals.
#' @return An ordinary data frame with focused coordinates. Its
#'   `ggchord_track_focus` attribute reports output-to-input rows and omitted
#'   input rows; `.source_*` columns retain earliest source coordinates.
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' library(ggchord)
#' library(dplyr)
#' data(plasmid_example_pUC19)
#' track_seq <- plasmid_example_pUC19 |>
#'   select(accver, length)
#' starts <- seq.int(1L, track_seq$length, by = 100L)
#' gc_fraction <- function(sequence) {
#'   bases <- strsplit(sequence, "")[[1]]
#'   mean(bases %in% c("G", "C"))
#' }
#' bins <- tibble(
#'   accver = track_seq$accver,
#'   start = starts,
#'   end = pmin(starts + 99L, track_seq$length)
#' ) |>
#'   rowwise() |>
#'   mutate(gc = gc_fraction(substr(plasmid_example_pUC19$sequence, start, end))) |>
#'   ungroup()
#' loci <- tibble(accver = track_seq$accver, start = 501, end = 2000)
#' focus_ggchord_track(bins, track_seq, loci, type = "interval")
focus_ggchord_track <- function(data, seq_data, loci, type = c("point", "interval"),
    boundary = c("trim", "drop")) {
  type <- match.arg(type); boundary <- match.arg(boundary)
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
  windows <- ggchord_focus_validate_loci(seq_data, loci, 0, caller)
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
  out
}
