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
#' data(seq_data_example)
#' data(ribbon_data_example)
#' data(gene_data_example)
#' seqs <- seq_data_example[c("accver", "length")]
#' links <- ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
#' # The fixture preserves biological strands separately from its old display strands.
#' genes <- transform(gene_data_example, strand = source_strand)[c("accver", "start", "end", "strand")]
#' # Optional appearance choices; none of these columns or scales is required input.
#' seq_palette <- setNames(c("#487A70", "#B48A61", "#637E9D", "#A16E7D"), seqs$accver)
#' gene_palette <- c("+" = "#56887C", "-" = "#B58169", "." = "#AAB9AE")
#' feature_palette <- c(CDS = "#6E9A8A", promoter = "#8798B2", rep_origin = "#C4A46A",
#'   primer_bind = "#A28EAA", misc_feature = "#CDD5CE", protein_bind = "#A9B9C4")
#' doc_theme <- theme_ggchord_publication(axis.text = element_text(size = 9),
#'   axis.text.check.overlap = TRUE) + theme(legend.position = "bottom")
#' data(plasmid_example_pUC19)
#' track_seq <- plasmid_example_pUC19[c("accver", "length")]
#' starts <- seq.int(1L, track_seq$length, by = 100L)
#' bins <- data.frame(accver = track_seq$accver, start = starts,
#'   end = pmin(starts + 99L, track_seq$length))
#' bins$gc <- vapply(seq_len(nrow(bins)), function(i) {
#'   bases <- strsplit(substr(plasmid_example_pUC19$sequence, bins$start[i], bins$end[i]), "")[[1]]
#'   mean(bases %in% c("G", "C"))
#' }, numeric(1))
#' p_tracks <- ggchord(seq_data = track_seq) + geom_seq() + geom_seq_center_label() +
#'   geom_track_bar(aes(track_value = gc, track_fill = gc), data = bins,
#'     position = position_track("gc", width = unit(6, "mm"))) +
#'   scale_track_value_continuous("gc", limits = c(0, 1), unit = "GC fraction") +
#'   scale_track_fill_gradientn("gc", colours = c("#EDF1EA", "#8AAE9C", "#345F58")) +
#'   coord_circular() + doc_theme
#' loci <- data.frame(accver = track_seq$accver, start = 501, end = 2000)
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
      row$.truncated_left <- isTRUE(row$.truncated_left) || a < lo
      row$.truncated_right <- isTRUE(row$.truncated_right) || b > hi
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
