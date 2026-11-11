# BEDPE format description (no third-party parser code copied):
# https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format

#' Read paired BED intervals as ggchord links
#'
#' Reads tab-separated BEDPE with six required columns, optional name, score,
#' strand1 and strand2 in standard order, and any additional fields. Each file
#' must have a consistent field count. Blank/comment lines and tab-free
#' track/browser headers are skipped. Accession identifiers and the order of
#' the two anchors are retained: the first anchor is query, the second subject.
#'
#' Zero-based half-open intervals become one-based inclusive intervals.
#' Negative-strand anchor endpoints decrease, matching [as_ggchord_link()] for
#' paired interval objects. Unspecified strands remain `"."`; relative anchor
#' direction is known only when both anchors and both strands are known.
#' Anchor direction does not establish a pairwise alignment. No identity,
#' alignment length, read support or full sequence length is inferred.
#'
#' Original standard fields are retained as character columns, including
#' arbitrary score strings. Extra fields are retained in the `bedpe_extra`
#' list-column as character vectors named by their physical column number.
#' Names are labels, not guaranteed unique HSP IDs. Unknown chromosome `"."`
#' and coordinate `-1` are explicit missing-value sentinels. The default stops
#' on incomplete anchors; `unmapped = "keep"` retains them for auditing with
#' missing canonical fields and `q_mapped`/`s_mapped` flags. Filter for both
#' flags before adaptation or plotting. Invalid known coordinates, including
#' zero-width intervals, are always errors.
#'
#' @param file Optional local path vector or glob patterns; read before `files`.
#' @param files Optional additional local paths or glob patterns.
#' @param source_file Logical. Add a `.source_file` column, default FALSE.
#' @param unmapped `"error"` (default) or `"keep"` for incomplete anchors.
#' @return A data frame with canonical q/s endpoints, anchor strands, mapping
#'   flags, relative direction, original BEDPE fields, `bedpe_extra` and physical
#'   `.source_line`. The `ggchord_import` attribute records files, row counts
#'   and coordinate conventions. The keep policy does not make missing anchors
#'   drawable and does not manufacture source positions.
#' @seealso [as_ggchord_link()], [ggchord_data_contract],
#'   <https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format>
#' @md
#' @export
#' @examples
#' library(dplyr)
#' data(ribbon_data_example)
#' # Re-encode observed intervals for a format round trip, not paired-read evidence.
#' paired <- ribbon_data_example |>
#'   transmute(
#'     chrom1 = qaccver, start1 = pmin(qstart, qend) - 1, end1 = pmax(qstart, qend),
#'     chrom2 = saccver, start2 = pmin(sstart, send) - 1, end2 = pmax(sstart, send)
#'   )
#' path <- tempfile(fileext = ".bedpe")
#' write.table(paired, path,
#'   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
#' )
#' read_bedpe(path)
#' unlink(path)
read_bedpe <- function(file = NULL, files = NULL, source_file = FALSE,
                       unmapped = c("error", "keep")) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  caller <- "read_bedpe()"
  unmapped <- match.arg(unmapped)
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ": source_file must be TRUE or FALSE")
  paths <- ggchord_resolve_files(file, files, "read_bedpe")
  standard <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2",
    "name", "score", "strand1", "strand2")
  tables <- lapply(paths, function(path) {
    lines <- ggchord_read_text(path, caller)
    header <- !grepl("\t", lines, fixed = TRUE) &
      grepl("^[[:space:]]*(track|browser)[[:space:]]", lines)
    kept <- which(nzchar(trimws(lines)) & !grepl("^[[:space:]]*#", lines) & !header)
    if (!length(kept)) ggchord_stop(caller, ": no data rows in ", path)
    fields <- lapply(kept, function(i) {
      # A sentinel preserves trailing empty metadata fields in strsplit().
      value <- strsplit(paste0(lines[i], "\t__end__"), "\t", fixed = TRUE)[[1L]]
      utils::head(value, -1L)
    })
    widths <- lengths(fields)
    if (length(unique(widths)) != 1L || any(widths < 6L))
      ggchord_stop(caller, ": each file needs a consistent field count of at least 6 tab-separated columns: ", path)
    nstandard <- min(10L, widths[1L])
    raw <- as.data.frame(do.call(rbind, lapply(fields, function(x) x[seq_len(nstandard)])),
      stringsAsFactors = FALSE)
    names(raw) <- standard[seq_len(nstandard)]
    bad <- function(message, row) ggchord_stop(caller, ": ", message, " in ", path,
      " at line ", kept[row])
    number <- function(name, minimum) {
      value <- raw[[name]]
      parsed <- suppressWarnings(as.numeric(value))
      invalid <- !grepl("^(-1|[0-9]+)$", value) | !is.finite(parsed) |
        (parsed != -1 & (parsed < minimum | parsed > 2^53 - 1))
      if (any(invalid)) bad(paste0(name, " must be an integral boundary or -1"), which(invalid)[1L])
      parsed
    }
    coordinate <- list(start1 = number("start1", 0), end1 = number("end1", 1),
      start2 = number("start2", 0), end2 = number("end2", 1))
    anchor <- function(index) {
      chrom <- raw[[paste0("chrom", index)]]
      start <- coordinate[[paste0("start", index)]]; end <- coordinate[[paste0("end", index)]]
      if (any(!nzchar(chrom))) bad("empty chromosome identifier", which(!nzchar(chrom))[1L])
      invalid <- start >= 0 & end >= 0 & end <= start
      if (any(invalid)) bad("end must exceed start; zero-width anchors are not intervals", which(invalid)[1L])
      strand <- raw[[paste0("strand", index)]] %||% rep(".", nrow(raw))
      if (any(!strand %in% c("+", "-", "."))) bad("strand must be '+', '-' or '.'", which(!strand %in% c("+", "-", "."))[1L])
      mapped <- chrom != "." & start >= 0 & end >= 1
      low <- ifelse(start >= 0, start + 1, NA_real_)
      high <- ifelse(end >= 1, end, NA_real_)
      list(accver = ifelse(chrom == ".", NA_character_, chrom), strand = strand, mapped = mapped,
        start = ifelse(strand == "-", high, low), end = ifelse(strand == "-", low, high))
    }
    q <- anchor(1L); s <- anchor(2L)
    incomplete <- !q$mapped | !s$mapped
    if (any(incomplete) && unmapped == "error") bad(
      "incomplete anchor; use unmapped='keep' for auditing and filter q_mapped & s_mapped before plotting",
      which(incomplete)[1L])
    known <- q$mapped & s$mapped & q$strand != "." & s$strand != "."
    out <- data.frame(qaccver = q$accver, saccver = s$accver,
      qstart = q$start, qend = q$end, sstart = s$start, send = s$end,
      qstrand = q$strand, sstrand = s$strand, q_mapped = q$mapped, s_mapped = s$mapped,
      orientation_known = known,
      direction = ifelse(known, ifelse(q$strand == s$strand, "same", "reverse"), NA_character_),
      stringsAsFactors = FALSE)
    for (name in names(raw)) out[[name]] <- raw[[name]]
    out$bedpe_extra <- I(lapply(fields, function(x) {
      if (length(x) <= 10L) return(stats::setNames(character(), character()))
      stats::setNames(x[seq.int(11L, length(x))], paste0("column_", seq.int(11L, length(x))))
    }))
    out$.source_line <- kept
    ggchord_add_source_file(out, path, source_file)
  })
  out <- ggchord_format_import(ggchord_rbind_fill(tables), paths, "BEDPE", vapply(tables, nrow, integer(1)))
  attr(out, "ggchord_import")$unmapped <- unmapped
  out
}
