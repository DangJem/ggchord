# BEDPE format description (no third-party parser code copied):
# https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format

#' Read BEDPE interval pairs
#'
#' Convert each placed anchor from the file's zero-based half-open convention
#' to ordered one-based closed endpoints. Preserve strand evidence and the
#' original score. The input's first anchor is query and its second is
#' subject. Unplaced rows are recorded in the data report.
#'
#' @param file Local BEDPE paths.
#' @param sequence_unit The native coordinate unit, `"bp"` or `"aa"`.
#' @param source_file Whether to include a row-level source path column.
#' @return An interval relationship table. [data_report()] records source
#'   paths, original fields, and unplaced records.
#' @export
#' @examples
#' f <- tempfile(fileext = ".bedpe")
#' writeLines("A	99	240	B	149	290	pair1	10	+	-", f)
#' paired <- read_bedpe(f)
#' unlink(f)
read_bedpe <- function(file, sequence_unit = "bp", source_file = FALSE) {
  caller <- "read_bedpe()"
  sequence_unit <- ggchord_match_arg(sequence_unit, c("bp", "aa"))
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ": source_file must be TRUE or FALSE")
  paths <- ggchord_resolve_files(file, NULL, "read_bedpe")
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
      head(value, -1L)
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
    known <- q$mapped & s$mapped & q$strand != "." & s$strand != "."
    out <- data.frame(qaccver = q$accver, saccver = s$accver,
      qstart = q$start, qend = q$end, sstart = s$start, send = s$end,
      qstrand = q$strand, sstrand = s$strand, q_mapped = q$mapped, s_mapped = s$mapped,
      orientation_known = known,
      direction = ifelse(known, ifelse(q$strand == s$strand, "same", "reverse"), NA_character_),
      stringsAsFactors = FALSE)
    for (name in names(raw)) out[[name]] <- raw[[name]]
    out$bedpe_extra <- I(lapply(fields, function(x) {
      if (length(x) <= 10L) return(setNames(character(), character()))
      setNames(x[seq.int(11L, length(x))], paste0("column_", seq.int(11L, length(x))))
    }))
    out$.source_line <- kept
    ggchord_add_source_file(out, path, source_file)
  })
  out <- ggchord_format_import(ggchord_rbind_fill(tables), paths, "BEDPE", vapply(tables, nrow, integer(1)))
  out <- ggchord_order_link_endpoints(out, caller)
  out$sequence_unit <- rep(sequence_unit, nrow(out))
  placed <- out$q_mapped & out$s_mapped
  attr(out, "ggchord_import")$unplaced <- out[!placed, , drop = FALSE]
  out[placed, , drop = FALSE]
}
