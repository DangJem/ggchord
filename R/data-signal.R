# File coordinates belong to data. Track headers never override plot grammar.
# Specification: https://genome.ucsc.edu/goldenPath/help/bedgraph.html

#' Read bedGraph intervals as quantitative track data
#'
#' Reads four columns, converts zero-based half-open coordinates to one-based
#' inclusive positions, and retains every finite signed value and source row.
#' Missing intervals remain absent: no zero filling, interpolation, sorting or
#' aggregation is performed. Overlapping rows are retained as observations,
#' not silently resolved into a single value. Choose any summarisation explicitly.
#'
#' A file contains one signal. More than one track definition is rejected;
#' split different browser tracks explicitly rather than mixing their domains.
#' Comments and display headers are retained in import provenance. Their colour,
#' smoothing and other display instructions never configure a ggchord plot.
#' Sequence identifiers, including literal track/browser, are retained verbatim.
#' Neither value units nor sequencing coverage are inferred from the file name.
#'
#' @param file Optional input path vector or glob patterns, before files.
#' @param files Optional additional paths or wildcard patterns for bedGraph files.
#' @param source_file Logical. Add .source_file to each row, default FALSE.
#' @return A data frame with accver, start, end, value, original chromStart,
#'   chromEnd, dataValue and physical .source_line. The ggchord_import attribute
#'   records conventions, per-file row counts, and unchanged headers. When
#'   combining different signals, keep source_file = TRUE and select each
#'   signal's layers/scales explicitly. Coordinate units follow source metadata.
#' @seealso [geom_track()], [scale_track_value_continuous()]
#' @export
#' @examples
#' # Actual 100-bp-window GC percentages; not measured read coverage.
#' path <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
#' signal <- read_bedgraph(path)
#' head(signal)
read_bedgraph <- function(file = NULL, files = NULL, source_file = FALSE) {
  caller <- 'read_bedgraph()'
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, files, 'read_bedgraph')
  input <- lapply(paths, function(path) {
    parsed <- ggchord_format_lines(path, caller, bed = TRUE)
    if (any(lengths(parsed$fields) != 4L)) ggchord_stop(caller,
      ': each data row needs exactly four fields: ', path)
    if (sum(grepl('^track([[:space:]]|$)', trimws(parsed$headers$text))) > 1L)
      ggchord_stop(caller, ': each file must contain one signal; split multiple track definitions explicitly: ', path)
    fields <- do.call(rbind, parsed$fields)
    start <- ggchord_format_integer(fields[, 2], 'chromStart', caller)
    end <- ggchord_format_integer(fields[, 3], 'chromEnd', caller)
    if (any(!nzchar(fields[, 1])) || any(end <= start))
      ggchord_stop(caller, ': require nonempty identifiers and chromEnd > chromStart: ', path)
    value <- suppressWarnings(as.numeric(fields[, 4]))
    number <- '^[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][+-]?[0-9]+)?$'
    bad <- which(!grepl(number, fields[, 4]) | !is.finite(value))
    if (length(bad)) ggchord_stop(caller, ': dataValue must be a finite real number in ',
      path, ' at line ', parsed$line[bad[1L]])
    out <- data.frame(accver = fields[, 1], start = start+1, end = end, value = value,
      chromStart = start, chromEnd = end, dataValue = value,
      .source_line = parsed$line, stringsAsFactors = FALSE)
    list(data = ggchord_add_source_file(out, path, source_file), headers = parsed$headers)
  })
  tables <- lapply(input, `[[`, 'data')
  out <- ggchord_format_import(ggchord_rbind_fill(tables), paths, 'bedGraph', vapply(tables, nrow, integer(1)))
  provenance <- attr(out, 'ggchord_import')
  provenance$headers <- stats::setNames(lapply(input, `[[`, 'headers'), paths)
  provenance$values <- 'unchanged finite interval values; no imputation, aggregation or inferred units'
  attr(out, 'ggchord_import') <- provenance
  out
}

# Validate selection in canonical coordinates; the decoder owns binary parsing.
ggchord_bigwig_selection <- function(which, information, caller) {
  if (is.null(which)) return(NULL)
  if (inherits(which, 'GRanges')) {
    query <- which
    d <- data.frame(accver = as.character(GenomicRanges::seqnames(query)),
      start = as.numeric(GenomicRanges::start(query)), end = as.numeric(GenomicRanges::end(query)))
  } else if (is.data.frame(which)) {
    ggchord_require_columns(which, c('accver', 'start', 'end'), caller)
    d <- which[c('accver', 'start', 'end')]
    query <- NULL
  } else ggchord_stop(caller, ': which must be NULL, a GRanges, or an accver/start/end data frame')
  if (!is.character(d$accver) || anyNA(d$accver) || any(!nzchar(d$accver)) ||
      !is.numeric(d$start) || !is.numeric(d$end) || any(!is.finite(c(d$start, d$end))) ||
      any(d$start < 1 | d$end < d$start) || any(c(d$start, d$end) != floor(c(d$start, d$end))) ||
      any(d$end > .Machine$integer.max)) ggchord_stop(caller,
        ': queries need ordered positive integer coordinates within the GRanges range')
  at <- match(d$accver, information$accver)
  if (anyNA(at)) ggchord_stop(caller, ': query accession absent from the BigWig header: ',
    paste(unique(d$accver[is.na(at)]), collapse = ', '))
  if (any(d$end > information$length[at])) ggchord_stop(caller, ': query exceeds a stored sequence length')
  if (is.null(query)) query <- GenomicRanges::GRanges(d$accver,
    ranges = IRanges::IRanges(as.integer(d$start), as.integer(d$end)))
  GenomicRanges::reduce(query, ignore.strand = TRUE)
}

#' Read indexed BigWig signals through an optional mature decoder
#'
#' Uses optional rtracklayer and its GRanges import. Those coordinates are
#' already one-based inclusive and are not shifted again. Stored sequence
#' lengths come from the file header, not the last nonmissing interval.
#' Scores retain the binary file's floating-point representation. They are
#' signal values, not original observation counts, identity or read support.
#'
#' Optional selection is an unstranded coordinate union, so overlapping query
#' windows do not duplicate a source signal. Returned runs are clipped to
#' that union. The original query and effective union are retained
#' separately; no pre-clipping window boundaries are invented. Missing ranges
#' stay absent, and imported missing scores remain missing. No zero filling,
#' smoothing, summarisation or value-unit inference is performed.
#' The optional decoder does not support BigWig on Windows; use an existing
#' bedGraph file with read_bedgraph() there instead.
#'
#' @param file Optional local input path vector or glob patterns, before files.
#' @param files Optional additional paths or wildcards for local BigWig files.
#' @param which NULL for the whole file, or a GRanges / data frame with
#'   character accver and numeric start/end in one-based inclusive coordinates.
#'   Unknown accessions and out-of-bounds intervals fail explicitly. An empty
#'   selection returns no signal rows. The same selection applies to each file.
#' @param source_file Logical. Add .source_file to each row, default FALSE.
#' @return A data frame with accver, start, end, value and original imported
#'   score. Import provenance includes decoder/version, each file's sequence
#'   dictionary, original query and effective selected intervals. No physical
#'   source lines or original measurement IDs are manufactured from a binary file.
#' @seealso [read_bedgraph()], [geom_track()]
#' @export
#' @examples
#' if (.Platform$OS.type != "windows" && requireNamespace("rtracklayer", quietly = TRUE)) {
#'   # Packaged real GC signal, not measured read coverage.
#'   path <- system.file("extdata", "pUC19-gc-100.bw", package = "ggchord")
#'   signal <- read_bigwig(path)
#'   local <- read_bigwig(path,
#'     which = data.frame(accver = "pUC19", start = 50, end = 150))
#' }
read_bigwig <- function(file = NULL, files = NULL, which = NULL, source_file = FALSE) {
  caller <- 'read_bigwig()'
  if (.Platform$OS.type == 'windows') ggchord_stop(caller,
    ': optional rtracklayer BigWig decoding is unsupported on Windows; use read_bedgraph() with an existing bedGraph file')
  if (!requireNamespace('rtracklayer', quietly = TRUE)) ggchord_stop(caller,
    ": requires optional package 'rtracklayer'; install it or read an existing bedGraph with read_bedgraph()")
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, files, 'read_bigwig')
  input <- lapply(paths, function(path) {
    resource <- rtracklayer::BigWigFile(path)
    information <- ggchord_conversion_seqinfo(GenomicRanges::seqinfo(resource), caller)
    if (any(!is.finite(information$length) | information$length < 1 | information$length != floor(information$length)))
      ggchord_stop(caller, ': decoder cannot supply complete positive sequence lengths from the BigWig header')
    selection <- ggchord_bigwig_selection(which, information, caller)
    ranges <- if (is.null(selection)) rtracklayer::import(resource, as = 'GRanges') else
      if (!length(selection)) GenomicRanges::GRanges() else
        rtracklayer::import(resource, which = selection, as = 'GRanges')
    scores <- S4Vectors::mcols(ranges)$score
    if (length(ranges) && !is.numeric(scores)) ggchord_stop(caller, ': decoder did not supply numeric scores')
    if (!length(ranges)) scores <- numeric()
    # Enforce the requested union independently of decoder-version clipping.
    # Each fragment retains its imported score, with no new observation ID.
    if (!is.null(selection) && length(ranges)) {
      hits <- GenomicRanges::findOverlaps(ranges, selection, ignore.strand = TRUE)
      source <- S4Vectors::queryHits(hits)
      target <- S4Vectors::subjectHits(hits)
      clipped_start <- pmax(GenomicRanges::start(ranges)[source], GenomicRanges::start(selection)[target])
      clipped_end <- pmin(GenomicRanges::end(ranges)[source], GenomicRanges::end(selection)[target])
      ranges <- ranges[source]
      GenomicRanges::ranges(ranges) <- IRanges::IRanges(clipped_start, clipped_end)
      scores <- scores[source]
    }
    out <- data.frame(accver = as.character(GenomicRanges::seqnames(ranges)),
      start = as.numeric(GenomicRanges::start(ranges)), end = as.numeric(GenomicRanges::end(ranges)),
      value = scores, score = scores, stringsAsFactors = FALSE)
    chosen <- if (is.null(selection)) NULL else data.frame(
      accver = as.character(GenomicRanges::seqnames(selection)),
      start = as.numeric(GenomicRanges::start(selection)), end = as.numeric(GenomicRanges::end(selection)))
    list(data = ggchord_add_source_file(out, path, source_file), information = information, selection = chosen)
  })
  tables <- lapply(input, `[[`, 'data')
  out <- ggchord_format_import(ggchord_rbind_fill(tables), paths, 'BigWig', vapply(tables, nrow, integer(1)))
  provenance <- attr(out, 'ggchord_import')
  provenance$decoder <- list(package = 'rtracklayer', version = as.character(utils::packageVersion('rtracklayer')))
  provenance$adapter_coordinates <- 'rtracklayer GRanges already 1-based-inclusive; no further shift'
  provenance$sequence_info <- stats::setNames(lapply(input, `[[`, 'information'), paths)
  provenance$query_input <- which
  provenance$selected_intervals <- stats::setNames(lapply(input, `[[`, 'selection'), paths)
  provenance$values <- 'native imported floating-point signal; no imputation, aggregation or inferred units'
  attr(out, 'ggchord_import') <- provenance
  out
}
