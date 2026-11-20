# File coordinates belong to data. Track headers never override plot grammar.
# Specification: https://genome.ucsc.edu/goldenPath/help/bedgraph.html

#' Read bedGraph signal intervals
#'
#' Convert half-open spans to one-based closed source intervals without
#' computing values at unobserved positions.
#'
#' @param file Local bedGraph paths.
#' @param source_file Whether to include row-level source paths.
#' @return A track table with accver, start, end, value and source provenance.
#' @export
#' @examples
#' f <- tempfile(fileext = ".bedgraph")
#' writeLines(c("A	0	100	0.4", "A	100	200	0.6"), f)
#' signal_intervals <- read_bedgraph(f)
#' unlink(f)
read_bedgraph <- function(file, source_file = FALSE) {
  caller <- 'read_bedgraph()'
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, NULL, 'read_bedgraph')
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

#' Read a BigWig signal track
#'
#' Read source intervals using the optional rtracklayer decoder. Missing
#' intervals remain missing. Signal units retain their source meaning.
#' The source report retains sequence lengths and the actual selection.
#'
#' @param file Local BigWig paths.
#' @param which Optional genomic selection as a supported interval object
#'   or an accver/start/end table.
#' @param source_file Whether to include row-level source paths.
#' @return An accver/start/end/value track table with a decoder and selection report.
#' @export
#' @examples
#' if (requireNamespace("rtracklayer", quietly = TRUE)) {
#'   file <- system.file("extdata", "pUC19-gc-100.bw", package = "ggchord")
#'   read_bigwig(file)
#' }
read_bigwig <- function(file, which = NULL, source_file = FALSE) {
  caller <- 'read_bigwig()'
  if (.Platform$OS.type == 'windows') ggchord_stop(caller,
    ': optional rtracklayer BigWig decoding is unsupported on Windows; use read_bedgraph() with an existing bedGraph file')
  if (!requireNamespace('rtracklayer', quietly = TRUE)) ggchord_stop(caller,
    ": requires optional package 'rtracklayer'; install it or read an existing bedGraph with read_bedgraph()")
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, NULL, 'read_bigwig')
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
