# Explicit zero-based file adapters. Specifications (not copied parser code):
# https://genome.ucsc.edu/FAQ/FAQformat.html#format1
# https://github.com/lh3/miniasm/blob/master/PAF.md

ggchord_format_lines <- function(file, caller, bed = FALSE) {
  lines <- readLines(file, warn = FALSE)
  keep <- nzchar(trimws(lines)) & !grepl("^[[:space:]]*#", lines)
  if (bed) keep <- keep & !grepl("^[[:space:]]*(track|browser)([[:space:]]|$)", lines)
  index <- which(keep)
  if (!length(index)) ggchord_stop(caller, ": no data rows in ", file)
  values <- lapply(lines[index], function(line) {
    if (grepl("\t", line, fixed = TRUE)) strsplit(line, "\t", fixed = TRUE)[[1L]] else {
      if (!bed) ggchord_stop(caller, ": PAF fields must be tab-separated in ", file)
      strsplit(trimws(line), "[[:space:]]+")[[1L]]
    }
  })
  list(fields = values, line = index)
}

ggchord_format_integer <- function(x, name, caller, minimum = 0) {
  value <- suppressWarnings(as.numeric(x))
  if (any(!grepl("^[0-9]+$", x)) || any(!is.finite(value)) ||
      any(value < minimum | value > 2^53 - 1)) {
    ggchord_stop(caller, ": ", name, " must contain integers >= ", minimum,
      " within the exactly representable numeric range")
  }
  value
}

ggchord_format_import <- function(out, paths, format, rows_per_file,
                                   input_coordinates = "0-based-half-open") {
  rownames(out) <- NULL
  attr(out, "ggchord_import") <- list(format = format, files = paths,
    rows_per_file = rows_per_file, input_coordinates = input_coordinates,
    output_coordinates = if (input_coordinates == "sequence-lengths") "sequence-lengths" else "1-based-inclusive")
  out
}

#' Read BED intervals as ggchord features
#'
#' Reads standard BED3 through BED12. Starts are converted from zero-based
#' half-open intervals to one-based inclusive intervals; accession identifiers
#' are kept verbatim. Blank lines, comments, and custom-track header lines are
#' skipped. Every data row in a file must have the same number of fields.
#' Zero-width insertion records cannot be represented as nonempty genomic
#' intervals and produce an explicit error.
#'
#' BED names are labels, not guaranteed unique feature IDs. Unspecified strand
#' remains `"."` and `directionality = "none"`; use [geom_feature()] for these
#' records. Optional fields retain their BED names and coordinate conventions,
#' including zero-based `thickStart`, `thickEnd`, and relative `blockStarts`.
#' BED12 additionally creates a `segments` list-column of one-based inclusive
#' block intervals after checking block counts, boundaries and non-overlap.
#' Blocks are retained in genomic order on both strands. No colour or shape is
#' inferred from BED display fields. BED extensions beyond column 12 require an
#' explicit data-frame conversion with [as_ggchord_feature()].
#'
#' @param file Optional path to one BED file.
#' @param files Optional paths or wildcard patterns for multiple BED files.
#' @param source_file Logical. Add a `.source_file` column, default `FALSE`.
#' @return A feature data frame with `accver`, `start`, `end`, `strand`,
#'   `directionality`, original `chromStart`/`chromEnd`, optional BED fields,
#'   and `.source_line` (the physical file line). A `ggchord_import` attribute
#'   records the input/output coordinate conventions and files.
#' @seealso <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
#' @md
#' @export
#' @examples
#' bed <- tempfile(fileext = ".bed")
#' writeLines("host\t0\t100\tregulatory region\t0\t.", bed)
#' features <- read_bed(bed)
#' features[, c("accver", "start", "end", "strand")]
#' unlink(bed)
read_bed <- function(file = NULL, files = NULL, source_file = FALSE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  caller <- "read_bed()"
  paths <- ggchord_resolve_files(file, files, "read_bed")
  tables <- lapply(paths, function(path) {
    input <- ggchord_format_lines(path, caller, bed = TRUE)
    widths <- lengths(input$fields)
    if (length(unique(widths)) != 1L || any(widths < 3L | widths > 12L)) {
      ggchord_stop(caller, ": each file must have a consistent 3 to 12 BED fields: ", path)
    }
    raw <- as.data.frame(do.call(rbind, input$fields), stringsAsFactors = FALSE)
    names(raw) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
      "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")[seq_len(ncol(raw))]
    for (field in intersect(c("chromStart", "chromEnd", "score", "thickStart",
        "thickEnd", "blockCount"), names(raw))) {
      raw[[field]] <- ggchord_format_integer(raw[[field]], field, caller)
    }
    if (any(!nzchar(raw$chrom)) || any(raw$chromEnd <= raw$chromStart)) {
      ggchord_stop(caller, ": nonempty sequence IDs and chromEnd > chromStart are required; ",
        "zero-width insertion records cannot be drawn as intervals")
    }
    if ("score" %in% names(raw) && any(raw$score > 1000)) {
      ggchord_stop(caller, ": BED score must be between 0 and 1000")
    }
    strand <- if ("strand" %in% names(raw)) raw$strand else rep(".", nrow(raw))
    if (any(!strand %in% c("+", "-", "."))) ggchord_stop(caller, ": BED strand must be '+', '-' or '.'")
    for (field in intersect(c("thickStart", "thickEnd"), names(raw))) {
      if (any(raw[[field]] < raw$chromStart | raw[[field]] > raw$chromEnd)) {
        ggchord_stop(caller, ": ", field, " must lie within chromStart/chromEnd")
      }
    }
    if (all(c("thickStart", "thickEnd") %in% names(raw)) && any(raw$thickEnd < raw$thickStart)) {
      ggchord_stop(caller, ": thickEnd must be >= thickStart")
    }
    out <- data.frame(accver = raw$chrom, start = raw$chromStart + 1,
      end = raw$chromEnd, strand = strand,
      directionality = unname(c("+" = "forward", "-" = "reverse", "." = "none")[strand]),
      stringsAsFactors = FALSE)
    if ("name" %in% names(raw)) out$anno <- raw$name
    for (field in setdiff(names(raw), c("chrom", "strand"))) out[[field]] <- raw[[field]]
    if ("blockCount" %in% names(raw) && any(raw$blockCount < 1)) {
      ggchord_stop(caller, ": blockCount must be positive")
    }
    parse_blocks <- function(value, count, name) {
      values <- strsplit(sub(",$", "", value), ",", fixed = TRUE)[[1L]]
      result <- ggchord_format_integer(values, name, caller, if (name == "blockSizes") 1 else 0)
      if (length(result) != count) ggchord_stop(caller, ": ", name, " count differs from blockCount")
      result
    }
    if ("blockSizes" %in% names(raw)) {
      sizes <- lapply(seq_len(nrow(raw)), function(i) parse_blocks(raw$blockSizes[i], raw$blockCount[i], "blockSizes"))
      if ("blockStarts" %in% names(raw)) {
        out$segments <- I(lapply(seq_len(nrow(raw)), function(i) {
          offsets <- parse_blocks(raw$blockStarts[i], raw$blockCount[i], "blockStarts")
          ends <- offsets + sizes[[i]]
          width <- raw$chromEnd[i] - raw$chromStart[i]
          if (offsets[1L] != 0 || tail(ends, 1L) != width || any(ends > width) ||
              (length(offsets) > 1L && any(offsets[-1L] < head(ends, -1L)))) {
            ggchord_stop(caller, ": BED12 blocks must be ordered, non-overlapping, ",
              "within the feature, and span its first/last bases")
          }
          data.frame(start = raw$chromStart[i] + offsets + 1,
            end = raw$chromStart[i] + ends, segment_index = seq_along(offsets),
            segment_type = "standard", stringsAsFactors = FALSE)
        }))
      }
    }
    out$.source_line <- input$line
    ggchord_add_source_file(out, path, source_file)
  })
  ggchord_format_import(ggchord_rbind_fill(tables), paths, "BED", vapply(tables, nrow, integer(1)))
}

#' Read PAF mappings as ggchord links
#'
#' Converts zero-based half-open PAF query/target intervals to one-based
#' inclusive endpoints. Query endpoints increase; target endpoints decrease
#' for relative strand `"-"`. This display convention preserves relative
#' orientation and does not infer a biological strand for either sequence.
#'
#' `length` is the reported alignment block length, including gaps, and
#' `pident` is `100 * nmatch / length` using the reported PAF fields. Both
#' counts may be approximate when the producer did not compute an alignment.
#' The importer does not recompute an alignment, validate a CIGAR, or claim
#' exact base identity. `pident_basis` records this ratio's provenance.
#' Mapping quality 255 means unavailable: `mapq` becomes `NA`, while
#' `source_mapq` retains the original value. Optional typed tags are retained
#' verbatim in the `tags` list-column, including unknown tags and their order.
#'
#' @param file Optional path to one tab-delimited PAF file.
#' @param files Optional paths or wildcard patterns for multiple PAF files.
#' @param source_file Logical. Add a `.source_file` column, default `FALSE`.
#' @return A link data frame with canonical endpoints, `qlen`, `slen`,
#'   `nmatch`, `length`, `pident`, `mapq`, relative orientation, original
#'   zero-based coordinates, tags and `.source_line`. The `ggchord_import`
#'   attribute records input/output coordinate conventions and files.
#' @seealso <https://github.com/lh3/miniasm/blob/master/PAF.md>
#' @md
#' @export
#' @examples
#' paf <- tempfile(fileext = ".paf")
#' writeLines("query\t100\t0\t90\t-\ttarget\t200\t20\t110\t85\t90\t255", paf)
#' links <- read_paf(paf)
#' links[, c("qstart", "qend", "sstart", "send", "mapq")]
#' unlink(paf)
read_paf <- function(file = NULL, files = NULL, source_file = FALSE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  caller <- "read_paf()"
  paths <- ggchord_resolve_files(file, files, "read_paf")
  tables <- lapply(paths, function(path) {
    input <- ggchord_format_lines(path, caller)
    if (any(lengths(input$fields) < 12L)) ggchord_stop(caller, ": PAF requires at least 12 tab-separated fields: ", path)
    raw <- as.data.frame(do.call(rbind, lapply(input$fields, head, 12L)), stringsAsFactors = FALSE)
    names(raw) <- c("qaccver", "qlen", "query_start0", "query_end0", "relative_strand",
      "saccver", "slen", "target_start0", "target_end0", "nmatch", "length", "source_mapq")
    for (field in setdiff(names(raw), c("qaccver", "saccver", "relative_strand"))) {
      raw[[field]] <- ggchord_format_integer(raw[[field]], field, caller)
    }
    if (any(!nzchar(raw$qaccver) | !nzchar(raw$saccver)) ||
        any(!raw$relative_strand %in% c("+", "-"))) {
      ggchord_stop(caller, ": PAF requires nonempty sequence IDs and relative strand '+' or '-' ")
    }
    if (any(raw$qlen < 1 | raw$slen < 1 | raw$query_end0 <= raw$query_start0 |
        raw$target_end0 <= raw$target_start0 | raw$query_end0 > raw$qlen |
        raw$target_end0 > raw$slen)) {
      ggchord_stop(caller, ": PAF intervals must be nonempty and within their declared sequence lengths")
    }
    if (any(raw$length < 1 | raw$nmatch > raw$length | raw$source_mapq > 255)) {
      ggchord_stop(caller, ": PAF requires 0 <= nmatch <= positive block length and mapq in [0, 255]")
    }
    reverse <- raw$relative_strand == "-"
    out <- data.frame(qaccver = raw$qaccver, saccver = raw$saccver,
      qstart = raw$query_start0 + 1, qend = raw$query_end0,
      sstart = ifelse(reverse, raw$target_end0, raw$target_start0 + 1),
      send = ifelse(reverse, raw$target_start0 + 1, raw$target_end0),
      stringsAsFactors = FALSE)
    for (field in setdiff(names(raw), c("qaccver", "saccver"))) out[[field]] <- raw[[field]]
    out$pident <- 100 * out$nmatch / out$length
    out$pident_basis <- "reported_matches/block_length"
    out$mapq <- ifelse(out$source_mapq == 255, NA_real_, out$source_mapq)
    out$orientation_known <- TRUE
    out$direction <- ifelse(reverse, "reverse", "same")
    out$tags <- I(lapply(input$fields, function(fields) fields[-seq_len(12L)]))
    out$.source_line <- input$line
    ggchord_add_source_file(out, path, source_file)
  })
  ggchord_format_import(ggchord_rbind_fill(tables), paths, "PAF", vapply(tables, nrow, integer(1)))
}
