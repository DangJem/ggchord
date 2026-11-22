# Explicit zero-based file adapters. Specifications (not copied parser code):
# https://genome.ucsc.edu/FAQ/FAQformat.html#format1
# https://github.com/lh3/miniasm/blob/master/PAF.md

ggchord_format_header <- function(lines) {
  text <- trimws(lines)
  header <- grepl("^track([[:space:]]+[^[:space:]=]+=|$)", text) |
    grepl("^browser([[:space:]]+(position|hide|dense|squish|pack|full)([[:space:]]|$)|$)", text)
  # A tab-separated interval wins over a header-like accession prefix.
  for (i in which(header & grepl("\t", lines, fixed = TRUE))) {
    fields <- strsplit(lines[i], "\t", fixed = TRUE)[[1L]]
    if (length(fields) >= 3L && all(grepl("^[0-9]+$", fields[2:3]))) header[i] <- FALSE
  }
  header
}

ggchord_format_lines <- function(file, caller, bed = FALSE) {
  lines <- ggchord_read_text(file, caller)
  keep <- nzchar(trimws(lines)) & !grepl("^[[:space:]]*#", lines)
  header <- if (bed) ggchord_format_header(lines) else rep(FALSE, length(lines))
  if (bed) keep <- keep & !header
  index <- which(keep)
  if (!length(index)) ggchord_stop(caller, ": no data rows in ", file)
  values <- lapply(lines[index], function(line) {
    if (grepl("\t", line, fixed = TRUE)) strsplit(line, "\t", fixed = TRUE)[[1L]] else {
      if (!bed) ggchord_stop(caller, ": PAF fields must be tab-separated in ", file)
      strsplit(trimws(line), "[[:space:]]+")[[1L]]
    }
  })
  metadata <- which(header | grepl("^[[:space:]]*#", lines))
  list(fields = values, line = index,
    headers = data.frame(.source_line = metadata, text = lines[metadata], stringsAsFactors = FALSE))
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

ggchord_canonical_read_links <- function(parts, caller, directional_order = FALSE) {
  links <- ggchord_rbind_fill(parts) %||% data.frame(qaccver = character(),
    saccver = character(), qstart = numeric(), qend = numeric(), sstart = numeric(), send = numeric())
  ggchord_order_link_endpoints(links, caller, directional_order)
}

#' Read BED features and their observed segments
#'
#' Read standard BED3 through BED12. Convert each nonempty half-open range
#' to a one-based closed interval. BED12 blocks are preserved as segments.
#' A BED name provides a label, and a missing strand stays unknown.
#'
#' @param file Local BED paths.
#' @param source_file Whether to include row-level source paths.
#' @return A feature data frame and file provenance available from
#'   [data_report()].
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' f <- tempfile(fileext = ".bed")
#' writeLines("A	99	240	region1	0	+", f)
#' bed <- read_bed(f)
#' ggchord(seqs) + geom_feature(data = bed)
#' unlink(f)
read_bed <- function(file, source_file = FALSE) {
  ggchord_local_console_errors()
  caller <- "read_bed()"
  paths <- ggchord_resolve_files(file, NULL, "read_bed")
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
    if ("name" %in% names(raw)) out$label <- raw$name
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

#' Read PAF relationships
#'
#' Preserve query/subject definitions, declared lengths, directional evidence,
#' and optional tags. Numeric position fields are converted from zero-based
#' half-open ranges to one-based closed intervals. Identity is labelled as
#' reported match count divided by alignment block length.
#'
#' @param file Local PAF paths.
#' @param source_file Whether to include row-level source paths.
#' @return A relationship table with ordered endpoints, declared reference
#'   lengths, matching positions, alignment block length, mapq and tags.
#' @export
#' @examples
#' f <- tempfile(fileext=".paf")
#' writeLines("A	1000	99	240	+	B	800	149	290	138	141	60", f)
#' paf <- read_paf(f)
#' unlink(f)
read_paf <- function(file, source_file = FALSE) {
  ggchord_local_console_errors()
  caller <- "read_paf()"
  paths <- ggchord_resolve_files(file, NULL, "read_paf")
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
    out$match_count <- out$nmatch
    out$alignment_block_length <- out$length
    out$qstrand <- "+"
    out$sstrand <- raw$relative_strand
    out <- ggchord_order_link_endpoints(out, caller)
    out$pident <- 100 * out$nmatch / out$length
    out$pident_basis <- "reported_matches/block_length"
    out$mapq <- ifelse(out$source_mapq == 255, NA_real_, out$source_mapq)
    out$orientation_known <- TRUE
    out$direction <- ifelse(reverse, "reverse", "same")
    out$tags <- I(lapply(input$fields, function(fields) fields[-seq_len(12L)]))
    out$.source_line <- input$line
    ggchord_add_source_file(out, path, source_file)
  })
  out <- ggchord_format_import(ggchord_rbind_fill(tables), paths, "PAF", vapply(tables, nrow, integer(1)))
  attr(out, "ggchord_import")$sequences <- ggchord_alignment_dictionary(list(
    data.frame(accver = out$qaccver, length = out$qlen),
    data.frame(accver = out$saccver, length = out$slen)), "bp", caller)
  out
}
