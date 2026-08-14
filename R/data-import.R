# data-import.R - lightweight import helpers for common bioinformatics formats
#
# v0.8.0: read_blast(), read_gff3() and read_fasta_lengths() turn common file
# formats into the tidy data frames ggchord() expects, with no tidyverse
# dependency.  External command-line preparation (BLAST, seqkit, ...) stays
# outside the package; these helpers cover the "read the file in R" step.

#' Read BLAST tabular output into ribbon_data format
#'
#' Parses BLAST \code{-outfmt 6} or \code{-outfmt 7} tabular output into a data
#' frame that can be passed directly to \code{\link{ggchord}()} as
#' \code{ribbon_data}. The standard 12-column layout
#' (\code{qaccver saccver pident length mismatch gapopen qstart qend sstart
#' send evalue bitscore}) and the 17-column extension used by the README
#' workflow (\code{...+ qcovs qlen slen sstrand stitle}) are recognised
#' automatically. Optional columns such as \code{evalue}, \code{bitscore},
#' \code{qcovs}, \code{qlen}, \code{slen}, \code{sstrand} and \code{stitle} are
#' preserved for filtering and hover diagnostics.
#'
#' @param file Path to a BLAST tabular output file.
#' @param format Character. \code{"auto"} (default) detects the column layout
#'   from the number of columns; \code{"outfmt6"} / \code{"outfmt7"} require the
#'   standard 12/17-column layouts; \code{"custom"} requires \code{col_names}.
#' @param col_names Optional character vector naming the columns in the file,
#'   used with \code{format = "custom"} or to override auto-detection.
#' @param comment Character comment character, default \code{"#"} (BLAST
#'   outfmt 7 header lines start with \code{#}).
#' @param ... Additional arguments passed to \code{\link[utils]{read.table}}
#'   (e.g. \code{na.strings}).
#'
#' @return A data.frame with the required ribbon columns first
#'   (\code{qaccver}, \code{saccver}, \code{length}, \code{pident},
#'   \code{qstart}, \code{qend}, \code{sstart}, \code{send}) followed by any
#'   preserved optional columns.
#' @export
#'
#' @examples
#' \donttest{
#' blast_file <- tempfile(fileext = ".o7")
#' writeLines(c(
#'   "# BLASTN 2.13.0+",
#'   "# Query: seqA",
#'   "seqA\tseqB\t98.5\t1200\t18\t0\t1\t1200\t1\t1200\t1e-180\t2400\t100\t5000\t4800\tplus\tseqB",
#'   "seqA\tseqC\t95.0\t800\t40\t2\t1300\t2100\t50\t850\t1e-100\t1500\t80\t5000\t6000\tminus\tseqC"
#' ), blast_file)
#' rb <- read_blast(blast_file)
#' head(rb)
#' }
read_blast <- function(file, format = c("auto", "outfmt6", "outfmt7", "custom"),
                       col_names = NULL, comment = "#", ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  format <- match.arg(format)
  if (is.character(file) && length(file) == 1 && !file.exists(file)) {
    ggchord_stop("read_blast(): file not found: ", file, call. = FALSE)
  }

  raw <- utils::read.table(
    file, sep = "\t", header = FALSE, quote = "",
    comment.char = comment, fill = TRUE,
    stringsAsFactors = FALSE, ...
  )
  if (ncol(raw) == 0 || nrow(raw) == 0) {
    ggchord_stop("read_blast(): no data rows parsed from ", file, call. = FALSE)
  }

  std12 <- c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen",
             "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  ext17 <- c(std12, "qcovs", "qlen", "slen", "sstrand", "stitle")

  if (is.null(col_names)) {
    col_names <- if (format == "custom") {
      ggchord_stop("read_blast(): col_names is required when format = 'custom'",
           call. = FALSE)
    } else if (ncol(raw) == 12) {
      std12
    } else if (ncol(raw) == 17) {
      ext17
    } else {
      ggchord_stop(sprintf(paste0(
        "read_blast(): cannot auto-detect the column layout (%d columns found). ",
        "Expected 12 (outfmt 6) or 17 (outfmt 7 with qcovs/qlen/slen/sstrand/stitle). ",
        "Pass col_names and format = 'custom' for other layouts."),
        ncol(raw)), call. = FALSE)
    }
  } else {
    if (length(col_names) != ncol(raw)) {
      ggchord_stop(sprintf("read_blast(): col_names has %d entries but the file has %d columns",
                   length(col_names), ncol(raw)), call. = FALSE)
    }
  }
  names(raw) <- col_names

  req <- c("qaccver", "saccver", "length", "pident",
           "qstart", "qend", "sstart", "send")
  missing <- setdiff(req, col_names)
  if (length(missing) > 0) {
    ggchord_stop("read_blast(): required ribbon_data column(s) missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  num_cols <- intersect(c("pident", "length", "mismatch", "gapopen",
                          "qstart", "qend", "sstart", "send",
                          "evalue", "bitscore", "qcovs", "qlen", "slen"),
                        col_names)
  for (nm in num_cols) {
    if (!is.numeric(raw[[nm]])) {
      converted <- suppressWarnings(as.numeric(raw[[nm]]))
      if (anyNA(converted) && !anyNA(raw[[nm]])) {
        warning(sprintf("read_blast(): column '%s' contains non-numeric values that were coerced to NA",
                        nm), call. = FALSE)
      }
      raw[[nm]] <- converted
    }
  }

  extra <- setdiff(col_names, req)
  raw[, c(req, extra), drop = FALSE]
}

#' Read a GFF3 file into gene_data format
#'
#' Parses a GFF3 file (9 tab-separated columns) and returns a data frame that
#' can be passed to \code{\link{ggchord}()} as \code{gene_data}. Only features
#' whose \code{type} is in \code{feature_types} are kept. The annotation label
#' (\code{anno}) is extracted from the GFF3 attributes column by trying the
#' keys in \code{anno_from} in order (e.g. \code{product}, then \code{Name}).
#' The original useful columns (\code{type}, \code{source}, \code{score},
#' \code{phase}, \code{attributes}) are preserved.
#'
#' @param file Path to a GFF3 file.
#' @param feature_types Character vector of feature types to keep, default
#'   \code{"CDS"}. Common choices: \code{"gene"}, \code{"tRNA"},
#'   \code{"rRNA"}, \code{"repeat_region"}.
#' @param anno_from Character vector of GFF3 attribute keys, tried in order to
#'   fill the \code{anno} column, default \code{c("product", "Name", "gene",
#'   "ID")}. \code{anno} is \code{NA} when none of the keys is present.
#' @param unstranded Character, default \code{"plus"}. How to treat features
#'   whose strand is \code{"."} or \code{"?"}: \code{"plus"} maps them to
#'   \code{"+"}, \code{"drop"} removes them (ggchord requires \code{"+"} or
#'   \code{"-"}).
#'
#' @return A data.frame with \code{seq_id}, \code{start}, \code{end},
#'   \code{strand}, \code{anno} followed by \code{type}, \code{source},
#'   \code{score}, \code{phase} and \code{attributes}.
#' @export
#'
#' @examples
#' \donttest{
#' gff <- tempfile(fileext = ".gff3")
#' writeLines(c(
#'   "##gff-version 3",
#'   "seqA\tsource\tCDS\t101\t500\t.\t+\t0\tID=cds1;product=hypothetical protein",
#'   "seqA\tsource\tCDS\t600\t900\t.\t-\t0\tID=cds2;Name=integrase",
#'   "seqA\tsource\ttRNA\t1000\t1080\t.\t+\t0\tID=trna1"
#' ), gff)
#' gd <- read_gff3(gff)
#' gd
#' }
read_gff3 <- function(file, feature_types = "CDS",
                      anno_from = c("product", "Name", "gene", "ID"),
                      unstranded = c("plus", "drop")) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  unstranded <- match.arg(unstranded)
  if (!is.character(feature_types) || length(feature_types) == 0) {
    ggchord_stop("read_gff3(): feature_types must be a non-empty character vector",
         call. = FALSE)
  }
  if (is.character(file) && length(file) == 1 && !file.exists(file)) {
    ggchord_stop("read_gff3(): file not found: ", file, call. = FALSE)
  }

  raw <- utils::read.table(file, sep = "\t", header = FALSE, quote = "",
                           comment.char = "#", fill = TRUE,
                           stringsAsFactors = FALSE)
  if (ncol(raw) < 9) {
    ggchord_stop(sprintf("read_gff3(): expected 9 GFF3 columns but found %d in %s",
                 ncol(raw), file), call. = FALSE)
  }
  names(raw)[1:9] <- c("seqid", "source", "type", "start", "end",
                       "score", "strand", "phase", "attributes")

  raw$start <- suppressWarnings(as.numeric(raw$start))
  raw$end <- suppressWarnings(as.numeric(raw$end))
  if (anyNA(raw$start) || anyNA(raw$end)) {
    warning("read_gff3(): some start/end values could not be parsed as numbers",
            call. = FALSE)
  }

  sub <- raw[raw$type %in% feature_types, , drop = FALSE]
  if (nrow(sub) == 0) {
    ggchord_stop(sprintf(paste0("read_gff3(): no features of type(s) %s found in %s; ",
                        "use feature_types to select the features to keep"),
                 paste(feature_types, collapse = ", "), file), call. = FALSE)
  }

  if (unstranded == "plus") {
    sub$strand[!sub$strand %in% c("+", "-")] <- "+"
  } else {
    sub <- sub[sub$strand %in% c("+", "-"), , drop = FALSE]
  }
  if (nrow(sub) == 0) {
    ggchord_stop("read_gff3(): no features with a '+' or '-' strand remain after filtering",
         call. = FALSE)
  }

  anno <- vapply(sub$attributes, function(a) {
    extract_gff3_attr(a, anno_from)
  }, character(1), USE.NAMES = FALSE)

  data.frame(
    seq_id = sub$seqid,
    start = sub$start,
    end = sub$end,
    strand = sub$strand,
    anno = anno,
    type = sub$type,
    source = sub$source,
    score = sub$score,
    phase = sub$phase,
    attributes = sub$attributes,
    stringsAsFactors = FALSE
  )
}

#' Extract a GFF3 attribute value by key
#' @keywords internal
extract_gff3_attr <- function(attrs, keys) {
  if (is.na(attrs) || !nzchar(attrs)) return(NA_character_)
  parts <- strsplit(attrs, ";", fixed = TRUE)[[1]]
  for (k in keys) {
    pat <- paste0("^", k, "=")
    hit <- parts[grepl(pat, parts)]
    if (length(hit) > 0) {
      v <- sub(pat, "", hit[1])
      if (nzchar(v)) return(gff3_percent_decode(v))
    }
  }
  NA_character_
}

#' Decode GFF3 percent-encoding (\%XX) without touching literal '+'
#' @keywords internal
gff3_percent_decode <- function(x) {
  if (!grepl("%", x, fixed = TRUE)) return(x)
  vapply(x, function(s) {
    utils::URLdecode(gsub("%(?![0-9A-Fa-f]{2})", "%25", s, perl = TRUE))
  }, character(1), USE.NAMES = FALSE)
}

#' Read FASTA headers and sequence lengths
#'
#' Reads a FASTA file and returns a \code{seq_id}/\code{length} data frame
#' suitable for \code{\link{ggchord}()}. The sequence ID is the first
#' whitespace-delimited token of each header by default; pass
#' \code{header_delim} (e.g. \code{"|"}) to split NCBI-style headers further
#' and keep only the first field.
#'
#' @param file Path to a FASTA file.
#' @param header_delim Optional character. When given, each header is split at
#'   every occurrence of this delimiter and only the first piece is kept (e.g.
#'   \code{"|"} for NCBI headers such as \code{>NC_000001.1|cds|...}).
#'
#' @return A data.frame with columns \code{seq_id} and \code{length}. A warning
#'   is emitted when the file contains duplicate sequence IDs.
#' @export
#'
#' @examples
#' fasta <- tempfile(fileext = ".fna")
#' writeLines(c(
#'   ">seqA some description",
#'   "ACGTACGTACGTACGT",
#'   "ACGTACGT",
#'   ">seqB",
#'   "TTTTGGGG"
#' ), fasta)
#' read_fasta_lengths(fasta)
read_fasta_lengths <- function(file, header_delim = NULL) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (is.character(file) && length(file) == 1 && !file.exists(file)) {
    ggchord_stop("read_fasta_lengths(): file not found: ", file, call. = FALSE)
  }
  lines <- readLines(file, warn = FALSE)
  if (length(lines) == 0) {
    ggchord_stop("read_fasta_lengths(): file is empty", call. = FALSE)
  }
  is_header <- grepl("^>", lines)
  if (!any(is_header)) {
    ggchord_stop("read_fasta_lengths(): no FASTA headers (lines starting with '>') found",
         call. = FALSE)
  }
  ids <- sub("^>", "", lines[is_header])
  ids <- trimws(ids)
  if (!is.null(header_delim)) {
    ids <- vapply(strsplit(ids, header_delim, fixed = TRUE),
                  function(x) x[1], character(1))
  } else {
    ids <- vapply(strsplit(ids, "[[:space:]]+"),
                  function(x) x[1], character(1))
  }

  header_pos <- which(is_header)
  starts <- header_pos + 1L
  ends <- c(header_pos[-1] - 1L, length(lines))
  lens <- vapply(seq_along(header_pos), function(i) {
    if (starts[i] > ends[i]) return(0L)
    sum(nchar(gsub("[[:space:]]+", "", lines[starts[i]:ends[i]])))
  }, integer(1))

  out <- data.frame(seq_id = ids, length = as.numeric(lens),
                    stringsAsFactors = FALSE)
  dup <- unique(ids[duplicated(ids)])
  if (length(dup) > 0) {
    warning("read_fasta_lengths(): duplicate sequence IDs found: ",
            paste(dup, collapse = ", "), call. = FALSE)
  }
  out
}
