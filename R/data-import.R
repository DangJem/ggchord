# data-import.R - lightweight import helpers for common bioinformatics formats
#
# v0.8.0: read_blast(), read_gff3() and read_fasta_lengths() turn common file
# formats into the tidy data frames ggchord() expects, with no tidyverse
# dependency.  External command-line preparation (BLAST, seqkit, ...) stays
# outside the package; these helpers cover the "read the file in R" step.

# Resolve `file` and `files` into a character vector of paths. `files` may
# contain wildcard patterns (e.g. "examples/fasta/*.fna") and/or literal paths.
ggchord_resolve_files <- function(file, files, fun) {
  for (name in c("file", "files")) {
    value <- get(name)
    if (!is.null(value) && (!is.character(value) || anyNA(value) || any(!nzchar(value))))
      ggchord_stop(fun, "(): `", name, "` must contain non-missing, nonempty character paths")
  }
  paths <- c(file, files)
  if (!length(paths)) ggchord_stop(fun, "(): provide `file` or `files`")

  expanded <- unlist(lapply(paths, function(x) {
    if (grepl("[*?[]", x)) {
      g <- Sys.glob(x)
      if (length(g) == 0) {
        ggchord_stop(sprintf("%s(): no files match pattern: %s", fun, x),
             call. = FALSE)
      }
      g
    } else {
      if (!file.exists(x) || dir.exists(x)) {
        ggchord_stop(sprintf("%s(): file not found: %s", fun, x),
             call. = FALSE)
      }
      x
    }
  }), use.names = FALSE)

  if (any(dir.exists(expanded))) ggchord_stop(fun, "(): expected files, not directories")
  if (length(expanded) == 0) {
    ggchord_stop(sprintf("%s(): no files to read", fun), call. = FALSE)
  }
  expanded[!duplicated(normalizePath(expanded, mustWork = TRUE))]
}

# All text readers use UTF-8; reject invalid bytes instead of rewriting IDs.
ggchord_read_text <- function(file, caller) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  invalid <- which(is.na(iconv(lines, from = "UTF-8", to = "UTF-8", sub = NA)))
  if (length(invalid)) ggchord_stop(caller, ": invalid UTF-8 in ", file,
    " at line ", invalid[1L], "; convert the file to UTF-8 explicitly")
  if (length(lines)) lines[1L] <- sub("^\ufeff", "", lines[1L])
  lines
}

# Combine data.frames that may have different columns (e.g. BLAST outfmt 6 and
# outfmt 7 files), filling missing columns with typed NA values.
ggchord_rbind_fill <- function(xs) {
  xs <- Filter(function(x) !is.null(x), xs)
  if (length(xs) == 0) return(NULL)
  if (length(xs) == 1) return(xs[[1]])
  template <- xs
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(d) {
    miss <- setdiff(cols, names(d))
    if (length(miss) > 0) {
      for (m in miss) {
        src <- Find(function(x) m %in% names(x), template)
        d[[m]] <- rep(if (!is.null(src)) src[[m]][NA_integer_] else NA,
          length.out = nrow(d))
      }
    }
    d[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
}

# Translate labels emitted by BLAST outfmt 7's `# Fields:` header to the
# conventional ggchord/BLAST column names. Unknown fields are retained using a
# stable snake-case name rather than being discarded.
ggchord_blast_field_names <- function(file, lines = NULL) {
  if (is.null(lines)) lines <- ggchord_read_text(file, "read_blast()")
  headers <- grep("^# Fields:", lines, value = TRUE)
  if (length(headers) == 0) return(NULL)
  parsed <- lapply(headers, function(x) {
    trimws(strsplit(sub("^# Fields:\\s*", "", x), ",", fixed = TRUE)[[1]])
  })
  signatures <- vapply(parsed, paste, collapse = "\r", character(1))
  if (length(unique(signatures)) != 1) {
    ggchord_stop("read_blast(): inconsistent '# Fields:' headers in ", file)
  }
  fields <- tolower(parsed[[1]])
  lookup <- c(
    "query id" = "qaccver", "query acc." = "qaccver",
    "query acc.ver" = "qaccver", "subject id" = "saccver",
    "subject acc." = "saccver", "subject acc.ver" = "saccver",
    "% identity" = "pident", "alignment length" = "length",
    "mismatches" = "mismatch", "gap opens" = "gapopen",
    "q. start" = "qstart", "q. end" = "qend",
    "s. start" = "sstart", "s. end" = "send",
    "evalue" = "evalue", "bit score" = "bitscore",
    "% positives" = "ppos", "query coverage per subject" = "qcovs",
    "query length" = "qlen", "subject length" = "slen",
    "subject strand" = "sstrand", "subject title" = "stitle"
  )
  out <- unname(lookup[fields])
  unknown <- is.na(out)
  out[unknown] <- gsub("(^_+|_+$)", "",
                       gsub("[^a-z0-9]+", "_", fields[unknown]))
  make.unique(out, sep = "_")
}

ggchord_add_source_file <- function(data, file, source_file) {
  if (!is.logical(source_file) || length(source_file) != 1 || is.na(source_file)) {
    ggchord_stop("source_file must be TRUE or FALSE")
  }
  if (isTRUE(source_file)) data$.source_file <- as.character(file)
  data
}


# Internal single-file reader; the public read_blast() wrapper resolves
# `file` / `files` and combines the result.
read_blast_single <- function(file, format = c("auto", "outfmt6", "outfmt7", "custom"),
                               col_names = NULL, comment = "#", ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  format <- match.arg(format)
  if (is.character(file) && length(file) == 1 && !file.exists(file)) {
    ggchord_stop("read_blast(): file not found: ", file, call. = FALSE)
  }

  lines <- ggchord_read_text(file, "read_blast()")
  header_names <- ggchord_blast_field_names(file, lines)
  dots <- list(...)
  if (!is.character(comment) || length(comment) != 1L || is.na(comment) || nchar(comment) > 1L)
    ggchord_stop("read_blast(): comment must be one character or an empty string")
  content <- if (nzchar(comment)) vapply(strsplit(paste0(lines, "\t"), comment, fixed = TRUE),
    `[`, character(1), 1L) else lines
  source_line <- which(nzchar(trimws(content)))
  if (!is.null(dots$skip)) source_line <- source_line[source_line > dots$skip]
  if (!is.null(dots$nrows) && dots$nrows >= 0) source_line <- head(source_line, dots$nrows)
  if (!length(source_line)) ggchord_stop("read_blast(): no data rows parsed from ", file)
  widths <- lengths(strsplit(content[source_line], "\t", fixed = TRUE))
  if (length(unique(widths)) != 1L) ggchord_stop(
    "read_blast(): inconsistent field counts in ", file, " at source lines ",
    paste(source_line[widths != widths[1L]], collapse = ", "))
  classes <- dots$colClasses %||% "character"
  dots$colClasses <- NULL
  missing_values <- dots$na.strings %||% character()
  dots$na.strings <- NULL
  raw <- do.call(utils::read.table, c(list(
    text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE, quote = "",
    comment.char = comment, fill = FALSE,
    stringsAsFactors = FALSE, colClasses = classes, na.strings = missing_values), dots))
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
    } else if (format == "outfmt7") {
      if (is.null(header_names)) {
        ggchord_stop("read_blast(): format = 'outfmt7' requires a '# Fields:' header")
      }
      header_names
    } else if (format == "outfmt6") {
      if (ncol(raw) != 12) {
        ggchord_stop("read_blast(): format = 'outfmt6' requires the standard 12 columns; use format = 'custom' with col_names for another layout")
      }
      std12
    } else if (!is.null(header_names)) {
      header_names
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
  if (length(col_names) != ncol(raw)) {
    ggchord_stop(sprintf(
      "read_blast(): the declared field list has %d entries but data rows have %d columns",
      length(col_names), ncol(raw)))
  }
  if (!is.character(col_names) || anyNA(col_names) || any(!nzchar(col_names)) ||
      anyDuplicated(col_names) || ".source_line" %in% col_names)
    ggchord_stop("read_blast(): col_names must be unique nonempty names; .source_line is reserved for file provenance")
  names(raw) <- col_names
  raw$.source_line <- source_line[seq_len(nrow(raw))]

  req <- ggchord_ribbon_required_columns()
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
      value <- if (is.factor(raw[[nm]])) as.character(raw[[nm]]) else raw[[nm]]
      converted <- suppressWarnings(as.numeric(value))
      if (any(is.na(converted) & !is.na(raw[[nm]]))) {
        warning(sprintf("read_blast(): column '%s' contains non-numeric values that were coerced to NA",
                        nm), call. = FALSE)
      }
      raw[[nm]] <- converted
    }
  }
  for (nm in setdiff(col_names, c(num_cols, "qaccver", "saccver"))) {
    if (is.character(raw[[nm]])) raw[[nm]] <- utils::type.convert(raw[[nm]], as.is = TRUE,
      na.strings = character())
  }

  extra <- setdiff(names(raw), req)
  raw[, c(req, extra), drop = FALSE]
}

#' Read one or more BLAST tabular output files into ribbon_data format
#'
#' `file` reads a single BLAST tabular output file; `files` reads multiple
#' files at once and combines them. `files` accepts a character vector of
#' literal paths and/or wildcard patterns (e.g. `"examples/blastn/*.o7"`).
#'
#' @param file Optional path vector or glob patterns; read before `files`.
#' @param files Optional character vector of BLAST tabular output files
#'   (literal paths and/or wildcard patterns). All matched files are read and
#'   combined into one data.frame.
#' @param format Character. `"auto"` (default) detects the column layout
#'   from the file; `"outfmt6"` requires the standard 12-column layout,
#'   `"outfmt7"` parses its `# Fields:` declaration, and `"custom"` requires
#'   `col_names`.
#' @param col_names Optional character vector naming the columns in the file,
#'   used with `format = "custom"` or to override auto-detection.
#' @param comment Character comment character, default `"#"` (BLAST
#'   outfmt 7 header lines start with `#`).
#' @param source_file Logical. Add a `.source_file` column when reading one or
#'   more files, default `FALSE`.
#' @param ... Additional arguments passed to [utils::read.table()]
#'   (e.g. `na.strings`).
#'
#' @return A data.frame with the required ribbon columns first
#'   (`qaccver`, `saccver`, `qstart`, `qend`, `sstart`, `send`) followed by any preserved optional columns.
#'   A `ggchord_import` attribute records files, per-file output row counts and
#'   the one-based inclusive coordinate convention.
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#'
#' @examples
#' library(ggchord)
#' blast_file <- tempfile(fileext = ".o7")
#' writeLines(c(
#'   "# BLASTN 2.13.0+",
#'   "seqA	seqB	98.5	1200	18	0	1	1200	1	1200	1e-180	2400"
#' ), blast_file)
#' read_blast(blast_file)
read_blast <- function(file = NULL, files = NULL,
                       format = c("auto", "outfmt6", "outfmt7", "custom"),
                       col_names = NULL, comment = "#", source_file = FALSE,
                       ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  format <- match.arg(format)
  paths <- ggchord_resolve_files(file, files, "read_blast")
  out <- lapply(paths, function(f) {
    data <- read_blast_single(f, format = format, col_names = col_names,
                              comment = comment, ...)
    ggchord_add_source_file(data, f, source_file)
  })
  ggchord_format_import(ggchord_rbind_fill(out), paths, "BLAST",
    vapply(out, nrow, integer(1)), "1-based-inclusive")
}

# Internal single-file reader; the public read_gff3() wrapper resolves
# `file` / `files` and combines the result.
read_gff3_single <- function(file, feature_types = "CDS",
                              anno_from = c("product", "Name", "gene", "ID"),
                              unstranded = c("plus", "drop", "keep")) {
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

  lines <- ggchord_read_text(file, "read_gff3()")
  fasta <- which(lines == "##FASTA")
  if (length(fasta) > 0) {
    lines <- if (fasta[1] > 1) lines[seq_len(fasta[1] - 1L)] else character(0)
  }
  source_line <- which(nzchar(trimws(lines)) & !grepl("^#", lines))
  data_lines <- lines[source_line]
  if (length(data_lines) == 0) {
    ggchord_stop("read_gff3(): no feature rows found in ", file)
  }
  if (any(lengths(strsplit(data_lines, "\t", fixed = TRUE)) != 9L))
    ggchord_stop("read_gff3(): expected exactly 9 tab-separated fields in ", file)
  raw <- utils::read.table(text = paste(data_lines, collapse = "\n"),
                           sep = "\t", header = FALSE, quote = "",
                           comment.char = "", fill = TRUE,
                           stringsAsFactors = FALSE, colClasses = "character",
                           na.strings = character())
  raw$.source_line <- source_line
  if (ncol(raw) < 9) {
    ggchord_stop(sprintf("read_gff3(): expected 9 GFF3 columns but found %d in %s",
                 ncol(raw), file), call. = FALSE)
  }
  names(raw)[1:9] <- c("seqid", "source", "type", "start", "end",
                       "score", "strand", "phase", "attributes")

  raw$start <- suppressWarnings(as.numeric(raw$start))
  raw$end <- suppressWarnings(as.numeric(raw$end))
  for (column in c("score", "phase")) raw[[column]] <- utils::type.convert(
    raw[[column]], as.is = TRUE, na.strings = character())
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

  source_strand <- sub$strand
  if (unstranded == "plus") {
    sub$strand[!sub$strand %in% c("+", "-")] <- "+"
  } else if (unstranded == "drop") {
    keep <- sub$strand %in% c("+", "-")
    sub <- sub[keep, , drop = FALSE]
    source_strand <- source_strand[keep]
  } else {
    sub$strand[!sub$strand %in% c("+", "-")] <- "."
  }
  if (nrow(sub) == 0) {
    ggchord_stop("read_gff3(): no features with a '+' or '-' strand remain after filtering",
         call. = FALSE)
  }

  anno <- vapply(sub$attributes, function(a) {
    extract_gff3_attr(a, anno_from)
  }, character(1), USE.NAMES = FALSE)

  out <- data.frame(
    accver = sub$seqid,
    start = sub$start,
    end = sub$end,
    strand = sub$strand,
    anno = anno,
    type = sub$type,
    source = sub$source,
    score = sub$score,
    phase = sub$phase,
    attributes = sub$attributes,
    .source_line = sub$.source_line,
    stringsAsFactors = FALSE
  )
  out$source_strand <- source_strand
  out$feature_id <- vapply(sub$attributes, extract_gff3_attr, character(1), keys = "ID")
  # Split before decoding: an escaped comma is part of one identifier.
  out$parent_ids <- I(lapply(sub$attributes, function(a) {
    parts <- strsplit(a, ";", fixed = TRUE)[[1]]
    hit <- parts[startsWith(parts, "Parent=")]
    if (!length(hit) || !nzchar(substring(hit[1L], 8L))) return(character())
    vapply(strsplit(substring(hit[1L], 8L), ",", fixed = TRUE)[[1]],
      gff3_percent_decode, character(1), USE.NAMES = FALSE)
  }))
  out
}

#' Read one or more GFF3 files into gene_data format
#'
#' `file` reads a single GFF3 file; `files` reads multiple files at once and
#' combines them. `files` accepts a character vector of literal paths and/or
#' wildcard patterns (e.g. `"examples/gff3/*.gff3"`).
#'
#' @param file Optional path vector or glob patterns; read before `files`.
#' @param files Optional character vector of GFF3 files (literal paths and/or
#'   wildcard patterns). All matched files are read and combined.
#' @param feature_types Character vector of feature types to keep, default
#'   `"CDS"`.
#' @param anno_from Character vector of GFF3 attribute keys, tried in order to
#'   fill the `anno` column.
#' @param unstranded `"plus"` (legacy default) assigns unknown strands to `+`;
#'   `"drop"` removes them. Use `"keep"` with [geom_feature()] to retain them
#'   as nondirectional `.` intervals. `source_strand` always records the input.
#' @param source_file Logical. Add a `.source_file` column when reading one or
#'   more files, default `FALSE`.
#'
#' @return A data.frame with `accver`, `start`, `end`, `strand`, `anno`
#'   followed by `type`, `source`, `score`, `phase`, `attributes`,
#'   `source_strand`, `feature_id` (decoded GFF3 ID, or missing) and a
#'   `parent_ids` list-column. Repeated IDs for discontinuous features and
#'   multiple parents are retained; selecting feature types may omit parents.
#'   A `ggchord_import` attribute records files, per-file output row counts and
#'   the one-based inclusive coordinate convention.
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#'
#' @examples
#' library(ggchord)
#' gff <- tempfile(fileext = ".gff3")
#' writeLines(c(
#'   "##gff-version 3",
#'   paste0("seqA	source	CDS	101	500	.	+	0	",
#'          "ID=cds1;product=hypothetical protein")
#' ), gff)
#' read_gff3(gff)
#' read_gff3(gff, unstranded = "keep")
read_gff3 <- function(file = NULL, files = NULL, feature_types = "CDS",
                      anno_from = c("product", "Name", "gene", "ID"),
                      unstranded = c("plus", "drop", "keep"), source_file = FALSE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  unstranded <- match.arg(unstranded)
  paths <- ggchord_resolve_files(file, files, "read_gff3")
  out <- lapply(paths, function(f) {
    data <- read_gff3_single(f, feature_types = feature_types,
                             anno_from = anno_from, unstranded = unstranded)
    ggchord_add_source_file(data, f, source_file)
  })
  ggchord_format_import(ggchord_rbind_fill(out), paths, "GFF3",
    vapply(out, nrow, integer(1)), "1-based-inclusive")
}

#' Extract a GFF3 attribute value by key
#' @noRd
extract_gff3_attr <- function(attrs, keys) {
  if (is.na(attrs) || !nzchar(attrs)) return(NA_character_)
  parts <- strsplit(attrs, ";", fixed = TRUE)[[1]]
  for (k in keys) {
    prefix <- paste0(k, "=")
    hit <- parts[startsWith(parts, prefix)]
    if (length(hit) > 0) {
      v <- substring(hit[1], nchar(prefix) + 1L)
      if (nzchar(v)) return(gff3_percent_decode(v))
    }
  }
  NA_character_
}

#' Decode GFF3 percent-encoding (\%XX) without touching literal '+'
#' @noRd
gff3_percent_decode <- function(x) {
  if (!grepl("%", x, fixed = TRUE)) return(x)
  vapply(x, function(s) {
    utils::URLdecode(gsub("%(?![0-9A-Fa-f]{2})", "%25", s, perl = TRUE))
  }, character(1), USE.NAMES = FALSE)
}

# Internal single-file reader; the public read_fasta_lengths() wrapper
# resolves `file` / `files` and combines the result.
read_fasta_lengths_single <- function(file, header_delim = NULL) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (is.character(file) && length(file) == 1 && !file.exists(file)) {
    ggchord_stop("read_fasta_lengths(): file not found: ", file, call. = FALSE)
  }
  lines <- ggchord_read_text(file, "read_fasta_lengths()")
  if (length(lines) == 0) {
    ggchord_stop("read_fasta_lengths(): file is empty", call. = FALSE)
  }
  is_header <- grepl("^>", lines)
  if (!any(is_header)) {
    ggchord_stop("read_fasta_lengths(): no FASTA headers (lines starting with '>') found",
         call. = FALSE)
  }
  if (any(nzchar(trimws(lines[seq_len(which(is_header)[1L] - 1L)]))))
    ggchord_stop("read_fasta_lengths(): sequence text before the first FASTA header in ", file)
  ids <- sub("^>", "", lines[is_header])
  ids <- trimws(ids)
  if (!is.null(header_delim)) {
    ids <- vapply(strsplit(ids, header_delim, fixed = TRUE),
                  function(x) x[1], character(1))
  } else {
    ids <- vapply(strsplit(ids, "[[:space:]]+"),
                  function(x) x[1], character(1))
  }
  if (anyNA(ids) || any(!nzchar(ids))) ggchord_stop(
    "read_fasta_lengths(): empty sequence identifier in ", file)

  header_pos <- which(is_header)
  starts <- header_pos + 1L
  ends <- c(header_pos[-1] - 1L, length(lines))
  lens <- vapply(seq_along(header_pos), function(i) {
    if (starts[i] > ends[i]) return(0L)
    sum(nchar(gsub("[[:space:]]+", "", lines[starts[i]:ends[i]])))
  }, integer(1))

  out <- data.frame(accver = ids, length = as.numeric(lens),
                    .source_line = header_pos,
                    stringsAsFactors = FALSE)
  dup <- unique(ids[duplicated(ids)])
  if (length(dup) > 0) {
    warning("read_fasta_lengths(): duplicate sequence IDs found: ",
            paste(dup, collapse = ", "), call. = FALSE)
  }
  out
}

#' Read one or more FASTA files into seq_data format
#'
#' `file` reads a single FASTA file; `files` reads multiple files at once and
#' combines them. `files` accepts a character vector of literal paths and/or
#' wildcard patterns (e.g. `"examples/fasta/*.fna"`).
#'
#' @param file Optional path vector or glob patterns; read before `files`.
#' @param files Optional character vector of FASTA files (literal paths and/or
#'   wildcard patterns). All matched files are read and combined.
#' @param header_delim Optional character. When given, each header is split at
#'   every occurrence of this delimiter and only the first piece is kept.
#' @param source_file Logical. Add a `.source_file` column when reading one or
#'   more files, default `FALSE`.
#'
#' @return A data.frame with columns `accver` and `length`.
#'   A `ggchord_import` attribute records files and per-file output row counts.
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#'
#' @examples
#' library(ggchord)
#' fasta <- tempfile(fileext = ".fna")
#' writeLines(c(">seqA some description", "ACGTACGTACGTACGT", "ACGTACGT",
#'              ">seqB", "TTTTGGGG"), fasta)
#' read_fasta_lengths(fasta)
read_fasta_lengths <- function(file = NULL, files = NULL, header_delim = NULL,
                               source_file = FALSE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  paths <- ggchord_resolve_files(file, files, "read_fasta_lengths")
  out <- lapply(paths, function(f) {
    data <- read_fasta_lengths_single(f, header_delim = header_delim)
    ggchord_add_source_file(data, f, source_file)
  })
  res <- ggchord_rbind_fill(out)
  if (!is.null(res) && anyDuplicated(res$accver)) {
    dup <- unique(res$accver[duplicated(res$accver)])
    warning("read_fasta_lengths(): duplicate sequence IDs found: ",
            paste(dup, collapse = ", "), call. = FALSE)
  }
  ggchord_format_import(res, paths, "FASTA", vapply(out, nrow, integer(1)),
    "sequence-lengths")
}
