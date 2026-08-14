# data-import.R - lightweight import helpers for common bioinformatics formats
#
# v0.8.0: read_blast(), read_gff3() and read_fasta_lengths() turn common file
# formats into the tidy data frames ggchord() expects, with no tidyverse
# dependency.  External command-line preparation (BLAST, seqkit, ...) stays
# outside the package; these helpers cover the "read the file in R" step.

# Resolve `file` and `files` into a character vector of paths. `files` may
# contain wildcard patterns (e.g. "examples/fasta/*.fna") and/or literal paths.
ggchord_resolve_files <- function(file, files, fun) {
  if (is.null(files)) {
    if (is.null(file) || length(file) == 0) {
      ggchord_stop(sprintf("%s(): provide `file` or `files`", fun),
           call. = FALSE)
    }
    if (!is.character(file)) {
      ggchord_stop(sprintf("%s(): `file` must be a character path", fun),
           call. = FALSE)
    }
    paths <- as.character(file)
  } else {
    if (!is.character(files)) {
      ggchord_stop(sprintf("%s(): `files` must be a character vector", fun),
           call. = FALSE)
    }
    paths <- as.character(files)
  }

  expanded <- unlist(lapply(paths, function(x) {
    if (grepl("[*?[]", x)) {
      g <- Sys.glob(x)
      if (length(g) == 0) {
        ggchord_stop(sprintf("%s(): no files match pattern: %s", fun, x),
             call. = FALSE)
      }
      g
    } else {
      if (!file.exists(x)) {
        ggchord_stop(sprintf("%s(): file not found: %s", fun, x),
             call. = FALSE)
      }
      x
    }
  }), use.names = FALSE)

  # Allow both `file` and `files` to be supplied; `file` is read first.
  if (!is.null(file) && length(file) > 0 && !is.null(files)) {
    expanded <- unique(c(as.character(file), expanded))
  }
  if (length(expanded) == 0) {
    ggchord_stop(sprintf("%s(): no files to read", fun), call. = FALSE)
  }
  expanded
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
        d[[m]] <- if (!is.null(src)) src[[m]][NA_integer_] else NA
      }
    }
    d[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
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

#' Read one or more BLAST tabular output files into ribbon_data format
#'
#' `file` reads a single BLAST tabular output file; `files` reads multiple
#' files at once and combines them. `files` accepts a character vector of
#' literal paths and/or wildcard patterns (e.g. `"examples/blastn/*.o7"`).
#'
#' @param file Optional path to a single BLAST tabular output file.
#' @param files Optional character vector of BLAST tabular output files
#'   (literal paths and/or wildcard patterns). All matched files are read and
#'   combined into one data.frame.
#' @param format Character. `"auto"` (default) detects the column layout
#'   from the number of columns; `"outfmt6"` / `"outfmt7"` require the
#'   standard 12/17-column layouts; `"custom"` requires `col_names`.
#' @param col_names Optional character vector naming the columns in the file,
#'   used with `format = "custom"` or to override auto-detection.
#' @param comment Character comment character, default `"#"` (BLAST
#'   outfmt 7 header lines start with `#`).
#' @param ... Additional arguments passed to [utils::read.table()]
#'   (e.g. `na.strings`).
#'
#' @return A data.frame with the required ribbon columns first
#'   (`qaccver`, `saccver`, `length`, `pident`, `qstart`, `qend`,
#'   `sstart`, `send`) followed by any preserved optional columns.
#' @export
read_blast <- function(file = NULL, files = NULL,
                       format = c("auto", "outfmt6", "outfmt7", "custom"),
                       col_names = NULL, comment = "#", ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  format <- match.arg(format)
  paths <- ggchord_resolve_files(file, files, "read_blast")
  out <- lapply(paths, function(f) {
    read_blast_single(f, format = format, col_names = col_names,
                      comment = comment, ...)
  })
  ggchord_rbind_fill(out)
}

# Internal single-file reader; the public read_gff3() wrapper resolves
# `file` / `files` and combines the result.
read_gff3_single <- function(file, feature_types = "CDS",
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

#' Read one or more GFF3 files into gene_data format
#'
#' `file` reads a single GFF3 file; `files` reads multiple files at once and
#' combines them. `files` accepts a character vector of literal paths and/or
#' wildcard patterns (e.g. `"examples/gff3/*.gff3"`).
#'
#' @param file Optional path to a single GFF3 file.
#' @param files Optional character vector of GFF3 files (literal paths and/or
#'   wildcard patterns). All matched files are read and combined.
#' @param feature_types Character vector of feature types to keep, default
#'   `"CDS"`.
#' @param anno_from Character vector of GFF3 attribute keys, tried in order to
#'   fill the `anno` column.
#' @param unstranded Character, default `"plus"`.
#'
#' @return A data.frame with `seq_id`, `start`, `end`, `strand`, `anno`
#'   followed by `type`, `source`, `score`, `phase` and `attributes`.
#' @export
read_gff3 <- function(file = NULL, files = NULL, feature_types = "CDS",
                      anno_from = c("product", "Name", "gene", "ID"),
                      unstranded = c("plus", "drop")) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  unstranded <- match.arg(unstranded)
  paths <- ggchord_resolve_files(file, files, "read_gff3")
  out <- lapply(paths, function(f) {
    read_gff3_single(f, feature_types = feature_types,
                     anno_from = anno_from, unstranded = unstranded)
  })
  ggchord_rbind_fill(out)
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

# Internal single-file reader; the public read_fasta_lengths() wrapper
# resolves `file` / `files` and combines the result.
read_fasta_lengths_single <- function(file, header_delim = NULL) {
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

#' Read one or more FASTA files into seq_data format
#'
#' `file` reads a single FASTA file; `files` reads multiple files at once and
#' combines them. `files` accepts a character vector of literal paths and/or
#' wildcard patterns (e.g. `"examples/fasta/*.fna"`).
#'
#' @param file Optional path to a single FASTA file.
#' @param files Optional character vector of FASTA files (literal paths and/or
#'   wildcard patterns). All matched files are read and combined.
#' @param header_delim Optional character. When given, each header is split at
#'   every occurrence of this delimiter and only the first piece is kept.
#'
#' @return A data.frame with columns `seq_id` and `length`.
#' @export
read_fasta_lengths <- function(file = NULL, files = NULL, header_delim = NULL) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  paths <- ggchord_resolve_files(file, files, "read_fasta_lengths")
  out <- lapply(paths, function(f) {
    read_fasta_lengths_single(f, header_delim = header_delim)
  })
  res <- ggchord_rbind_fill(out)
  if (!is.null(res) && anyDuplicated(res$seq_id)) {
    dup <- unique(res$seq_id[duplicated(res$seq_id)])
    warning("read_fasta_lengths(): duplicate sequence IDs found: ",
            paste(dup, collapse = ", "), call. = FALSE)
  }
  res
}
