# validate.R - structured data validation for ggchord inputs
#
# v0.7.0: validate_ggchord_data() returns a structured "ggchord_validation"
# object (valid flag, errors, warnings, per-category counts, data summary,
# original row numbers and cleanable issues). ggchord() calls it internally
# (validate = "warn" | "error" | "none") and caches the result on the plot.

# ---------------------------------------------------------------------------
# Internal issue collector
# ---------------------------------------------------------------------------

new_validation_collector <- function() {
  empty <- data.frame(
    table = character(0), category = character(0),
    row = integer(0), column = character(0),
    message = character(0), stringsAsFactors = FALSE
  )
  list(errors = empty, warnings = empty)
}

#' Add one issue (error or warning) to a validation collector
#' @keywords internal
add_validation_issue <- function(col, table, category, rows = NA_integer_,
                                 column = NA_character_, message,
                                 severity = c("error", "warning")) {
  severity <- match.arg(severity)
  df <- data.frame(
    table = rep(table, length(rows)),
    category = rep(category, length(rows)),
    row = as.integer(rows),
    column = rep(column, length(rows)),
    message = rep(message, length(rows)),
    stringsAsFactors = FALSE
  )
  col[[if (severity == "error") "errors" else "warnings"]] <-
    rbind(col[[if (severity == "error") "errors" else "warnings"]], df)
  col
}

# ---------------------------------------------------------------------------
# Per-table validators
# ---------------------------------------------------------------------------

validate_seq_data <- function(seq_data, col) {
  if (!is.data.frame(seq_data)) {
    return(add_validation_issue(
      col, "seq", "not_dataframe", NA_integer_, NA_character_,
      "seq_data must be a data.frame", "error"))
  }
  if (nrow(seq_data) == 0) {
    return(add_validation_issue(
      col, "seq", "empty", NA_integer_, NA_character_,
      "seq_data is empty; at least one sequence is required", "error"))
  }
  required <- c("seq_id", "length")
  missing <- setdiff(required, colnames(seq_data))
  if (length(missing) > 0) {
    return(add_validation_issue(
      col, "seq", "missing_columns", NA_integer_, NA_character_,
      sprintf("seq_data is missing required column(s): %s",
              paste(missing, collapse = ", ")), "error"))
  }

  sid <- seq_data$seq_id
  if (anyNA(sid)) {
    rows <- which(is.na(sid))
    col <- add_validation_issue(
      col, "seq", "missing_id", rows, "seq_id",
      "seq_id contains NA values", "error")
  }
  empty_ids <- !is.na(sid) & !nzchar(as.character(sid))
  if (any(empty_ids)) {
    rows <- which(empty_ids)
    col <- add_validation_issue(
      col, "seq", "empty_id", rows, "seq_id",
      "seq_id contains empty strings", "error")
  }
  dup <- unique(sid[duplicated(sid)])
  dup <- dup[!is.na(dup) & nzchar(as.character(dup))]
  if (length(dup) > 0) {
    rows <- which(sid %in% dup)
    col <- add_validation_issue(
      col, "seq", "duplicate_id", rows, "seq_id",
      sprintf("seq_id is duplicated (values: %s); seq_id must be unique",
              paste(dup, collapse = ", ")), "error")
  }

  if (!is.numeric(seq_data$length)) {
    col <- add_validation_issue(
      col, "seq", "invalid_length", NA_integer_, "length",
      "seq_data$length must be numeric", "error")
  } else {
    bad <- !is.finite(seq_data$length) | seq_data$length <= 0
    if (any(bad)) {
      rows <- which(bad)
      col <- add_validation_issue(
        col, "seq", "invalid_length", rows, "length",
        "seq_data$length must be finite positive numbers (NA/NaN/Inf/<= 0 found)",
        "error")
    }
    frac <- is.finite(seq_data$length) & seq_data$length > 0 &
      seq_data$length != round(seq_data$length)
    if (any(frac)) {
      rows <- which(frac)
      col <- add_validation_issue(
        col, "seq", "fractional_length", rows, "length",
        "seq_data$length is not an integer (fractional lengths may produce odd arc scales)",
        "warning")
    }
    ok <- is.finite(seq_data$length) & seq_data$length > 0
    if (sum(ok) >= 2) {
      lens <- seq_data$length[ok]
      ratio <- max(lens) / min(lens)
      if (ratio > 1000) {
        col <- add_validation_issue(
          col, "seq", "extreme_length", NA_integer_, "length",
          sprintf(paste0("Extreme sequence length differences (max/min = %.1f); ",
                         "short sequences will occupy tiny arcs"),
                  ratio), "warning")
      }
    }
  }
  col
}

validate_ribbon_data <- function(ribbon_data, seq_data, col,
                                 check_coordinates, check_duplicates,
                                 check_self_links) {
  if (is.null(ribbon_data)) return(col)
  if (!is.data.frame(ribbon_data)) {
    return(add_validation_issue(
      col, "ribbon", "not_dataframe", NA_integer_, NA_character_,
      "ribbon_data must be a data.frame", "error"))
  }
  required <- c("qaccver", "saccver", "length", "pident",
                "qstart", "qend", "sstart", "send")
  missing <- setdiff(required, colnames(ribbon_data))
  if (length(missing) > 0) {
    return(add_validation_issue(
      col, "ribbon", "missing_columns", NA_integer_, NA_character_,
      sprintf("ribbon_data is missing required column(s): %s",
              paste(missing, collapse = ", ")), "error"))
  }
  if (nrow(ribbon_data) == 0) return(col)

  num_cols <- c("length", "pident", "qstart", "qend", "sstart", "send")
  for (nm in num_cols) {
    if (!is.numeric(ribbon_data[[nm]])) {
      col <- add_validation_issue(
        col, "ribbon", "non_numeric", NA_integer_, nm,
        sprintf("ribbon_data$%s must be numeric", nm), "error")
    } else if (any(!is.finite(ribbon_data[[nm]]))) {
      rows <- which(!is.finite(ribbon_data[[nm]]))
      col <- add_validation_issue(
        col, "ribbon", "non_finite", rows, nm,
        sprintf("ribbon_data$%s contains NA/NaN/Inf values; all alignment coordinates must be finite", nm),
        "error")
    }
  }

  for (nm in c("qaccver", "saccver")) {
    x <- ribbon_data[[nm]]
    if (anyNA(x)) {
      rows <- which(is.na(x))
      col <- add_validation_issue(
        col, "ribbon", "missing_id", rows, nm,
        sprintf("ribbon_data$%s contains NA sequence IDs", nm), "error")
    }
    empty <- !is.na(x) & !nzchar(as.character(x))
    if (any(empty)) {
      rows <- which(empty)
      col <- add_validation_issue(
        col, "ribbon", "empty_id", rows, nm,
        sprintf("ribbon_data$%s contains empty sequence IDs", nm), "error")
    }
  }

  if (any(ribbon_data$length <= 0)) {
    rows <- which(ribbon_data$length <= 0)
    col <- add_validation_issue(
      col, "ribbon", "invalid_length", rows, "length",
      "ribbon_data$length must be positive", "error")
  }

  known_ids <- seq_data$seq_id
  if (is.data.frame(seq_data) && "seq_id" %in% colnames(seq_data)) {
    for (nm in c("qaccver", "saccver")) {
      x <- as.character(ribbon_data[[nm]])
      bad <- !is.na(x) & nzchar(x) & !(x %in% known_ids)
      if (any(bad)) {
        rows <- which(bad)
        unk <- sort(unique(x[bad]))
        col <- add_validation_issue(
          col, "ribbon", "unknown_id", rows, nm,
          sprintf(paste0("ribbon_data$%s references sequence IDs not present in ",
                         "seq_data (%s); these rows are silently skipped when drawing"),
                  nm, paste(unk, collapse = ", ")), "error")
      }
    }
  }

  if (isTRUE(check_self_links) && nrow(ribbon_data) > 0) {
    self <- !is.na(ribbon_data$qaccver) & !is.na(ribbon_data$saccver) &
      as.character(ribbon_data$qaccver) == as.character(ribbon_data$saccver)
    if (any(self)) {
      rows <- which(self)
      col <- add_validation_issue(
        col, "ribbon", "self_link", rows, "qaccver/saccver",
        "Self-links (qaccver == saccver) connect a sequence to itself and are not drawn",
        "error")
    }
  }

  # Reversed intervals are a legitimate representation of reverse-complement
  # hits in BLAST output (sstart > send), so they are a warning, not an error.
  rev <- ribbon_data$qstart > ribbon_data$qend |
    ribbon_data$sstart > ribbon_data$send
  rev[is.na(rev)] <- FALSE
  if (any(rev)) {
    rows <- which(rev)
    col <- add_validation_issue(
      col, "ribbon", "reversed_interval", rows, "qstart/qend/sstart/send",
      "start > end in ribbon coordinates (may represent a reverse-complement hit; clean_ggchord_data() can sort them for stable drawing)",
      "warning")
  }

  pid <- ribbon_data$pident
  pid_bad <- !is.na(pid) & (pid < 0 | pid > 100)
  if (any(pid_bad)) {
    rows <- which(pid_bad)
    col <- add_validation_issue(
      col, "ribbon", "invalid_pident", rows, "pident",
      "pident values outside [0, 100] (clean_ggchord_data(invalid_pident = 'clip') clamps them)",
      "warning")
  }

  if (isTRUE(check_coordinates) && is.data.frame(seq_data) &&
      all(c("seq_id", "length") %in% colnames(seq_data))) {
    lens <- stats::setNames(seq_data$length, seq_data$seq_id)
    coords <- list(
      qstart = "qaccver", qend = "qaccver",
      sstart = "saccver", send = "saccver"
    )
    for (nm in names(coords)) {
      idcol <- coords[[nm]]
      ids <- as.character(ribbon_data[[idcol]])
      vals <- ribbon_data[[nm]]
      known <- !is.na(ids) & ids %in% names(lens)
      if (!any(known)) next
      bad <- known & (is.na(vals) | vals < 1 | vals > lens[ids])
      bad[is.na(bad)] <- FALSE
      if (any(bad)) {
        rows <- which(bad)
        col <- add_validation_issue(
          col, "ribbon", "out_of_range", rows, nm,
          sprintf("ribbon_data$%s is outside [1, sequence length] (clean_ggchord_data(out_of_range = 'clip') clamps it)", nm),
          "error")
      }
    }
  }

  if (isTRUE(check_duplicates)) {
    dups <- validation_ribbon_duplicates(ribbon_data)
    if (nrow(dups$exact) > 0) {
      col <- add_validation_issue(
        col, "ribbon", "exact_duplicate", dups$exact$row, NA_character_,
        "Fully duplicated alignment rows (identical pair and coordinates)",
        "warning")
    }
    if (nrow(dups$near) > 0) {
      col <- add_validation_issue(
        col, "ribbon", "near_duplicate", dups$near$row, NA_character_,
        "Near-duplicate alignment rows (same pair, coordinates within 5 bp)",
        "warning")
    }
    if (nrow(dups$overlap) > 0) {
      col <- add_validation_issue(
        col, "ribbon", "overlapping", dups$overlap$row, NA_character_,
        "Highly overlapping alignment blocks (reciprocal overlap >= 0.9 on both sequences)",
        "warning")
    }
  }

  # Optional-field diagnostics (informative only; the columns are not required)
  if ("evalue" %in% colnames(ribbon_data)) {
    e <- ribbon_data$evalue
    if (is.numeric(e) && any(!is.na(e) & e < 0)) {
      rows <- which(!is.na(e) & e < 0)
      col <- add_validation_issue(
        col, "ribbon", "evalue_diag", rows, "evalue",
        "evalue should be non-negative (negative E-values are unusual)",
        "warning")
    }
  }
  if ("bitscore" %in% colnames(ribbon_data)) {
    b <- ribbon_data$bitscore
    if (is.numeric(b) && any(!is.na(b) & b < 0)) {
      rows <- which(!is.na(b) & b < 0)
      col <- add_validation_issue(
        col, "ribbon", "bitscore_diag", rows, "bitscore",
        "bitscore should be non-negative", "warning")
    }
  }
  if ("qcovs" %in% colnames(ribbon_data)) {
    qc <- ribbon_data$qcovs
    if (is.numeric(qc) && any(!is.na(qc) & (qc < 0 | qc > 100))) {
      rows <- which(!is.na(qc) & (qc < 0 | qc > 100))
      col <- add_validation_issue(
        col, "ribbon", "qcovs_diag", rows, "qcovs",
        "qcovs (query coverage) should be within [0, 100]", "warning")
    }
  }
  for (nm in c("qlen", "slen")) {
    if (nm %in% colnames(ribbon_data)) {
      v <- ribbon_data[[nm]]
      if (is.numeric(v) && any(!is.na(v) & v <= 0)) {
        rows <- which(!is.na(v) & v <= 0)
        col <- add_validation_issue(
          col, "ribbon", "qlen_slen_diag", rows, nm,
          sprintf("%s should be positive", nm), "warning")
      }
    }
  }
  col
}

# Detect exact / near / highly-overlapping duplicate ribbon rows.
# Returns a list(exact, near, overlap), each a data.frame(row, dup_of).
# Near/overlap checks run pairwise inside (qaccver, saccver) groups and are
# skipped for very large groups to keep the default validation fast.
validation_ribbon_duplicates <- function(ribbon_data, near_tol = 5,
                                         overlap_ratio = 0.9,
                                         max_group = 1000) {
  exact <- data.frame(row = integer(0), dup_of = integer(0))
  near <- data.frame(row = integer(0), dup_of = integer(0))
  overlap <- data.frame(row = integer(0), dup_of = integer(0))

  key <- paste(ribbon_data$qaccver, ribbon_data$saccver,
               ribbon_data$qstart, ribbon_data$qend,
               ribbon_data$sstart, ribbon_data$send, sep = "\r")
  dup_idx <- which(duplicated(key))
  if (length(dup_idx) > 0) {
    first <- match(key[dup_idx], key)
    exact <- data.frame(row = dup_idx, dup_of = first)
  }

  # pairwise checks within (q, s) groups
  pair <- paste(ribbon_data$qaccver, ribbon_data$saccver, sep = "\r")
  grp <- split(seq_len(nrow(ribbon_data)), pair)
  near_rows <- list()
  overlap_rows <- list()
  for (g in grp) {
    if (length(g) < 2) next
    if (length(g) > max_group) next
    qs <- ribbon_data$qstart[g]
    qe <- ribbon_data$qend[g]
    ss <- ribbon_data$sstart[g]
    se <- ribbon_data$send[g]
    finite_rows <- is.finite(qs) & is.finite(qe) &
      is.finite(ss) & is.finite(se)
    g <- g[finite_rows]
    qs <- qs[finite_rows]
    qe <- qe[finite_rows]
    ss <- ss[finite_rows]
    se <- se[finite_rows]
    if (length(g) < 2) next
    for (i in seq_len(length(g) - 1L)) {
      for (j in (i + 1L):length(g)) {
        d <- max(abs(qs[i] - qs[j]), abs(qe[i] - qe[j]),
                 abs(ss[i] - ss[j]), abs(se[i] - se[j]))
        if (d <= near_tol) {
          near_rows[[length(near_rows) + 1L]] <- c(g[j], g[i])
        }
        qr <- interval_recip_overlap(qs[i], qe[i], qs[j], qe[j])
        sr <- interval_recip_overlap(ss[i], se[i], ss[j], se[j])
        if (is.finite(qr) && is.finite(sr) &&
            qr >= overlap_ratio && sr >= overlap_ratio) {
          overlap_rows[[length(overlap_rows) + 1L]] <- c(g[j], g[i])
        }
      }
    }
  }
  if (length(near_rows) > 0) {
    m <- do.call(rbind, near_rows)
    near <- data.frame(row = m[, 1], dup_of = m[, 2])
  }
  if (length(overlap_rows) > 0) {
    m <- do.call(rbind, overlap_rows)
    overlap <- data.frame(row = m[, 1], dup_of = m[, 2])
  }
  list(exact = exact, near = near, overlap = overlap)
}

#' Reciprocal overlap of two 1-based closed intervals
#' @keywords internal
interval_recip_overlap <- function(a1, a2, b1, b2) {
  lo <- max(a1, b1)
  hi <- min(a2, b2)
  ovl <- max(0, hi - lo + 1)
  denom <- min(abs(a2 - a1) + 1, abs(b2 - b1) + 1)
  if (denom <= 0) return(0)
  ovl / denom
}

validate_gene_data <- function(gene_data, seq_data, col,
                               check_coordinates, check_duplicates) {
  if (is.null(gene_data)) return(col)
  if (!is.data.frame(gene_data)) {
    return(add_validation_issue(
      col, "gene", "not_dataframe", NA_integer_, NA_character_,
      "gene_data must be a data.frame", "error"))
  }
  required <- c("seq_id", "start", "end", "strand", "anno")
  missing <- setdiff(required, colnames(gene_data))
  if (length(missing) > 0) {
    return(add_validation_issue(
      col, "gene", "missing_columns", NA_integer_, NA_character_,
      sprintf("gene_data is missing required column(s): %s",
              paste(missing, collapse = ", ")), "error"))
  }
  if (nrow(gene_data) == 0) return(col)

  if (!is.numeric(gene_data$start) || !is.numeric(gene_data$end)) {
    col <- add_validation_issue(
      col, "gene", "non_numeric", NA_integer_, "start/end",
      "gene_data$start and $end must be numeric", "error")
  } else {
    for (nm in c("start", "end")) {
      if (any(!is.finite(gene_data[[nm]]))) {
        rows <- which(!is.finite(gene_data[[nm]]))
        col <- add_validation_issue(
          col, "gene", "non_finite", rows, nm,
          sprintf("gene_data$%s contains NA/NaN/Inf values", nm), "error")
      }
    }
  }

  sid <- gene_data$seq_id
  if (anyNA(sid)) {
    rows <- which(is.na(sid))
    col <- add_validation_issue(
      col, "gene", "missing_id", rows, "seq_id",
      "gene_data$seq_id contains NA values", "error")
  }
  empty <- !is.na(sid) & !nzchar(as.character(sid))
  if (any(empty)) {
    rows <- which(empty)
    col <- add_validation_issue(
      col, "gene", "empty_id", rows, "seq_id",
      "gene_data$seq_id contains empty strings", "error")
  }
  if (anyNA(gene_data$strand)) {
    rows <- which(is.na(gene_data$strand))
    col <- add_validation_issue(
      col, "gene", "invalid_strand", rows, "strand",
      "gene_data$strand contains NA; strand must be '+' or '-'", "error")
  }
  bad_strand <- !is.na(gene_data$strand) & !gene_data$strand %in% c("+", "-")
  if (any(bad_strand)) {
    rows <- which(bad_strand)
    col <- add_validation_issue(
      col, "gene", "invalid_strand", rows, "strand",
      sprintf("gene_data$strand contains invalid values (%s); strand must be '+' or '-'",
              paste(unique(gene_data$strand[bad_strand]), collapse = ", ")), "error")
  }

  if (is.data.frame(seq_data) && "seq_id" %in% colnames(seq_data)) {
    x <- as.character(sid)
    bad <- !is.na(x) & nzchar(x) & !(x %in% seq_data$seq_id)
    if (any(bad)) {
      rows <- which(bad)
      unk <- sort(unique(x[bad]))
      col <- add_validation_issue(
        col, "gene", "unknown_id", rows, "seq_id",
        sprintf(paste0("gene_data$seq_id references sequence IDs not present in ",
                       "seq_data (%s); these features are skipped when drawing"),
                paste(unk, collapse = ", ")), "warning")
    }
  }

  rev <- gene_data$start > gene_data$end
  rev[is.na(rev)] <- FALSE
  if (any(rev)) {
    rows <- which(rev)
    col <- add_validation_issue(
      col, "gene", "reversed_interval", rows, "start/end",
      "start > end in gene_data (drawn with min/max; clean_ggchord_data(reversed_interval = 'sort') can sort them)",
      "warning")
  }

  if (isTRUE(check_coordinates) && is.data.frame(seq_data) &&
      all(c("seq_id", "length") %in% colnames(seq_data))) {
    lens <- stats::setNames(seq_data$length, seq_data$seq_id)
    ids <- as.character(sid)
    known <- !is.na(ids) & ids %in% names(lens)
    for (nm in c("start", "end")) {
      vals <- gene_data[[nm]]
      bad <- known & (is.na(vals) | vals < 1 | vals > lens[ids])
      bad[is.na(bad)] <- FALSE
      if (any(bad)) {
        rows <- which(bad)
        col <- add_validation_issue(
          col, "gene", "out_of_range", rows, nm,
          sprintf("gene_data$%s is outside [1, sequence length] (clean_ggchord_data(out_of_range = 'clip') clamps it)", nm),
          "error")
      }
    }
  }

  anno <- gene_data$anno
  miss_anno <- is.na(anno) | !nzchar(as.character(anno))
  if (any(miss_anno)) {
    rows <- which(miss_anno)
    col <- add_validation_issue(
      col, "gene", "missing_anno", rows, "anno",
      "gene_data$anno is missing or empty (clean_ggchord_data(empty_annotation = 'replace') fills it)",
      "warning")
  }

  if (isTRUE(check_duplicates) && nrow(gene_data) > 0) {
    key <- paste(gene_data$seq_id, gene_data$start, gene_data$end,
                 gene_data$strand, gene_data$anno, sep = "\r")
    dup_idx <- which(duplicated(key))
    if (length(dup_idx) > 0) {
      first <- match(key[dup_idx], key)
      rows <- unique(c(dup_idx, first))
      col <- add_validation_issue(
        col, "gene", "exact_duplicate", rows, NA_character_,
        "Fully duplicated gene features (same seq_id, coordinates, strand and anno)",
        "warning")
    }
    # highly overlapping features on the same sequence
    grp <- split(seq_len(nrow(gene_data)), gene_data$seq_id)
    ov_rows <- integer(0)
    for (g in grp) {
      if (length(g) < 2 || length(g) > 1000) next
      st <- gene_data$start[g]
      en <- gene_data$end[g]
      for (i in seq_len(length(g) - 1L)) {
        for (j in (i + 1L):length(g)) {
          r <- interval_recip_overlap(st[i], en[i], st[j], en[j])
          if (is.finite(r) && r >= 0.9) ov_rows <- c(ov_rows, g[j])
        }
      }
    }
    if (length(ov_rows) > 0) {
      col <- add_validation_issue(
        col, "gene", "overlapping", unique(ov_rows), NA_character_,
        "Highly overlapping gene features on the same sequence (reciprocal overlap >= 0.9)",
        "warning")
    }
  }
  col
}

# ---------------------------------------------------------------------------
# Summary / data summary helpers
# ---------------------------------------------------------------------------

validation_summary_table <- function(errors, warnings) {
  counts <- function(df, sev) {
    if (nrow(df) == 0) return(NULL)
    agg <- stats::aggregate(seq_len(nrow(df)),
                            by = list(table = df$table, category = df$category),
                            FUN = length)
    names(agg)[3] <- "n"
    agg$severity <- sev
    agg[, c("table", "category", "severity", "n")]
  }
  tab <- rbind(counts(errors, "error"), counts(warnings, "warning"))
  if (is.null(tab)) {
    return(data.frame(table = character(0), category = character(0),
                      severity = character(0), n = integer(0)))
  }
  tab[order(tab$table, tab$category), , drop = FALSE]
}

validation_invalid_rows <- function(errors, warnings = NULL) {
  all_issues <- rbind(
    if (!is.null(errors) && nrow(errors) > 0) errors else NULL,
    if (!is.null(warnings) && nrow(warnings) > 0) warnings else NULL
  )
  if (is.null(all_issues) || nrow(all_issues) == 0) return(list())
  keep <- !is.na(all_issues$row)
  if (!any(keep)) return(list())
  issues <- all_issues[keep, , drop = FALSE]
  split(issues$row, issues$category)
}

validation_cleanable <- function(errors, warnings) {
  cats <- rbind(errors, warnings)
  if (nrow(cats) == 0) {
    return(data.frame(table = character(0), category = character(0),
                      rows = I(list()), suggestion = character(0)))
  }
  suggestions <- c(
    unknown_id = "clean_ggchord_data(unknown_id = 'drop' | 'error' | 'keep')",
    out_of_range = "clean_ggchord_data(out_of_range = 'clip' | 'drop' | 'error' | 'keep')",
    reversed_interval = "clean_ggchord_data(reversed_interval = 'sort' | 'drop' | 'error' | 'keep')",
    invalid_pident = "clean_ggchord_data(invalid_pident = 'clip' | 'drop' | 'error' | 'keep')",
    missing_anno = "clean_ggchord_data(empty_annotation = 'replace' | 'drop' | 'keep')",
    exact_duplicate = "deduplicate_ggchord_ribbons(by = 'exact')",
    near_duplicate = "deduplicate_ggchord_ribbons(by = 'coordinates')",
    overlapping = "deduplicate_ggchord_ribbons(by = 'overlap') or filter manually"
  )
  cats$suggestion <- suggestions[cats$category]
  cats <- cats[!is.na(cats$suggestion), , drop = FALSE]
  if (nrow(cats) == 0) {
    return(data.frame(table = character(0), category = character(0),
                      rows = I(list()), suggestion = character(0)))
  }
  rows <- lapply(seq_len(nrow(cats)), function(i) {
    r <- cats$row[i]
    if (is.na(r)) integer(0) else r
  })
  data.frame(table = cats$table, category = cats$category,
             rows = I(rows), suggestion = cats$suggestion,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Validate ggchord input data before plotting
#'
#' Performs structured validation of \code{seq_data}, \code{ribbon_data} and
#' \code{gene_data} so that problems can be found, understood and fixed before
#' plotting. The result is a \code{"ggchord_validation"} object with a
#' \code{valid} flag, \code{errors} (severe problems that make the plot
#' misleading), \code{warnings} (drawable but noteworthy issues), per-category
#' counts, a data summary, the original row numbers of every problem, and a
#' list of automatically fixable issues.
#'
#' @param seq_data data.frame/tibble, required. Basic sequence information
#'   (columns \code{seq_id}, \code{length}).
#' @param ribbon_data data.frame/tibble, optional. Alignment results (columns
#'   \code{qaccver}, \code{saccver}, \code{length}, \code{pident},
#'   \code{qstart}, \code{qend}, \code{sstart}, \code{send}).
#' @param gene_data data.frame/tibble, optional. Gene annotation data (columns
#'   \code{seq_id}, \code{start}, \code{end}, \code{strand}, \code{anno}).
#' @param strict Logical. When \code{TRUE}, stop with an error as soon as any
#'   severe problem is found. When \code{FALSE} (default), return the full
#'   diagnostic report without stopping.
#' @param check_coordinates Logical, default \code{TRUE}. Whether to check that
#'   ribbon/gene coordinates stay inside \code{[1, sequence length]}.
#' @param check_duplicates Logical, default \code{TRUE}. Whether to look for
#'   fully duplicated, near-duplicated and highly overlapping records.
#' @param check_self_links Logical, default \code{TRUE}. Whether to flag
#'   alignment rows where \code{qaccver == saccver}.
#'
#' @return A \code{"ggchord_validation"} object (a list) with at least:
#'   \describe{
#'     \item{\code{valid}}{Logical: \code{TRUE} when there are no severe errors.}
#'     \item{\code{errors}}{data.frame of severe issues (table, category, row,
#'       column, message).}
#'     \item{\code{warnings}}{data.frame of non-severe issues (same columns).}
#'     \item{\code{summary}}{Per-category counts (table, category, severity, n).}
#'     \item{\code{data_summary}}{Counts of sequences/ribbons/genes, unknown
#'       IDs, out-of-range rows, etc.}
#'     \item{\code{invalid_rows}}{Named list of original row numbers per
#'       problem category.}
#'     \item{\code{cleanable}}{data.frame of fixable issues with suggested
#'       actions.}
#'   }
#'   \code{print()} and \code{summary()} methods are provided.
#' @export
#'
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' data(gene_data_example)
#'
#' res <- validate_ggchord_data(seq_data_example, ribbon_data_example,
#'                              gene_data_example)
#' res$valid
#' print(res)
#' summary(res)
#'
#' # Introduce a problem: unknown sequence ID in the ribbons
#' bad <- transform(ribbon_data_example, saccver = "NOT_A_SEQUENCE")
#' v <- validate_ggchord_data(seq_data_example, bad)
#' v$invalid_rows$ribbon_unknown_id
validate_ggchord_data <- function(seq_data,
                                  ribbon_data = NULL,
                                  gene_data = NULL,
                                  strict = FALSE,
                                  check_coordinates = TRUE,
                                  check_duplicates = TRUE,
                                  check_self_links = TRUE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!is.logical(strict) || length(strict) != 1 || is.na(strict)) {
    ggchord_stop("validate_ggchord_data(): strict must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(check_coordinates) || length(check_coordinates) != 1 ||
      is.na(check_coordinates)) {
    ggchord_stop("validate_ggchord_data(): check_coordinates must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(check_duplicates) || length(check_duplicates) != 1 ||
      is.na(check_duplicates)) {
    ggchord_stop("validate_ggchord_data(): check_duplicates must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(check_self_links) || length(check_self_links) != 1 ||
      is.na(check_self_links)) {
    ggchord_stop("validate_ggchord_data(): check_self_links must be TRUE or FALSE", call. = FALSE)
  }

  col <- new_validation_collector()
  col <- validate_seq_data(seq_data, col)
  col <- validate_ribbon_data(ribbon_data, seq_data, col,
                              check_coordinates, check_duplicates,
                              check_self_links)
  col <- validate_gene_data(gene_data, seq_data, col,
                            check_coordinates, check_duplicates)

  errors <- col$errors
  warnings <- col$warnings
  valid <- nrow(errors) == 0

  ds <- list(
    n_seq = if (is.data.frame(seq_data)) nrow(seq_data) else NA_integer_,
    n_ribbon = if (is.data.frame(ribbon_data)) nrow(ribbon_data) else NA_integer_,
    n_gene = if (is.data.frame(gene_data)) nrow(gene_data) else NA_integer_,
    n_unknown_ribbon_ids = sum(errors$category == "unknown_id" &
                                 errors$table == "ribbon" & !is.na(errors$row)),
    n_unknown_gene_ids = sum(warnings$category == "unknown_id" &
                               warnings$table == "gene" & !is.na(warnings$row)),
    n_self_links = sum(errors$category == "self_link" & !is.na(errors$row)),
    n_out_of_range_ribbon = sum(errors$category == "out_of_range" &
                                  errors$table == "ribbon" & !is.na(errors$row)),
    n_out_of_range_gene = sum(errors$category == "out_of_range" &
                                errors$table == "gene" & !is.na(errors$row)),
    n_reversed_ribbon = sum(warnings$category == "reversed_interval" &
                              warnings$table == "ribbon" & !is.na(warnings$row)),
    n_reversed_gene = sum(warnings$category == "reversed_interval" &
                            warnings$table == "gene" & !is.na(warnings$row)),
    n_duplicates_ribbon = sum(errors$category %in% c("exact_duplicate", "near_duplicate",
                                                     "overlapping") &
                                errors$table == "ribbon" & !is.na(errors$row)) +
      sum(warnings$category %in% c("exact_duplicate", "near_duplicate", "overlapping") &
            warnings$table == "ribbon" & !is.na(warnings$row)),
    n_duplicates_gene = sum(warnings$category %in% c("exact_duplicate", "overlapping") &
                              warnings$table == "gene" & !is.na(warnings$row)),
    n_bad_strand = sum(errors$category == "invalid_strand" & !is.na(errors$row))
  )

  out <- structure(list(
    valid        = valid,
    errors       = errors,
    warnings     = warnings,
    summary      = validation_summary_table(errors, warnings),
    data_summary = ds,
    invalid_rows = validation_invalid_rows(errors, warnings),
    cleanable    = validation_cleanable(errors, warnings)
  ), class = "ggchord_validation")

  if (strict && !valid) {
    ggchord_stop(sprintf(paste0(
      "validate_ggchord_data(): %d severe error(s) found in the input data ",
      "(first error: %s). Run validate_ggchord_data(..., strict = FALSE) to ",
      "obtain the full report."),
      nrow(errors), errors$message[1]), call. = FALSE)
  }
  out
}

#' @export
print.ggchord_validation <- function(x, ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  cat("ggchord data validation\n")
  cat("========================\n")
  status <- if (x$valid) "VALID" else "INVALID"
  cat(sprintf("Result: %s (%d error(s), %d warning(s))\n",
              status, nrow(x$errors), nrow(x$warnings)))
  ds <- x$data_summary
  cat(sprintf("Sequences: %s | Ribbons: %s | Genes: %s\n",
              ifelse(is.na(ds$n_seq), "NA", ds$n_seq),
              ifelse(is.na(ds$n_ribbon), "NA", ds$n_ribbon),
              ifelse(is.na(ds$n_gene), "NA", ds$n_gene)))
  if (nrow(x$errors) > 0) {
    cat("\nSevere errors (fix before plotting):\n")
    print(utils::head(x$errors[, c("table", "row", "column", "message")], 10),
          row.names = FALSE)
    if (nrow(x$errors) > 10) {
      cat(sprintf("... and %d more error(s)\n", nrow(x$errors) - 10))
    }
  }
  if (nrow(x$warnings) > 0) {
    cat("\nWarnings (drawable but noteworthy):\n")
    print(utils::head(x$warnings[, c("table", "row", "column", "message")], 10),
          row.names = FALSE)
    if (nrow(x$warnings) > 10) {
      cat(sprintf("... and %d more warning(s)\n", nrow(x$warnings) - 10))
    }
  }
  if (!x$valid) {
    cat("\nRun summary(x) for per-category counts and clean_ggchord_data() to\n",
        "automatically fix the fixable issues.\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.ggchord_validation <- function(object, ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  out <- object$summary
  attr(out, "valid") <- object$valid
  attr(out, "data_summary") <- object$data_summary
  class(out) <- c("ggchord_validation_summary", class(out))
  out
}

#' @export
print.ggchord_validation_summary <- function(x, ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  cat("ggchord validation summary\n")
  cat(sprintf("Valid: %s\n", if (isTRUE(attr(x, "valid"))) "TRUE" else "FALSE"))
  ds <- attr(x, "data_summary")
  cat(sprintf("Sequences: %s | Ribbons: %s | Genes: %s\n",
              ifelse(is.na(ds$n_seq), "NA", ds$n_seq),
              ifelse(is.na(ds$n_ribbon), "NA", ds$n_ribbon),
              ifelse(is.na(ds$n_gene), "NA", ds$n_gene)))
  if (nrow(x) == 0) {
    cat("No issues found.\n")
  } else {
    print.data.frame(as.data.frame(x), row.names = FALSE)
  }
  invisible(x)
}
