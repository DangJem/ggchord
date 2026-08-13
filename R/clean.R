# clean.R - conservative, report-driven data cleaning for ggchord inputs
#
# v0.7.0: clean_ggchord_data() applies explicit, per-problem policies to
# seq/ribbon/gene data and returns the cleaned tables plus a full report of
# every change (original row number, reason, original/new values, action).
# The user's input objects are never modified.

#' Clean ggchord input data with explicit, report-driven policies
#'
#' Applies a conservative, predictable set of cleaning policies to the three
#' ggchord tables and returns the cleaned copies plus a full report. The
#' original data frames are never modified. Nothing is dropped silently: every
#' change (including every dropped row) is recorded in \code{report} with the
#' original row number, the reason, the original value(s) and the new value(s).
#'
#' @param seq_data data.frame/tibble, required. Must contain \code{seq_id} and
#'   \code{length}; used as the coordinate reference for clipping.
#' @param ribbon_data data.frame/tibble, optional. Alignment results.
#' @param gene_data data.frame/tibble, optional. Gene annotation data.
#' @param unknown_id Character, default \code{"drop"}. Policy for rows whose
#'   sequence ID is missing, empty or unknown: \code{"drop"} removes them,
#'   \code{"error"} stops with a message, \code{"keep"} leaves them unchanged.
#' @param out_of_range Character, default \code{"clip"}. Policy for coordinates
#'   outside \code{[1, sequence length]}: \code{"clip"} clamps them,
#'   \code{"drop"} removes the row, \code{"error"} stops, \code{"keep"} leaves
#'   them unchanged. After clipping, intervals that become invalid (degenerate)
#'   are removed explicitly and reported.
#' @param reversed_interval Character, default \code{"sort"}. Policy for
#'   intervals with \code{start > end}: \code{"sort"} swaps the endpoints so
#'   the feature draws stably (the original direction is recorded in the
#'   report), \code{"drop"} removes the row, \code{"error"} stops,
#'   \code{"keep"} leaves them unchanged.
#' @param invalid_pident Character, default \code{"clip"}. Policy for
#'   \code{pident} outside \code{[0, 100]}: \code{"clip"} clamps them,
#'   \code{"drop"} removes the row, \code{"error"} stops, \code{"keep"} leaves
#'   them unchanged.
#' @param empty_annotation Character, default \code{"keep"}. Policy for
#'   missing/empty \code{anno} in gene data: \code{"replace"} fills them with
#'   \code{replacement_annotation}, \code{"drop"} removes the row,
#'   \code{"keep"} leaves them unchanged.
#' @param replacement_annotation Character, default \code{"unannotated"}.
#'   Replacement annotation used when \code{empty_annotation = "replace"}.
#'
#' @return A list with four components: \code{seq_data}, \code{ribbon_data},
#'   \code{gene_data} (cleaned copies) and \code{report} (a data.frame with
#'   columns \code{table}, \code{row} (original row number), \code{column},
#'   \code{reason}, \code{original_value}, \code{new_value} and
#'   \code{action}).
#' @export
#'
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' data(gene_data_example)
#'
#' # Introduce a few typical problems
#' bad_r <- transform(ribbon_data_example,
#'                    qstart = pmin(qstart, 1),
#'                    pident = pmin(pident, 150))
#' bad_g <- transform(gene_data_example,
#'                    anno = ifelse(seq_len(nrow(gene_data_example)) == 1,
#'                                  NA_character_, anno))
#'
#' out <- clean_ggchord_data(seq_data_example, bad_r, bad_g)
#' head(out$report)
#' # The cleaned tables are ready for ggchord()
#' \donttest{
#' p <- ggchord(out$seq_data, out$ribbon_data, out$gene_data) +
#'   geom_seq() + geom_ribbon() + geom_gene()
#' }
clean_ggchord_data <- function(
    seq_data,
    ribbon_data = NULL,
    gene_data = NULL,
    unknown_id = c("drop", "error", "keep"),
    out_of_range = c("clip", "drop", "error", "keep"),
    reversed_interval = c("sort", "drop", "error", "keep"),
    invalid_pident = c("clip", "drop", "error", "keep"),
    empty_annotation = c("keep", "drop", "replace"),
    replacement_annotation = "unannotated") {
  unknown_id <- match.arg(unknown_id)
  out_of_range <- match.arg(out_of_range)
  reversed_interval <- match.arg(reversed_interval)
  invalid_pident <- match.arg(invalid_pident)
  empty_annotation <- match.arg(empty_annotation)

  if (!is.character(replacement_annotation) ||
      length(replacement_annotation) != 1 ||
      is.na(replacement_annotation)) {
    stop("clean_ggchord_data(): replacement_annotation must be a single non-NA character string",
         call. = FALSE)
  }

  # seq_data is the coordinate reference; it must itself be usable.
  if (!is.data.frame(seq_data) || nrow(seq_data) == 0) {
    stop("clean_ggchord_data(): seq_data must be a non-empty data.frame",
         call. = FALSE)
  }
  if (!all(c("seq_id", "length") %in% colnames(seq_data))) {
    stop("clean_ggchord_data(): seq_data must contain the columns seq_id and length",
         call. = FALSE)
  }
  if (anyNA(seq_data$seq_id) || any(!nzchar(as.character(seq_data$seq_id))) ||
      anyDuplicated(seq_data$seq_id)) {
    stop("clean_ggchord_data(): seq_data$seq_id must be non-missing, non-empty and unique",
         call. = FALSE)
  }
  if (!is.numeric(seq_data$length) || any(!is.finite(seq_data$length)) ||
      any(seq_data$length <= 0)) {
    stop("clean_ggchord_data(): seq_data$length must contain finite positive numbers",
         call. = FALSE)
  }
  seq_lens <- stats::setNames(seq_data$length, seq_data$seq_id)

  # Copies: the user's objects are never modified.
  seq_out <- seq_data
  ribbon_out <- if (is.data.frame(ribbon_data)) ribbon_data else NULL
  gene_out <- if (is.data.frame(gene_data)) gene_data else NULL

  report <- data.frame(
    table = character(0), row = integer(0), column = character(0),
    reason = character(0), original_value = character(0),
    new_value = character(0), action = character(0),
    stringsAsFactors = FALSE
  )
  report_chunk <- function(table, row, column, reason,
                           original_value, new_value, action) {
    report <<- rbind(report, data.frame(
      table = table, row = as.integer(row), column = column,
      reason = reason, original_value = original_value,
      new_value = new_value, action = action,
      stringsAsFactors = FALSE
    ))
  }

  # -------------------------------------------------------------------------
  # ribbon_data
  # -------------------------------------------------------------------------
  if (!is.null(ribbon_out)) {
    if (!is.data.frame(ribbon_out)) {
      stop("clean_ggchord_data(): ribbon_data must be a data.frame", call. = FALSE)
    }
    req <- c("qaccver", "saccver", "length", "pident",
             "qstart", "qend", "sstart", "send")
    missing <- setdiff(req, colnames(ribbon_out))
    if (length(missing) > 0) {
      stop("clean_ggchord_data(): ribbon_data is missing required column(s): ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    num_cols <- c("length", "pident", "qstart", "qend", "sstart", "send")
    for (nm in num_cols) {
      if (!is.numeric(ribbon_out[[nm]])) {
        stop("clean_ggchord_data(): ribbon_data$", nm, " must be numeric",
             call. = FALSE)
      }
      if (any(!is.finite(ribbon_out[[nm]]))) {
        stop("clean_ggchord_data(): ribbon_data$", nm,
             " contains NA/NaN/Inf values; clean non-finite values before cleaning",
             call. = FALSE)
      }
    }
    if (any(ribbon_out$length <= 0)) {
      stop("clean_ggchord_data(): ribbon_data$length must be positive",
           call. = FALSE)
    }
    if (anyNA(ribbon_out$qaccver) || anyNA(ribbon_out$saccver) ||
        any(!nzchar(as.character(ribbon_out$qaccver))) ||
        any(!nzchar(as.character(ribbon_out$saccver)))) {
      stop("clean_ggchord_data(): ribbon_data sequence IDs must be non-missing and non-empty",
           call. = FALSE)
    }

    drop <- rep(FALSE, nrow(ribbon_out))
    n <- nrow(ribbon_out)

    for (colname in c("qaccver", "saccver")) {
      id <- as.character(ribbon_out[[colname]])
      bad <- !id %in% names(seq_lens)
      if (any(bad)) {
        rows <- which(bad)
        if (unknown_id == "error") {
          stop(sprintf(paste0("clean_ggchord_data(): unknown sequence ID(s) in ",
                              "ribbon_data row(s) %s (column %s): %s; choose ",
                              "unknown_id = 'drop' or 'keep'"),
                       paste(rows, collapse = ", "), colname,
                       paste(unique(id[bad]), collapse = ", ")), call. = FALSE)
        } else if (unknown_id == "drop") {
          drop[bad] <- TRUE
          report_chunk("ribbon", rows, colname, "unknown sequence ID",
                       id[bad], NA_character_, "drop")
        } else {
          report_chunk("ribbon", rows, colname, "unknown sequence ID (kept)",
                       id[bad], id[bad], "keep")
        }
      }
    }

    # Coordinate out-of-range / reversed / degenerate handling, row by row.
    for (i in which(!drop)) {
      row <- ribbon_out[i, ]
      qid <- as.character(row$qaccver)
      sid_ <- as.character(row$saccver)
      if (!qid %in% names(seq_lens) || !sid_ %in% names(seq_lens)) {
        # unknown IDs were kept by policy; their coordinates cannot be
        # validated or clipped against a sequence length, so leave the row
        # untouched (it was already reported under unknown_id = 'keep').
        next
      }
      coords <- c(qstart = row$qstart, qend = row$qend,
                  sstart = row$sstart, send = row$send)
      idof <- c(qstart = qid, qend = qid, sstart = sid_, send = sid_)
      new_coords <- coords
      clipped <- FALSE
      for (nm in names(coords)) {
        len <- seq_lens[[idof[[nm]]]]
        val <- coords[[nm]]
        if (val < 1 || val > len) {
          if (out_of_range == "error") {
            stop(sprintf(paste0("clean_ggchord_data(): ribbon_data$%s row %d is ",
                                "outside [1, %d] (value %s); choose ",
                                "out_of_range = 'clip', 'drop' or 'keep'"),
                         nm, i, len, val), call. = FALSE)
          } else if (out_of_range == "drop") {
            drop[i] <- TRUE
            report_chunk("ribbon", i, nm, "coordinate out of range",
                         sprintf("%s=%s", nm, val), NA_character_, "drop")
            break
          } else if (out_of_range == "clip") {
            new_coords[[nm]] <- min(max(val, 1), len)
            clipped <- TRUE
          } else {
            report_chunk("ribbon", i, nm, "coordinate out of range (kept)",
                         sprintf("%s=%s", nm, val), sprintf("%s=%s", nm, val),
                         "keep")
          }
        }
      }
      if (drop[i]) next
      if (clipped) {
        report_chunk("ribbon", i, paste(names(coords)[new_coords != coords],
                                        collapse = ","),
                     "coordinates clipped to [1, sequence length]",
                     paste(sprintf("%s=%s", names(coords), coords), collapse = ", "),
                     paste(sprintf("%s=%s", names(new_coords), new_coords),
                           collapse = ", "),
                     "clip")
      }
      coords <- new_coords

      # Reversed intervals
      if (coords[["qstart"]] > coords[["qend"]] ||
          coords[["sstart"]] > coords[["send"]]) {
        if (reversed_interval == "error") {
          stop(sprintf(paste0("clean_ggchord_data(): reversed interval in ",
                              "ribbon_data row %d (start > end); choose ",
                              "reversed_interval = 'sort', 'drop' or 'keep'"),
                       i), call. = FALSE)
        } else if (reversed_interval == "drop") {
          drop[i] <- TRUE
          report_chunk("ribbon", i, "qstart/qend/sstart/send",
                       "reversed interval (start > end)",
                       paste(sprintf("%s=%s", names(coords), coords),
                             collapse = ", "),
                       NA_character_, "drop")
          next
        } else if (reversed_interval == "sort") {
          report_chunk("ribbon", i, "qstart/qend/sstart/send",
                       "reversed interval sorted (original direction recorded)",
                       paste(sprintf("%s=%s", names(coords), coords),
                             collapse = ", "),
                       paste(sprintf("%s=%s", names(coords), coords),
                             collapse = ", "),
                       "sort")
          if (coords[["qstart"]] > coords[["qend"]]) {
            tmp <- coords[["qstart"]]; coords[["qstart"]] <- coords[["qend"]]
            coords[["qend"]] <- tmp
          }
          if (coords[["sstart"]] > coords[["send"]]) {
            tmp <- coords[["sstart"]]; coords[["sstart"]] <- coords[["send"]]
            coords[["send"]] <- tmp
          }
        } else {
          report_chunk("ribbon", i, "qstart/qend/sstart/send",
                       "reversed interval (kept)",
                       paste(sprintf("%s=%s", names(coords), coords),
                             collapse = ", "),
                       paste(sprintf("%s=%s", names(coords), coords),
                             collapse = ", "),
                       "keep")
        }
      }

      # Degenerate interval after clipping (cannot be drawn meaningfully).
      if (coords[["qstart"]] == coords[["qend"]] ||
          coords[["sstart"]] == coords[["send"]]) {
        drop[i] <- TRUE
        report_chunk("ribbon", i, "qstart/qend/sstart/send",
                     "degenerate (zero-length) interval after clipping",
                     paste(sprintf("%s=%s", names(coords), coords),
                           collapse = ", "),
                     NA_character_, "drop")
        next
      }

      # pident
      p <- row$pident
      if (p < 0 || p > 100) {
        if (invalid_pident == "error") {
          stop(sprintf("clean_ggchord_data(): ribbon_data$pident row %d is %s; choose invalid_pident = 'clip', 'drop' or 'keep'",
                       i, p), call. = FALSE)
        } else if (invalid_pident == "drop") {
          drop[i] <- TRUE
          report_chunk("ribbon", i, "pident", "pident outside [0, 100]",
                       as.character(p), NA_character_, "drop")
          next
        } else if (invalid_pident == "clip") {
          report_chunk("ribbon", i, "pident", "pident clipped to [0, 100]",
                       as.character(p), as.character(min(max(p, 0), 100)),
                       "clip")
          p <- min(max(p, 0), 100)
        } else {
          report_chunk("ribbon", i, "pident", "pident outside [0, 100] (kept)",
                       as.character(p), as.character(p), "keep")
        }
      }

      # write back
      for (nm in names(coords)) ribbon_out[i, nm] <- coords[[nm]]
      ribbon_out[i, "pident"] <- p
    }

    ribbon_out <- ribbon_out[!drop, , drop = FALSE]
  }

  # -------------------------------------------------------------------------
  # gene_data
  # -------------------------------------------------------------------------
  if (!is.null(gene_out)) {
    if (!is.data.frame(gene_out)) {
      stop("clean_ggchord_data(): gene_data must be a data.frame", call. = FALSE)
    }
    req <- c("seq_id", "start", "end", "strand", "anno")
    missing <- setdiff(req, colnames(gene_out))
    if (length(missing) > 0) {
      stop("clean_ggchord_data(): gene_data is missing required column(s): ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    if (!is.numeric(gene_out$start) || !is.numeric(gene_out$end) ||
        any(!is.finite(gene_out$start)) || any(!is.finite(gene_out$end))) {
      stop("clean_ggchord_data(): gene_data$start and $end must contain finite numbers",
           call. = FALSE)
    }
    if (anyNA(gene_out$strand) || any(!gene_out$strand %in% c("+", "-"))) {
      stop("clean_ggchord_data(): gene_data$strand can only be '+' or '-'",
           call. = FALSE)
    }
    if (anyNA(gene_out$seq_id) || any(!nzchar(as.character(gene_out$seq_id)))) {
      stop("clean_ggchord_data(): gene_data$seq_id must be non-missing and non-empty",
           call. = FALSE)
    }

    drop <- rep(FALSE, nrow(gene_out))
    n <- nrow(gene_out)

    id <- as.character(gene_out$seq_id)
    bad <- !id %in% names(seq_lens)
    if (any(bad)) {
      rows <- which(bad)
      if (unknown_id == "error") {
        stop(sprintf(paste0("clean_ggchord_data(): unknown sequence ID(s) in ",
                            "gene_data row(s) %s (column seq_id): %s; choose ",
                            "unknown_id = 'drop' or 'keep'"),
                     paste(rows, collapse = ", "),
                     paste(unique(id[bad]), collapse = ", ")), call. = FALSE)
      } else if (unknown_id == "drop") {
        drop[bad] <- TRUE
        report_chunk("gene", rows, "seq_id", "unknown sequence ID",
                     id[bad], NA_character_, "drop")
      } else {
        report_chunk("gene", rows, "seq_id", "unknown sequence ID (kept)",
                     id[bad], id[bad], "keep")
      }
    }

    for (i in which(!drop)) {
      row <- gene_out[i, ]
      sid_ <- as.character(row$seq_id)
      len <- seq_lens[[sid_]]
      st <- row$start
      en <- row$end
      changed <- FALSE

      if (st < 1 || st > len) {
        if (out_of_range == "error") {
          stop(sprintf("clean_ggchord_data(): gene_data$start row %d is %s (outside [1, %d]); choose out_of_range = 'clip', 'drop' or 'keep'",
                       i, st, len), call. = FALSE)
        } else if (out_of_range == "drop") {
          drop[i] <- TRUE
          report_chunk("gene", i, "start", "coordinate out of range",
                       as.character(st), NA_character_, "drop")
          next
        } else if (out_of_range == "clip") {
          report_chunk("gene", i, "start", "coordinate clipped to [1, length]",
                       as.character(st), as.character(min(max(st, 1), len)),
                       "clip")
          st <- min(max(st, 1), len)
          changed <- TRUE
        } else {
          report_chunk("gene", i, "start", "coordinate out of range (kept)",
                       as.character(st), as.character(st), "keep")
        }
      }
      if (en < 1 || en > len) {
        if (out_of_range == "error") {
          stop(sprintf("clean_ggchord_data(): gene_data$end row %d is %s (outside [1, %d]); choose out_of_range = 'clip', 'drop' or 'keep'",
                       i, en, len), call. = FALSE)
        } else if (out_of_range == "drop") {
          drop[i] <- TRUE
          report_chunk("gene", i, "end", "coordinate out of range",
                       as.character(en), NA_character_, "drop")
          next
        } else if (out_of_range == "clip") {
          report_chunk("gene", i, "end", "coordinate clipped to [1, length]",
                       as.character(en), as.character(min(max(en, 1), len)),
                       "clip")
          en <- min(max(en, 1), len)
          changed <- TRUE
        } else {
          report_chunk("gene", i, "end", "coordinate out of range (kept)",
                       as.character(en), as.character(en), "keep")
        }
      }
      if (drop[i]) next

      if (st > en) {
        if (reversed_interval == "error") {
          stop(sprintf("clean_ggchord_data(): reversed interval in gene_data row %d (start > end); choose reversed_interval = 'sort', 'drop' or 'keep'",
                       i), call. = FALSE)
        } else if (reversed_interval == "drop") {
          drop[i] <- TRUE
          report_chunk("gene", i, "start/end", "reversed interval (start > end)",
                       sprintf("start=%s, end=%s", st, en), NA_character_, "drop")
          next
        } else if (reversed_interval == "sort") {
          report_chunk("gene", i, "start/end",
                       "reversed interval sorted (original direction recorded)",
                       sprintf("start=%s, end=%s", st, en),
                       sprintf("start=%s, end=%s", en, st), "sort")
          tmp <- st; st <- en; en <- tmp
          changed <- TRUE
        } else {
          report_chunk("gene", i, "start/end", "reversed interval (kept)",
                       sprintf("start=%s, end=%s", st, en),
                       sprintf("start=%s, end=%s", st, en), "keep")
        }
      }
      if (st == en) {
        drop[i] <- TRUE
        report_chunk("gene", i, "start/end",
                     "degenerate (zero-length) feature after clipping",
                     sprintf("start=%s, end=%s", st, en), NA_character_, "drop")
        next
      }

      anno <- as.character(row$anno)
      if (is.na(anno) || !nzchar(anno)) {
        if (empty_annotation == "drop") {
          drop[i] <- TRUE
          report_chunk("gene", i, "anno", "missing or empty annotation",
                       ifelse(is.na(anno), NA_character_, ""),
                       NA_character_, "drop")
          next
        } else if (empty_annotation == "replace") {
          report_chunk("gene", i, "anno", "missing or empty annotation replaced",
                       ifelse(is.na(anno), NA_character_, ""),
                       replacement_annotation, "replace")
          anno <- replacement_annotation
          changed <- TRUE
        } else {
          report_chunk("gene", i, "anno", "missing or empty annotation (kept)",
                       ifelse(is.na(anno), NA_character_, ""),
                       ifelse(is.na(anno), NA_character_, ""), "keep")
        }
      }

      if (changed) {
        gene_out[i, "start"] <- st
        gene_out[i, "end"] <- en
        gene_out[i, "anno"] <- anno
      }
    }

    gene_out <- gene_out[!drop, , drop = FALSE]
  }

  list(
    seq_data = seq_out,
    ribbon_data = ribbon_out,
    gene_data = gene_out,
    report = report
  )
}
