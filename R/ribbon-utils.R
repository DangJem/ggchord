# ribbon-utils.R - ribbon filtering, deduplication and merging
#
# v0.8.0: filter_ggchord_ribbons(), deduplicate_ggchord_ribbons() and
# merge_ggchord_ribbons() prepare alignment tables before plotting. All three
# keep the original column order plus extra columns, return a report, and
# attach the original row numbers as the "source_rows" attribute of the data.

#' Filter alignment ribbons before plotting
#'
#' Keeps the ribbon rows that satisfy all requested criteria and optionally
#' sorts them. The returned data frame keeps every extra column and the
#' original column order; the original row numbers are attached as the
#' \code{source_rows} attribute.
#'
#' @param ribbon_data data.frame with at least \code{qaccver} and
#'   \code{saccver} columns (normally the full ribbon_data format).
#' @param seq_ids Optional character vector. Keep only rows where both the
#'   query and the subject are in \code{seq_ids}.
#' @param min_pident,max_pident Optional numeric. Lower/upper bounds on
#'   \code{pident}.
#' @param min_length Optional numeric. Lower bound on \code{length}.
#' @param max_evalue Optional numeric. Upper bound on \code{evalue} (requires
#'   an \code{evalue} column).
#' @param min_bitscore Optional numeric. Lower bound on \code{bitscore}
#'   (requires a \code{bitscore} column).
#' @param min_query_coverage Optional numeric (0-100). Lower bound on the query
#'   coverage. Uses the \code{qcovs} column when present; otherwise computed
#'   from \code{qlen} and the query interval (requires \code{qlen}).
#' @param min_subject_coverage Optional numeric (0-100). Lower bound on the
#'   subject coverage, computed from \code{slen} and the subject interval
#'   (requires an \code{slen} column).
#' @param keep_pairs Optional data.frame/list/matrix describing an undirected
#'   set of sequence pairs. A data.frame or matrix with the first two columns
#'   used as query/subject IDs, or a list of length-2 character vectors. A row
#'   is kept when its query/subject pair matches any pair in either direction.
#' @param drop_self_links Logical, default \code{TRUE}. Remove rows where
#'   \code{qaccver == saccver}.
#' @param sort_by Optional character vector of column names. Prefix a name
#'   with \code{"-"} or \code{"desc:"} for descending order (e.g.
#'   \code{c("pident", "-evalue")}).
#'
#' @return A list with \code{data} (the filtered data frame, with
#'   \code{source_rows} attribute) and \code{report} (n_input/n_kept/n_removed,
#'   removed_by_reason, removed_rows and kept_rows).
#' @export
#'
#' @examples
#' library(ggchord)
#' data(ribbon_data_example)
#' out <- filter_ggchord_ribbons(
#'   ribbon_data_example,
#'   min_pident = 95,
#'   drop_self_links = TRUE,
#'   sort_by = c("pident", "-length")
#' )
#' out$report
filter_ggchord_ribbons <- function(
    ribbon_data,
    seq_ids = NULL,
    min_pident = NULL,
    max_pident = NULL,
    min_length = NULL,
    max_evalue = NULL,
    min_bitscore = NULL,
    min_query_coverage = NULL,
    min_subject_coverage = NULL,
    keep_pairs = NULL,
    drop_self_links = TRUE,
    sort_by = NULL) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!is.data.frame(ribbon_data)) {
    ggchord_stop("filter_ggchord_ribbons(): ribbon_data must be a data.frame",
         call. = FALSE)
  }
  req <- c("qaccver", "saccver")
  missing <- setdiff(req, colnames(ribbon_data))
  if (length(missing) > 0) {
    ggchord_stop("filter_ggchord_ribbons(): ribbon_data is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  n_input <- nrow(ribbon_data)
  reason <- rep(NA_character_, n_input)

  require_col <- function(col, arg) {
    if (!col %in% colnames(ribbon_data)) {
      ggchord_stop(sprintf(paste0(
        "filter_ggchord_ribbons(): column '%s' is required for the '%s' ",
        "filter but is missing from ribbon_data"), col, arg), call. = FALSE)
    }
  }

  if (!is.null(seq_ids)) {
    idx <- which(!(ribbon_data$qaccver %in% seq_ids &
                     ribbon_data$saccver %in% seq_ids))
    if (length(idx) > 0) reason[idx] <- "seq_ids"
  }
  if (!is.null(min_pident)) {
    require_col("pident", "min_pident")
    idx <- which(is.na(ribbon_data$pident) | ribbon_data$pident < min_pident)
    if (length(idx) > 0) reason[idx] <- "min_pident"
  }
  if (!is.null(max_pident)) {
    require_col("pident", "max_pident")
    idx <- which(is.na(ribbon_data$pident) | ribbon_data$pident > max_pident)
    if (length(idx) > 0) reason[idx] <- "max_pident"
  }
  if (!is.null(min_length)) {
    require_col("length", "min_length")
    idx <- which(is.na(ribbon_data$length) | ribbon_data$length < min_length)
    if (length(idx) > 0) reason[idx] <- "min_length"
  }
  if (!is.null(max_evalue)) {
    require_col("evalue", "max_evalue")
    idx <- which(is.na(ribbon_data$evalue) | ribbon_data$evalue > max_evalue)
    if (length(idx) > 0) reason[idx] <- "max_evalue"
  }
  if (!is.null(min_bitscore)) {
    require_col("bitscore", "min_bitscore")
    idx <- which(is.na(ribbon_data$bitscore) | ribbon_data$bitscore < min_bitscore)
    if (length(idx) > 0) reason[idx] <- "min_bitscore"
  }
  if (!is.null(min_query_coverage)) {
    if ("qcovs" %in% colnames(ribbon_data)) {
      cov <- ribbon_data$qcovs
    } else {
      require_col("qlen", "min_query_coverage")
      cov <- ifelse(
        is.na(ribbon_data$qlen) | ribbon_data$qlen <= 0, NA_real_,
        (abs(ribbon_data$qend - ribbon_data$qstart) + 1) /
          ribbon_data$qlen * 100
      )
    }
    idx <- which(is.na(cov) | cov < min_query_coverage)
    if (length(idx) > 0) reason[idx] <- "min_query_coverage"
  }
  if (!is.null(min_subject_coverage)) {
    require_col("slen", "min_subject_coverage")
    cov <- ifelse(
      is.na(ribbon_data$slen) | ribbon_data$slen <= 0, NA_real_,
      (abs(ribbon_data$send - ribbon_data$sstart) + 1) /
        ribbon_data$slen * 100
    )
    idx <- which(is.na(cov) | cov < min_subject_coverage)
    if (length(idx) > 0) reason[idx] <- "min_subject_coverage"
  }
  if (isTRUE(drop_self_links)) {
    idx <- which(ribbon_data$qaccver == ribbon_data$saccver)
    if (length(idx) > 0) reason[idx] <- "self_link"
  }
  if (!is.null(keep_pairs)) {
    pairs <- normalize_keep_pairs(keep_pairs)
    keep_row <- vapply(seq_len(n_input), function(i) {
      q <- as.character(ribbon_data$qaccver[i])
      s <- as.character(ribbon_data$saccver[i])
      any((pairs$q == q & pairs$s == s) | (pairs$q == s & pairs$s == q))
    }, logical(1))
    idx <- which(!keep_row)
    if (length(idx) > 0) reason[idx] <- "keep_pairs"
  }

  kept <- is.na(reason)
  removed_rows <- which(!kept)
  out <- ribbon_data[kept, , drop = FALSE]
  kept_src <- which(kept)

  sort_info <- NULL
  if (!is.null(sort_by)) {
    if (!is.character(sort_by) || length(sort_by) == 0) {
      ggchord_stop("filter_ggchord_ribbons(): sort_by must be a non-empty character vector",
           call. = FALSE)
    }
    keys <- lapply(sort_by, function(nm) {
      desc <- FALSE
      col <- nm
      if (grepl("^-", nm)) {
        desc <- TRUE
        col <- sub("^-", "", nm)
      } else if (grepl("^desc:", nm)) {
        desc <- TRUE
        col <- sub("^desc:", "", nm)
      }
      if (!col %in% colnames(out)) {
        ggchord_stop("filter_ggchord_ribbons(): sort_by column not found: ", col,
             call. = FALSE)
      }
      list(col = col, desc = desc)
    })
    ord_args <- lapply(keys, function(k) {
      v <- out[[k$col]]
      if (is.character(v)) v <- factor(v, levels = unique(v))
      if (k$desc) -xtfrm(v) else xtfrm(v)
    })
    ord <- do.call(order, ord_args)
    out <- out[ord, , drop = FALSE]
    kept_src <- kept_src[ord]
    sort_info <- sort_by
  }

  attr(out, "source_rows") <- kept_src

  removed_reasons <- reason[!is.na(reason)]
  if (length(removed_reasons) == 0) {
    removed_by_reason <- data.frame(reason = character(0), n = integer(0))
  } else {
    removed_by_reason <- as.data.frame(table(removed_reasons),
                                       stringsAsFactors = FALSE)
    names(removed_by_reason) <- c("reason", "n")
  }

  list(
    data = out,
    report = list(
      n_input = n_input,
      n_kept = sum(kept),
      n_removed = sum(!kept),
      removed_by_reason = removed_by_reason,
      removed_rows = removed_rows,
      kept_rows = which(kept),
      sort_by = sort_info
    )
  )
}

#' Normalize the keep_pairs argument into a data.frame(q, s)
#' @keywords internal
normalize_keep_pairs <- function(pairs) {
  if (is.data.frame(pairs)) {
    if (ncol(pairs) < 2) {
      ggchord_stop("keep_pairs data.frame must have at least two columns (query, subject)",
           call. = FALSE)
    }
    return(data.frame(
      q = as.character(pairs[[1]]),
      s = as.character(pairs[[2]]),
      stringsAsFactors = FALSE
    ))
  }
  if (is.matrix(pairs)) {
    if (ncol(pairs) < 2) {
      ggchord_stop("keep_pairs matrix must have at least two columns (query, subject)",
           call. = FALSE)
    }
    return(data.frame(
      q = as.character(pairs[, 1]),
      s = as.character(pairs[, 2]),
      stringsAsFactors = FALSE
    ))
  }
  if (is.list(pairs)) {
    if (length(pairs) == 2 && all(vapply(pairs, function(x) {
      is.character(x) && length(x) == 1
    }, logical(1)))) {
      return(data.frame(q = pairs[[1]], s = pairs[[2]],
                        stringsAsFactors = FALSE))
    }
    ok <- all(vapply(pairs, function(x) {
      is.character(x) && length(x) == 2
    }, logical(1)))
    if (!ok) {
      ggchord_stop("keep_pairs must be a data.frame/matrix with query/subject columns, a list of length-2 character vectors, or a single length-2 character vector",
           call. = FALSE)
    }
    m <- do.call(rbind, pairs)
    return(data.frame(q = m[, 1], s = m[, 2], stringsAsFactors = FALSE))
  }
  if (is.character(pairs) && length(pairs) == 2) {
    return(data.frame(q = pairs[1], s = pairs[2], stringsAsFactors = FALSE))
  }
  ggchord_stop("keep_pairs must be a data.frame/matrix with query/subject columns, a list of length-2 character vectors, or a single length-2 character vector",
       call. = FALSE)
}

#' Deduplicate alignment ribbons
#'
#' Removes fully duplicated, coordinate-near-duplicated or highly overlapping
#' alignment blocks within each (query, subject) pair and keeps the best
#' representative per \code{keep}.
#'
#' @param ribbon_data data.frame in ribbon_data format.
#' @param tolerance Numeric, default 0. Maximum absolute difference (in bp)
#'   allowed on each of \code{qstart}/\code{qend}/\code{sstart}/\code{send} for
#'   two rows to count as coordinate duplicates (used with
#'   \code{by = "coordinates"}).
#' @param by Character, default \code{"exact"}. \code{"exact"} removes rows
#'   with identical pair and coordinates; \code{"coordinates"} additionally
#'   treats rows within \code{tolerance} bp as duplicates;
#'   \code{"overlap"} treats blocks whose reciprocal overlap is at least
#'   \code{min_reciprocal_overlap} on both sequences as duplicates.
#' @param keep Character, default \code{"best_pident"}. Which representative to
#'   keep: \code{"best_pident"} (highest pident), \code{"longest"} (longest
#'   length) or \code{"first"} (first occurrence in the input).
#' @param min_reciprocal_overlap Numeric (0-1), default 0.9. Reciprocal overlap
#'   threshold used with \code{by = "overlap"}.
#'
#' @return A list with \code{data} (deduplicated data frame with
#'   \code{source_rows} attribute) and \code{report} (n_input/n_kept/n_removed
#'   plus a data.frame of removed rows with \code{row}, \code{duplicate_of} and
#'   \code{reason}).
#' @export
#'
#' @examples
#' library(ggchord)
#' data(ribbon_data_example)
#' dup <- rbind(ribbon_data_example, ribbon_data_example[1, ])
#' out <- deduplicate_ggchord_ribbons(dup, by = "exact")
#' out$report
deduplicate_ggchord_ribbons <- function(
    ribbon_data,
    tolerance = 0,
    by = c("exact", "coordinates", "overlap"),
    keep = c("best_pident", "longest", "first"),
    min_reciprocal_overlap = 0.9) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  by <- match.arg(by)
  keep <- match.arg(keep)
  if (!is.data.frame(ribbon_data)) {
    ggchord_stop("deduplicate_ggchord_ribbons(): ribbon_data must be a data.frame",
         call. = FALSE)
  }
  req <- c("qaccver", "saccver", "length", "pident",
           "qstart", "qend", "sstart", "send")
  missing <- setdiff(req, colnames(ribbon_data))
  if (length(missing) > 0) {
    ggchord_stop("deduplicate_ggchord_ribbons(): ribbon_data is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1 || is.na(tolerance) ||
      tolerance < 0) {
    ggchord_stop("deduplicate_ggchord_ribbons(): tolerance must be a non-negative number",
         call. = FALSE)
  }
  if (!is.numeric(min_reciprocal_overlap) || length(min_reciprocal_overlap) != 1 ||
      is.na(min_reciprocal_overlap) || min_reciprocal_overlap < 0 ||
      min_reciprocal_overlap > 1) {
    ggchord_stop("deduplicate_ggchord_ribbons(): min_reciprocal_overlap must be in [0, 1]",
         call. = FALSE)
  }

  n <- nrow(ribbon_data)
  if (n == 0) {
    out <- ribbon_data
    attr(out, "source_rows") <- integer(0)
    return(list(
      data = out,
      report = list(
        n_input = 0L, n_kept = 0L, n_removed = 0L,
        removed = data.frame(row = integer(0), duplicate_of = integer(0),
                             reason = character(0))
      )
    ))
  }

  pair <- paste(ribbon_data$qaccver, ribbon_data$saccver, sep = "\r")
  grp <- split(seq_len(n), pair)
  keep_idx <- logical(n)
  removed <- data.frame(row = integer(0), duplicate_of = integer(0),
                        reason = character(0), stringsAsFactors = FALSE)

  is_dup <- function(a, b) {
    if (by == "exact") {
      return(identical(
        unname(unlist(ribbon_data[a, c("qstart", "qend", "sstart", "send")])),
        unname(unlist(ribbon_data[b, c("qstart", "qend", "sstart", "send")]))))
    }
    if (by == "coordinates") {
      d <- max(abs(ribbon_data$qstart[a] - ribbon_data$qstart[b]),
               abs(ribbon_data$qend[a] - ribbon_data$qend[b]),
               abs(ribbon_data$sstart[a] - ribbon_data$sstart[b]),
               abs(ribbon_data$send[a] - ribbon_data$send[b]))
      return(d <= tolerance)
    }
    qr <- interval_recip_overlap(ribbon_data$qstart[a], ribbon_data$qend[a],
                                 ribbon_data$qstart[b], ribbon_data$qend[b])
    sr <- interval_recip_overlap(ribbon_data$sstart[a], ribbon_data$send[a],
                                 ribbon_data$sstart[b], ribbon_data$send[b])
    is.finite(qr) && is.finite(sr) &&
      qr >= min_reciprocal_overlap && sr >= min_reciprocal_overlap
  }
  better <- function(a, b) {
    if (keep == "best_pident") {
      return(ribbon_data$pident[a] > ribbon_data$pident[b])
    }
    if (keep == "longest") {
      return(ribbon_data$length[a] > ribbon_data$length[b])
    }
    FALSE  # "first": never replace the earlier representative
  }

  for (g in grp) {
    if (length(g) <= 1) {
      keep_idx[g] <- TRUE
      next
    }
    ord_g <- g[order(ribbon_data$qstart[g], ribbon_data$sstart[g])]
    for (i in seq_along(ord_g)) {
      cur <- ord_g[i]
      if (keep_idx[cur]) next
      dup_of <- NA_integer_
      for (j in seq_len(i - 1L)) {
        prev <- ord_g[j]
        if (!keep_idx[prev]) next
        if (is_dup(prev, cur)) {
          dup_of <- prev
          break
        }
      }
      if (is.na(dup_of)) {
        keep_idx[cur] <- TRUE
      } else if (better(cur, dup_of)) {
        keep_idx[dup_of] <- FALSE
        removed <- rbind(removed, data.frame(
          row = dup_of, duplicate_of = cur,
          reason = paste0("duplicate_", by), stringsAsFactors = FALSE))
        keep_idx[cur] <- TRUE
      } else {
        removed <- rbind(removed, data.frame(
          row = cur, duplicate_of = dup_of,
          reason = paste0("duplicate_", by), stringsAsFactors = FALSE))
      }
    }
  }

  out <- ribbon_data[keep_idx, , drop = FALSE]
  attr(out, "source_rows") <- which(keep_idx)
  list(
    data = out,
    report = list(
      n_input = n,
      n_kept = sum(keep_idx),
      n_removed = sum(!keep_idx),
      removed = removed
    )
  )
}

#' Merge adjacent or overlapping alignment blocks of the same sequence pair
#'
#' Merges alignment blocks that belong to the same (query, subject) pair, are
#' adjacent (gap <= \code{max_gap}) or overlapping on both sequences, and are
#' compatible. The merged pident is weighted by alignment length. Merging is
#' deliberately conservative: blocks whose merged query and subject spans would
#' be inconsistent (unequal), whose pident differs by more than
#' \code{min_pident_difference}, or whose orientation differs (when
#' \code{require_same_orientation}) are left unmerged.
#'
#' @param ribbon_data data.frame in ribbon_data format.
#' @param max_gap Numeric, default 0. Maximum gap (in bp) allowed between two
#'   blocks on both sequences for them to be merged.
#' @param min_pident_difference Numeric, default 0. When > 0, two blocks are
#'   only merged when their pident difference is <= this value.
#' @param require_same_orientation Logical, default \code{TRUE}. Only merge
#'   blocks whose query/subject direction is the same (both ascending or both
#'   descending, i.e. collinear or both inverted in the same way).
#' @param group_by Character vector, default \code{c("qaccver", "saccver")}.
#'   Columns used to identify the same sequence pair.
#'
#' @return A list with \code{data} (merged data frame with \code{source_rows}
#'   attribute) and \code{report} (data.frame with \code{output_row},
#'   \code{from_rows} and \code{n_merged} for every output row).
#' @export
#'
#' @examples
#' library(ggchord)
#' # Two adjacent, collinear blocks of the same pair -> one merged block
#' rb <- data.frame(
#'   qaccver = c("A", "A"), saccver = c("B", "B"),
#'   length = c(100, 100), pident = c(95, 97),
#'   qstart = c(1, 101), qend = c(100, 200),
#'   sstart = c(501, 601), send = c(600, 700)
#' )
#' out <- merge_ggchord_ribbons(rb, max_gap = 0)
#' out$data
#' out$report
merge_ggchord_ribbons <- function(
    ribbon_data,
    max_gap = 0,
    min_pident_difference = 0,
    require_same_orientation = TRUE,
    group_by = c("qaccver", "saccver")) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!is.data.frame(ribbon_data)) {
    ggchord_stop("merge_ggchord_ribbons(): ribbon_data must be a data.frame",
         call. = FALSE)
  }
  req <- c("qaccver", "saccver", "length", "pident",
           "qstart", "qend", "sstart", "send")
  missing <- setdiff(req, colnames(ribbon_data))
  if (length(missing) > 0) {
    ggchord_stop("merge_ggchord_ribbons(): ribbon_data is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(group_by) == 0 || !all(group_by %in% colnames(ribbon_data))) {
    ggchord_stop("merge_ggchord_ribbons(): group_by must be a subset of the ribbon_data columns",
         call. = FALSE)
  }
  if (!is.numeric(max_gap) || length(max_gap) != 1 || is.na(max_gap) ||
      max_gap < 0) {
    ggchord_stop("merge_ggchord_ribbons(): max_gap must be a non-negative number",
         call. = FALSE)
  }
  if (!is.numeric(min_pident_difference) || length(min_pident_difference) != 1 ||
      is.na(min_pident_difference) || min_pident_difference < 0) {
    ggchord_stop("merge_ggchord_ribbons(): min_pident_difference must be a non-negative number",
         call. = FALSE)
  }

  n <- nrow(ribbon_data)
  if (n == 0) {
    out <- ribbon_data
    attr(out, "source_rows") <- integer(0)
    return(list(
      data = out,
      report = data.frame(output_row = integer(0),
                          from_rows = character(0),
                          n_merged = integer(0))
    ))
  }

  qA <- pmin(ribbon_data$qstart, ribbon_data$qend)
  qB <- pmax(ribbon_data$qstart, ribbon_data$qend)
  sA <- pmin(ribbon_data$sstart, ribbon_data$send)
  sB <- pmax(ribbon_data$sstart, ribbon_data$send)
  qdir <- sign(ribbon_data$qend - ribbon_data$qstart)
  sdir <- sign(ribbon_data$send - ribbon_data$sstart)
  orientation <- ifelse(qdir == sdir, "same", "opposite")

  interval_gap <- function(l1, h1, l2, h2) {
    if (h1 < l2) return(l2 - h1 - 1)
    if (h2 < l1) return(l1 - h2 - 1)
    0
  }

  grp_key <- do.call(paste, c(ribbon_data[group_by], sep = "\r"))
  grp <- split(seq_len(n), grp_key)

  out_chunks <- list()
  out_src <- list()
  report <- data.frame(output_row = integer(0), from_rows = character(0),
                       n_merged = integer(0), stringsAsFactors = FALSE)

  for (g in grp) {
    ord <- g[order(qA[g], sA[g])]
    if (isTRUE(require_same_orientation)) {
      sub_groups <- split(ord, orientation[ord])
    } else {
      sub_groups <- list(ord)
    }
    for (sg in sub_groups) {
      if (length(sg) == 0) next
      cur_rows <- sg[1]
      cur_qA <- qA[sg[1]]
      cur_qB <- qB[sg[1]]
      cur_sA <- sA[sg[1]]
      cur_sB <- sB[sg[1]]
      cur_qdir <- qdir[sg[1]]
      cur_sdir <- sdir[sg[1]]
      cur_wsum <- ribbon_data$pident[sg[1]] * ribbon_data$length[sg[1]]
      cur_lsum <- ribbon_data$length[sg[1]]
      cur_pident <- ribbon_data$pident[sg[1]]

      emit <- function() {
        r0 <- cur_rows[1]
        new_row <- ribbon_data[r0, , drop = FALSE]
        new_row$qstart <- if (cur_qdir >= 0) cur_qA else cur_qB
        new_row$qend <- if (cur_qdir >= 0) cur_qB else cur_qA
        new_row$sstart <- if (cur_sdir >= 0) cur_sA else cur_sB
        new_row$send <- if (cur_sdir >= 0) cur_sB else cur_sA
        new_row$length <- cur_qB - cur_qA + 1
        new_row$pident <- cur_wsum / cur_lsum
        out_chunks[[length(out_chunks) + 1L]] <<- new_row
        out_src[[length(out_src) + 1L]] <<- cur_rows
        report <<- rbind(report, data.frame(
          output_row = length(out_chunks),
          from_rows = paste(cur_rows, collapse = ","),
          n_merged = length(cur_rows),
          stringsAsFactors = FALSE))
      }

      if (length(sg) > 1) {
        for (i in seq.int(2L, length(sg))) {
          r <- sg[i]
          q_gap <- interval_gap(cur_qA, cur_qB, qA[r], qB[r])
          s_gap <- interval_gap(cur_sA, cur_sB, sA[r], sB[r])
          pid_ok <- isTRUE(min_pident_difference == 0) ||
            abs(cur_pident - ribbon_data$pident[r]) <= min_pident_difference
          if (q_gap <= max_gap && s_gap <= max_gap && pid_ok) {
            q_span <- max(cur_qB, qB[r]) - min(cur_qA, qA[r]) + 1
            s_span <- max(cur_sB, sB[r]) - min(cur_sA, sA[r]) + 1
            if (abs(q_span - s_span) <= max(max_gap, 1)) {
              cur_qA <- min(cur_qA, qA[r])
              cur_qB <- max(cur_qB, qB[r])
              cur_sA <- min(cur_sA, sA[r])
              cur_sB <- max(cur_sB, sB[r])
              cur_wsum <- cur_wsum + ribbon_data$pident[r] * ribbon_data$length[r]
              cur_lsum <- cur_lsum + ribbon_data$length[r]
              cur_pident <- cur_wsum / cur_lsum
              cur_rows <- c(cur_rows, r)
              next
            }
          }
          emit()
          cur_rows <- r
          cur_qA <- qA[r]
          cur_qB <- qB[r]
          cur_sA <- sA[r]
          cur_sB <- sB[r]
          cur_qdir <- qdir[r]
          cur_sdir <- sdir[r]
          cur_wsum <- ribbon_data$pident[r] * ribbon_data$length[r]
          cur_lsum <- ribbon_data$length[r]
          cur_pident <- ribbon_data$pident[r]
        }
      }
      emit()
    }
  }

  out <- do.call(rbind, out_chunks)
  attr(out, "source_rows") <- unlist(out_src)
  list(data = out, report = report)
}
