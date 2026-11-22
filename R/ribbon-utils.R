# ribbon-utils.R - ribbon filtering, deduplication and merging
#
# v0.8.0: ggchord_filter_link_rows(), ggchord_deduplicate_link_rows() and
# merge_ggchord_ribbons() prepare alignment tables before plotting. All three
# keep the original column order plus extra columns, return a report, and
# attach the original row numbers as the "source_rows" attribute of the data.

ggchord_ribbon_numeric_input <- function(data, columns, caller) {
  for (column in columns) {
    if (!is.numeric(data[[column]]) || any(!is.finite(data[[column]])))
      ggchord_stop(caller, ": ", column, " must contain finite numeric values")
  }
  for (column in intersect(c("qaccver", "saccver"), names(data)))
    ggchord_conversion_ids(data[[column]], caller)
  invisible(NULL)
}

ggchord_filter_link_rows <- function(
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
  if (!is.data.frame(ribbon_data)) {
    ggchord_stop("filter_ggchord_links(): ribbon_data must be a data.frame",
         call. = FALSE)
  }
  req <- c("qaccver", "saccver")
  missing <- setdiff(req, colnames(ribbon_data))
  if (length(missing) > 0) {
    ggchord_stop("filter_ggchord_links(): ribbon_data is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  thresholds <- list(min_pident = min_pident, max_pident = max_pident,
    min_length = min_length, max_evalue = max_evalue, min_bitscore = min_bitscore,
    min_query_coverage = min_query_coverage, min_subject_coverage = min_subject_coverage)
  for (name in names(thresholds)) {
    value <- thresholds[[name]]
    if (!is.null(value) && (!is.numeric(value) || length(value) != 1L || !is.finite(value)))
      ggchord_stop("filter_ggchord_links(): ", name, " must be NULL or one finite number")
  }
  if (!is.logical(drop_self_links) || length(drop_self_links) != 1L || is.na(drop_self_links))
    ggchord_stop("filter_ggchord_links(): drop_self_links must be TRUE or FALSE")
  n_input <- nrow(ribbon_data)
  reasons <- replicate(n_input, character(0), simplify = FALSE)
  mark_reason <- function(idx, value) {
    for (i in idx) reasons[[i]] <<- unique(c(reasons[[i]], value))
  }

  require_col <- function(col, arg) {
    if (!col %in% colnames(ribbon_data)) {
      ggchord_stop(sprintf(paste0(
        "filter_ggchord_links(): column '%s' is required for the '%s' ",
        "filter but is missing from ribbon_data"), col, arg), call. = FALSE)
    }
    if (!is.numeric(ribbon_data[[col]])) ggchord_stop(
      "filter_ggchord_links(): ", col, " must be numeric for ", arg)
  }

  if (!is.null(seq_ids)) {
    idx <- which(!(ribbon_data$qaccver %in% seq_ids &
                     ribbon_data$saccver %in% seq_ids))
    if (length(idx) > 0) mark_reason(idx, "seq_ids")
  }
  if (!is.null(min_pident)) {
    require_col("pident", "min_pident")
    idx <- which(is.na(ribbon_data$pident) | ribbon_data$pident < min_pident)
    if (length(idx) > 0) mark_reason(idx, "min_pident")
  }
  if (!is.null(max_pident)) {
    require_col("pident", "max_pident")
    idx <- which(is.na(ribbon_data$pident) | ribbon_data$pident > max_pident)
    if (length(idx) > 0) mark_reason(idx, "max_pident")
  }
  if (!is.null(min_length)) {
    require_col("length", "min_length")
    idx <- which(is.na(ribbon_data$length) | ribbon_data$length < min_length)
    if (length(idx) > 0) mark_reason(idx, "min_length")
  }
  if (!is.null(max_evalue)) {
    require_col("evalue", "max_evalue")
    idx <- which(is.na(ribbon_data$evalue) | ribbon_data$evalue > max_evalue)
    if (length(idx) > 0) mark_reason(idx, "max_evalue")
  }
  if (!is.null(min_bitscore)) {
    require_col("bitscore", "min_bitscore")
    idx <- which(is.na(ribbon_data$bitscore) | ribbon_data$bitscore < min_bitscore)
    if (length(idx) > 0) mark_reason(idx, "min_bitscore")
  }
  if (!is.null(min_query_coverage)) {
    if ("qcovs" %in% colnames(ribbon_data)) {
      require_col("qcovs", "min_query_coverage")
      cov <- ribbon_data$qcovs
    } else {
      require_col("qlen", "min_query_coverage")
      require_col("qstart", "min_query_coverage")
      require_col("qend", "min_query_coverage")
      cov <- ifelse(
        is.na(ribbon_data$qlen) | ribbon_data$qlen <= 0, NA_real_,
        (abs(ribbon_data$qend - ribbon_data$qstart) + 1) /
          ribbon_data$qlen * 100
      )
    }
    idx <- which(is.na(cov) | cov < min_query_coverage)
    if (length(idx) > 0) mark_reason(idx, "min_query_coverage")
  }
  if (!is.null(min_subject_coverage)) {
    require_col("slen", "min_subject_coverage")
    require_col("sstart", "min_subject_coverage")
    require_col("send", "min_subject_coverage")
    cov <- ifelse(
      is.na(ribbon_data$slen) | ribbon_data$slen <= 0, NA_real_,
      (abs(ribbon_data$send - ribbon_data$sstart) + 1) /
        ribbon_data$slen * 100
    )
    idx <- which(is.na(cov) | cov < min_subject_coverage)
    if (length(idx) > 0) mark_reason(idx, "min_subject_coverage")
  }
  if (isTRUE(drop_self_links)) {
    idx <- which(ribbon_data$qaccver == ribbon_data$saccver)
    if (length(idx) > 0) mark_reason(idx, "self_link")
  }
  if (!is.null(keep_pairs)) {
    pairs <- normalize_keep_pairs(keep_pairs)
    ggchord_conversion_ids(pairs$q, "ggchord_filter_link_rows() keep_pairs")
    ggchord_conversion_ids(pairs$s, "ggchord_filter_link_rows() keep_pairs")
    keep_row <- vapply(seq_len(n_input), function(i) {
      q <- as.character(ribbon_data$qaccver[i])
      s <- as.character(ribbon_data$saccver[i])
      any((pairs$q == q & pairs$s == s) | (pairs$q == s & pairs$s == q))
    }, logical(1))
    idx <- which(!keep_row)
    if (length(idx) > 0) mark_reason(idx, "keep_pairs")
  }

  kept <- lengths(reasons) == 0L
  removed_rows <- which(!kept)
  out <- ribbon_data[kept, , drop = FALSE]
  kept_src <- which(kept)

  sort_info <- NULL
  if (!is.null(sort_by)) {
    if (!is.character(sort_by) || length(sort_by) == 0 || anyNA(sort_by) || any(!nzchar(sort_by))) {
      ggchord_stop("filter_ggchord_links(): sort_by must be a non-empty character vector",
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
        ggchord_stop("filter_ggchord_links(): sort_by column not found: ", col,
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

  removed_detail <- if (length(removed_rows) == 0) {
    data.frame(row = integer(0), reason = character(0))
  } else {
    do.call(rbind, lapply(removed_rows, function(i) {
      data.frame(row = i, reason = reasons[[i]], stringsAsFactors = FALSE)
    }))
  }
  if (nrow(removed_detail) == 0) {
    removed_by_reason <- data.frame(reason = character(0), n = integer(0))
  } else {
    removed_by_reason <- as.data.frame(table(removed_detail$reason),
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
      removed_reasons = removed_detail,
      kept_rows = which(kept),
      sort_by = sort_info
    )
  )
}

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

ggchord_deduplicate_link_rows <- function(
    ribbon_data,
    tolerance = 0,
    by = c("exact", "coordinates", "overlap"),
    keep = c("best_pident", "longest", "first"),
    min_reciprocal_overlap = 0.9) {
  by <- match.arg(by)
  keep <- match.arg(keep)
  if (!is.data.frame(ribbon_data)) {
    ggchord_stop("deduplicate_ggchord_links(): ribbon_data must be a data.frame",
         call. = FALSE)
  }
  req <- c(ggchord_ribbon_required_columns(),
    switch(keep, first = character(), longest = "length", best_pident = "pident"))
  missing <- setdiff(req, colnames(ribbon_data))
  if (length(missing) > 0) {
    ggchord_stop("deduplicate_ggchord_links(): ribbon_data is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  ggchord_ribbon_numeric_input(ribbon_data, setdiff(req, c("qaccver", "saccver")),
    "ggchord_deduplicate_link_rows()")
  if (!is.numeric(tolerance) || length(tolerance) != 1 || is.na(tolerance) ||
      !is.finite(tolerance) || tolerance < 0) {
    ggchord_stop("deduplicate_ggchord_links(): tolerance must be a non-negative number",
         call. = FALSE)
  }
  if (!is.numeric(min_reciprocal_overlap) || length(min_reciprocal_overlap) != 1 ||
      is.na(min_reciprocal_overlap) || min_reciprocal_overlap < 0 ||
      min_reciprocal_overlap > 1) {
    ggchord_stop("deduplicate_ggchord_links(): min_reciprocal_overlap must be in [0, 1]",
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
  relative_direction <- ggchord_link_orientation(ribbon_data)
  keep_idx <- logical(n)
  removed <- data.frame(row = integer(0), duplicate_of = integer(0),
                        reason = character(0), stringsAsFactors = FALSE)

  is_dup <- function(a, b) {
    if (relative_direction[a] != relative_direction[b]) return(FALSE)
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
    # "first" means first in the user's input, not first after coordinate
    # sorting. Other policies retain coordinate order to avoid needless
    # comparisons between distant blocks.
    ord_g <- if (keep == "first") {
      g
    } else {
      g[order(ribbon_data$qstart[g], ribbon_data$sstart[g])]
    }
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
        # Earlier duplicates must point at the final representative rather
        # than at the representative that has just been replaced.
        removed$duplicate_of[removed$duplicate_of == dup_of] <- cur
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
