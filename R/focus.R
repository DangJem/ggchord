# Synchronized locus focusing for sequence, ribbon and annotation tables.

ggchord_focus_validate_loci <- function(seq_data, loci, expand, caller) {
  if (!is.data.frame(seq_data) ||
      !all(c("accver", "length") %in% names(seq_data))) {
    ggchord_stop(caller, ": seq_data must contain accver and length")
  }
  if (!is.data.frame(loci) ||
      !all(c("accver", "start", "end") %in% names(loci))) {
    ggchord_stop(caller, ": loci must contain accver, start and end")
  }
  if (nrow(loci) == 0L) ggchord_stop(caller, ": loci must not be empty")
  seq_ids <- as.character(seq_data$accver)
  if (anyNA(seq_ids) || any(!nzchar(seq_ids)) || anyDuplicated(seq_ids) ||
      !is.numeric(seq_data$length) || any(!is.finite(seq_data$length)) ||
      any(seq_data$length <= 0)) {
    ggchord_stop(caller, ": seq_data must contain unique IDs and positive lengths")
  }
  if (!is.numeric(loci$start) || !is.numeric(loci$end) ||
      any(!is.finite(loci$start)) || any(!is.finite(loci$end))) {
    ggchord_stop(caller, ": loci start and end must be finite numbers")
  }
  if (!is.numeric(expand) || anyNA(expand) || any(!is.finite(expand)) ||
      any(expand < 0) || !length(expand) ||
      !(length(expand) %in% c(1L, nrow(loci)))) {
    ggchord_stop(caller, ": expand must be one non-negative number or one per locus")
  }
  if (length(expand) == 1L) expand <- rep(expand, nrow(loci))

  out <- as.data.frame(loci, stringsAsFactors = FALSE)
  out$accver <- as.character(out$accver)
  unknown <- setdiff(unique(out$accver), seq_ids)
  if (length(unknown)) {
    ggchord_stop(caller, ": loci contain unknown sequence ID(s): ",
                 paste(unknown, collapse = ", "))
  }
  lens <- stats::setNames(as.numeric(seq_data$length), seq_ids)
  out$.source_accver <- out$accver
  out$.focus_start <- pmax(1, pmin(out$start, out$end) - expand)
  out$.focus_end <- pmin(lens[out$accver], pmax(out$start, out$end) + expand)
  if (any(out$.focus_end < out$.focus_start)) {
    ggchord_stop(caller, ": every expanded locus must overlap its sequence")
  }
  out$.locus_row <- seq_len(nrow(out))

  if ("locus_id" %in% names(out)) {
    output_id <- as.character(out$locus_id)
    if (anyNA(output_id) || any(!nzchar(output_id))) {
      ggchord_stop(caller, ": locus_id must be non-missing and non-empty")
    }
  } else {
    repeated <- duplicated(out$accver) | duplicated(out$accver, fromLast = TRUE)
    output_id <- out$accver
    output_id[repeated] <- paste0(
      out$accver[repeated], ":", out$.focus_start[repeated], "-",
      out$.focus_end[repeated]
    )
  }
  if (anyDuplicated(output_id)) {
    ggchord_stop(caller, ": focused sequence IDs must be unique; supply unique locus_id values")
  }
  out$.output_accver <- output_id
  out
}

ggchord_focus_interval <- function(a, b, lo, hi) {
  if (a == b) {
    if (a >= lo && a <= hi) return(c(0, 1))
    return(NULL)
  }
  t <- sort(c((lo - a) / (b - a), (hi - a) / (b - a)))
  range <- c(max(0, t[1]), min(1, t[2]))
  if (range[2] < range[1]) NULL else range
}

#' Focus ggchord data on one or more loci
#'
#' Crops sequence, gene and ribbon tables together and relocates retained
#' coordinates to start at one in each output locus. Multiple loci from the
#' same source sequence are supported; provide a unique `locus_id` column to
#' control their output sequence IDs.
#'
#' @param seq_data Sequence table containing `accver` and `length`.
#' @param ribbon_data Optional ggchord alignment table.
#' @param gene_data Optional gene or feature table.
#' @param loci Data frame containing `accver`, `start`, and `end`, and
#'   optionally `locus_id`.
#' @param expand Non-negative number of bases added on both sides, scalar or
#'   one value per locus.
#' @param boundary Boundary policy. `"trim"` clips partially overlapping genes
#'   and alignments; `"drop"` keeps only rows fully contained in selected loci.
#' @return A list of focused `seq_data`, `ribbon_data`, `gene_data`, normalized
#'   `loci`, and a compact `report`.
#' @export
focus_ggchord_data <- function(
    seq_data,
    ribbon_data = NULL,
    gene_data = NULL,
    loci,
    expand = 0,
    boundary = c("trim", "drop")) {
  seq_data <- ggchord_normalize_accver(seq_data)
  gene_data <- ggchord_normalize_accver(gene_data)
  loci <- ggchord_normalize_accver(loci)

  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  caller <- "focus_ggchord_data()"
  boundary <- match.arg(boundary)
  focus <- ggchord_focus_validate_loci(seq_data, loci, expand, caller)

  seq_index <- match(focus$.source_accver, as.character(seq_data$accver))
  seq_out <- as.data.frame(seq_data[seq_index, , drop = FALSE],
                           stringsAsFactors = FALSE)
  seq_out$.source_accver <- focus$.source_accver
  seq_out$.focus_start <- focus$.focus_start
  seq_out$.focus_end <- focus$.focus_end
  seq_out$.locus_row <- focus$.locus_row
  seq_out$accver <- focus$.output_accver
  seq_out$length <- focus$.focus_end - focus$.focus_start + 1

  gene_out <- NULL
  if (!is.null(gene_data)) {
    required <- c("accver", "start", "end", "strand")
    if (!is.data.frame(gene_data) || any(!required %in% names(gene_data)) ||
        !is.numeric(gene_data$start) || !is.numeric(gene_data$end)) {
      ggchord_stop(caller, ": gene_data must contain accver, numeric start/end, and strand")
    }
    pieces <- list()
    for (i in seq_len(nrow(gene_data))) {
      hits <- which(focus$.source_accver == as.character(gene_data$accver[i]))
      if (!length(hits)) next
      a <- gene_data$start[i]
      b <- gene_data$end[i]
      glo <- min(a, b)
      ghi <- max(a, b)
      for (j in hits) {
        lo <- focus$.focus_start[j]
        hi <- focus$.focus_end[j]
        keep <- if (boundary == "drop") glo >= lo && ghi <= hi else {
          ghi >= lo && glo <= hi
        }
        if (!keep) next
        clipped <- c(max(glo, lo), min(ghi, hi))
        row <- as.data.frame(gene_data[i, , drop = FALSE],
                             stringsAsFactors = FALSE)
        row$accver <- focus$.output_accver[j]
        if (a <= b) {
          row$start <- clipped[1] - lo + 1
          row$end <- clipped[2] - lo + 1
        } else {
          row$start <- clipped[2] - lo + 1
          row$end <- clipped[1] - lo + 1
        }
        row$.source_accver <- focus$.source_accver[j]
        row$.source_row <- i
        row$.locus_row <- focus$.locus_row[j]
        pieces[[length(pieces) + 1L]] <- row
      }
    }
    gene_out <- if (length(pieces)) ggchord_rbind_fill(pieces) else {
      empty <- gene_data[0, , drop = FALSE]
      empty$.source_accver <- character(0)
      empty$.source_row <- integer(0)
      empty$.locus_row <- integer(0)
      empty
    }
  }

  ribbon_out <- NULL
  if (!is.null(ribbon_data)) {
    ggchord_check_ribbon_tables(seq_data, ribbon_data, caller)
    pieces <- list()
    for (i in seq_len(nrow(ribbon_data))) {
      qhits <- which(focus$.source_accver == as.character(ribbon_data$qaccver[i]))
      shits <- which(focus$.source_accver == as.character(ribbon_data$saccver[i]))
      if (!length(qhits) || !length(shits)) next
      for (qj in qhits) for (sj in shits) {
        qr <- ggchord_focus_interval(
          ribbon_data$qstart[i], ribbon_data$qend[i],
          focus$.focus_start[qj], focus$.focus_end[qj]
        )
        sr <- ggchord_focus_interval(
          ribbon_data$sstart[i], ribbon_data$send[i],
          focus$.focus_start[sj], focus$.focus_end[sj]
        )
        if (is.null(qr) || is.null(sr)) next
        tr <- if (boundary == "drop") {
          if (qr[1] > 0 || qr[2] < 1 || sr[1] > 0 || sr[2] < 1) next
          c(0, 1)
        } else {
          c(max(qr[1], sr[1]), min(qr[2], sr[2]))
        }
        if (tr[2] < tr[1]) next
        row <- as.data.frame(ribbon_data[i, , drop = FALSE],
                             stringsAsFactors = FALSE)
        qcoord <- ribbon_data$qstart[i] + tr *
          (ribbon_data$qend[i] - ribbon_data$qstart[i])
        scoord <- ribbon_data$sstart[i] + tr *
          (ribbon_data$send[i] - ribbon_data$sstart[i])
        row$qaccver <- focus$.output_accver[qj]
        row$saccver <- focus$.output_accver[sj]
        row$qstart <- round(qcoord[1] - focus$.focus_start[qj] + 1)
        row$qend <- round(qcoord[2] - focus$.focus_start[qj] + 1)
        row$sstart <- round(scoord[1] - focus$.focus_start[sj] + 1)
        row$send <- round(scoord[2] - focus$.focus_start[sj] + 1)
        if ("length" %in% names(ribbon_data)) row$length <- max(1, round(ribbon_data$length[i] * diff(tr)))
        row$.source_row <- i
        row$.query_locus_row <- focus$.locus_row[qj]
        row$.subject_locus_row <- focus$.locus_row[sj]
        pieces[[length(pieces) + 1L]] <- row
      }
    }
    ribbon_out <- if (length(pieces)) ggchord_rbind_fill(pieces) else {
      empty <- ribbon_data[0, , drop = FALSE]
      empty$.source_row <- integer(0)
      empty$.query_locus_row <- integer(0)
      empty$.subject_locus_row <- integer(0)
      empty
    }
  }

  list(
    seq_data = seq_out,
    ribbon_data = ribbon_out,
    gene_data = gene_out,
    loci = focus,
    report = data.frame(
      component = c("sequence", "ribbon", "gene"),
      input_rows = c(nrow(seq_data), if (is.null(ribbon_data)) 0L else nrow(ribbon_data),
                     if (is.null(gene_data)) 0L else nrow(gene_data)),
      output_rows = c(nrow(seq_out), if (is.null(ribbon_out)) 0L else nrow(ribbon_out),
                      if (is.null(gene_out)) 0L else nrow(gene_out)),
      boundary = boundary,
      stringsAsFactors = FALSE
    )
  )
}
