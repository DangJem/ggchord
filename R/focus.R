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

# Focus transforms genomic coordinates, including nested annotation tables.
# Provenance is kept separately so repeated focusing still maps back to the
# original sequence rather than treating the previous window as a new genome.
ggchord_focus_source <- function(seq_data) {
  number <- function(name, fallback) {
    value <- seq_data[[name]]
    if (is.null(value)) return(fallback)
    value <- as.numeric(value)
    ifelse(is.finite(value), value, fallback)
  }
  origin <- number(".source_origin",
    number(".source_start", number(".focus_start", rep(1, nrow(seq_data)))))
  source_accver <- as.character(seq_data$.source_accver %||% seq_data$accver)
  source_length <- number(".source_length", if (".focus_start" %in% names(seq_data)) {
    rep(NA_real_, nrow(seq_data))
  } else as.numeric(seq_data$length))
  list(origin = origin, accver = source_accver, length = source_length)
}

ggchord_focus_original <- function(row, name, fallback) {
  value <- row[[name]]
  if (is.numeric(value) && length(value) == 1L && is.finite(value)) value else fallback
}

ggchord_focus_points <- function(value, lo, hi) {
  value <- suppressWarnings(as.numeric(value))
  value <- value[is.finite(value) & value >= lo & value <= hi]
  value - lo + 1
}

ggchord_focus_nested_intervals <- function(data, lo, hi, origin, caller,
                                           point_column = NULL) {
  if (is.null(data) || (is.data.frame(data) && !nrow(data))) return(data)
  if (!is.data.frame(data) || !all(c("start", "end") %in% names(data)) ||
      !is.numeric(data$start) || !is.numeric(data$end) ||
      any(!is.finite(data$start)) || any(!is.finite(data$end))) {
    ggchord_stop(caller, ": nested feature intervals must contain finite numeric start/end")
  }
  a <- data$start; b <- data$end
  lower <- pmin(a, b); upper <- pmax(a, b)
  keep <- upper >= lo & lower <= hi
  out <- as.data.frame(data[keep, , drop = FALSE], stringsAsFactors = FALSE)
  a <- a[keep]; b <- b[keep]
  lower <- lower[keep]; upper <- upper[keep]
  clipped_lower <- pmax(lower, lo); clipped_upper <- pmin(upper, hi)
  if (!".source_start" %in% names(out)) out$.source_start <- a + origin - 1
  if (!".source_end" %in% names(out)) out$.source_end <- b + origin - 1
  out$.truncated_left <- (out$.truncated_left %||% rep(FALSE, nrow(out))) |
    clipped_lower > lower
  out$.truncated_right <- (out$.truncated_right %||% rep(FALSE, nrow(out))) |
    clipped_upper < upper
  out$start <- ifelse(a <= b, clipped_lower, clipped_upper) - lo + 1
  out$end <- ifelse(a <= b, clipped_upper, clipped_lower) - lo + 1
  if ("length_bp" %in% names(out)) {
    if (!".source_length_bp" %in% names(out)) out$.source_length_bp <- out$length_bp
    out$length_bp <- clipped_upper - clipped_lower + 1
  }
  if (!is.null(point_column) && point_column %in% names(out)) {
    for (i in seq_len(nrow(out))) {
      points <- out[[point_column]][[i]]
      inside <- is.finite(points) & points >= lo & points <= hi
      out[[point_column]][i] <- list(ggchord_focus_points(points, lo, hi))
      if ("boundary_styles" %in% names(out)) {
        out$boundary_styles[i] <- list(out$boundary_styles[[i]][inside])
      }
    }
  }
  out$.lo <- out$.hi <- NULL
  out
}

ggchord_focus_feature_parts <- function(row, lo, hi, origin, caller) {
  if ("segments" %in% names(row) && !is.null(row$segments[[1L]]) &&
      !is.data.frame(row$segments[[1L]])) {
    ggchord_stop(caller, ": segments must be a list-column of data frames")
  }
  has_segments <- "segments" %in% names(row) &&
    is.data.frame(row$segments[[1L]]) && nrow(row$segments[[1L]]) > 0L
  if (has_segments) {
    segments <- ggchord_focus_nested_intervals(
      row$segments[[1L]], lo, hi, origin, caller)
    if (!nrow(segments)) return(NULL)
    if (!"segment_index" %in% names(segments)) segments$segment_index <- seq_len(nrow(segments))
    if (!"segment_type" %in% names(segments)) segments$segment_type <- "standard"
    if (!"color" %in% names(segments)) segments$color <- NA_character_
    if (!any(segments$segment_type != "gap")) return(NULL)
    row$segments <- I(list(segments))
    if ("segment_count" %in% names(row)) row$segment_count <- nrow(segments)
    # These fields derive from segments. Keeping a cached run would overwrite
    # the new top-level interval with a stale coordinate during rendering.
    row$.feature_runs <- row$.feature_boundaries <- row$.feature_boundary_styles <- NULL
  } else if (".feature_runs" %in% names(row) &&
      is.data.frame(row$.feature_runs[[1L]])) {
    runs <- ggchord_focus_nested_intervals(row$.feature_runs[[1L]], lo, hi,
      origin, caller, point_column = "boundaries")
    if (!nrow(runs)) return(NULL)
    row$.feature_runs <- I(list(runs))
  }
  if ("cleavage_arrows" %in% names(row)) {
    row$cleavage_arrows <- I(list(ggchord_focus_points(row$cleavage_arrows[[1L]], lo, hi)))
  }
  if (!has_segments && ".feature_boundaries" %in% names(row)) {
    points <- row$.feature_boundaries[[1L]]
    inside <- is.finite(points) & points >= lo & points <= hi
    row$.feature_boundaries <- I(list(ggchord_focus_points(points, lo, hi)))
    if (".feature_boundary_styles" %in% names(row)) {
      row$.feature_boundary_styles <- I(list(row$.feature_boundary_styles[[1L]][inside]))
    }
  }
  row
}

#' Focus ggchord data on one or more loci
#'
#' Crops sequence, gene and ribbon tables together and relocates retained
#' coordinates to start at one in each output locus. Multiple loci from the
#' same source sequence are supported; provide a unique `locus_id` column to
#' control their output sequence IDs.
#' Segmented annotations, cleavage positions and cached feature runs are
#' cropped and relocated together. Source intervals, source origins and
#' genomic left/right truncation flags are retained on the returned tables.
#' Alignment cropping is a linear display projection, not a new alignment:
#' identity is retained, and a scaled alignment length is marked as estimated.
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
#'   `loci`, and a compact `report` including retained/dropped source rows and
#'   truncated output rows. `.source_origin` is the original genomic position
#'   corresponding to local position one, including after repeated focusing.
#' @seealso \code{\link{ggchord_data_contract}}
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
  source <- ggchord_focus_source(seq_data)
  source_origin <- source$origin[seq_index] + focus$.focus_start - 1
  seq_out <- as.data.frame(seq_data[seq_index, , drop = FALSE],
                           stringsAsFactors = FALSE)
  seq_out$.source_accver <- source$accver[seq_index]
  seq_out$.source_start <- seq_out$.source_origin <- source_origin
  seq_out$.source_end <- source$origin[seq_index] + focus$.focus_end - 1
  seq_out$.source_length <- source$length[seq_index]
  seq_out$.truncated_left <- source_origin > 1
  seq_out$.truncated_right <- ifelse(is.finite(seq_out$.source_length),
    seq_out$.source_end < seq_out$.source_length, NA)
  seq_out$.focus_start <- focus$.focus_start
  seq_out$.focus_end <- focus$.focus_end
  seq_out$.locus_row <- focus$.locus_row
  seq_out$accver <- focus$.output_accver
  seq_out$length <- focus$.focus_end - focus$.focus_start + 1

  gene_out <- NULL
  gene_kept <- ribbon_kept <- integer()
  if (!is.null(gene_data)) {
    required <- c("accver", "start", "end", "strand")
    if (!is.data.frame(gene_data) || any(!required %in% names(gene_data)) ||
        !is.numeric(gene_data$start) || !is.numeric(gene_data$end) ||
        any(!is.finite(gene_data$start)) || any(!is.finite(gene_data$end))) {
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
        input_seq <- seq_index[j]
        row <- ggchord_focus_feature_parts(row, lo, hi, source$origin[input_seq], caller)
        if (is.null(row)) next
        row$accver <- focus$.output_accver[j]
        if (a <= b) {
          row$start <- clipped[1] - lo + 1
          row$end <- clipped[2] - lo + 1
        } else {
          row$start <- clipped[2] - lo + 1
          row$end <- clipped[1] - lo + 1
        }
        row$.source_start <- ggchord_focus_original(row, ".source_start",
          a + source$origin[input_seq] - 1)
        row$.source_end <- ggchord_focus_original(row, ".source_end",
          b + source$origin[input_seq] - 1)
        row$.source_accver <- source$accver[input_seq]
        row$.source_origin <- source_origin[j]
        row$.truncated_left <- isTRUE(row$.truncated_left) || clipped[1L] > glo
        row$.truncated_right <- isTRUE(row$.truncated_right) || clipped[2L] < ghi
        row$.source_row <- ggchord_focus_original(row, ".source_row", i)
        row$.locus_row <- focus$.locus_row[j]
        pieces[[length(pieces) + 1L]] <- row
        gene_kept <- c(gene_kept, i)
      }
    }
    gene_out <- if (length(pieces)) ggchord_rbind_fill(pieces) else {
      empty <- gene_data[0, , drop = FALSE]
      empty$.source_accver <- character(0)
      empty$.source_start <- empty$.source_end <- empty$.source_origin <- numeric()
      empty$.truncated_left <- empty$.truncated_right <- logical()
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
        query_source <- seq_index[qj]; subject_source <- seq_index[sj]
        for (end in c("q", "s")) {
          source_seq <- if (end == "q") query_source else subject_source
          locus_row <- if (end == "q") qj else sj
          original <- as.numeric(ribbon_data[i, paste0(end, c("start", "end"))])
          projected <- if (end == "q") qcoord else scoord
          for (coordinate in c("start", "end")) {
            name <- paste0(".source_", end, coordinate)
            row[[name]] <- ggchord_focus_original(row, name,
              original[if (coordinate == "start") 1L else 2L] +
                source$origin[source_seq] - 1)
          }
          row[[paste0(".source_", end, "accver")]] <- source$accver[source_seq]
          row[[paste0(".", end, "source_origin")]] <- source_origin[locus_row]
          row[[paste0(".", end, "truncated_left")]] <-
            isTRUE(row[[paste0(".", end, "truncated_left")]]) ||
              min(projected) > min(original)
          row[[paste0(".", end, "truncated_right")]] <-
            isTRUE(row[[paste0(".", end, "truncated_right")]]) ||
              max(projected) < max(original)
        }
        clipped <- tr[1L] > 0 || tr[2L] < 1
        row$.alignment_clipped <- isTRUE(row$.alignment_clipped) || clipped
        row$.alignment_projection <- if (row$.alignment_clipped) "linear" else "unchanged"
        row$.length_estimated <- isTRUE(row$.length_estimated) ||
          (clipped && "length" %in% names(ribbon_data))
        if ("length" %in% names(ribbon_data)) {
          row$.source_alignment_length <- ggchord_focus_original(row,
            ".source_alignment_length", ribbon_data$length[i])
        }
        if ("pident" %in% names(ribbon_data)) {
          row$.source_pident <- ggchord_focus_original(row, ".source_pident", ribbon_data$pident[i])
        }
        if ("length" %in% names(ribbon_data)) row$length <- max(1, round(ribbon_data$length[i] * diff(tr)))
        row$.source_row <- ggchord_focus_original(row, ".source_row", i)
        row$.query_locus_row <- focus$.locus_row[qj]
        row$.subject_locus_row <- focus$.locus_row[sj]
        pieces[[length(pieces) + 1L]] <- row
        ribbon_kept <- c(ribbon_kept, i)
      }
    }
    ribbon_out <- if (length(pieces)) ggchord_rbind_fill(pieces) else {
      empty <- ribbon_data[0, , drop = FALSE]
      empty$.source_row <- integer(0)
      for (end in c("q", "s")) {
        empty[[paste0(".source_", end, "accver")]] <- character()
        for (coordinate in c("start", "end")) {
          empty[[paste0(".source_", end, coordinate)]] <- numeric()
        }
        empty[[paste0(".", end, "source_origin")]] <- numeric()
        empty[[paste0(".", end, "truncated_left")]] <- logical()
        empty[[paste0(".", end, "truncated_right")]] <- logical()
      }
      empty$.alignment_clipped <- empty$.length_estimated <- logical()
      empty$.alignment_projection <- character()
      if ("length" %in% names(empty)) empty$.source_alignment_length <- numeric()
      if ("pident" %in% names(empty)) empty$.source_pident <- numeric()
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
      retained_input_rows = c(length(unique(seq_index)), length(unique(ribbon_kept)),
        length(unique(gene_kept))),
      dropped_input_rows = c(nrow(seq_data) - length(unique(seq_index)),
        if (is.null(ribbon_data)) 0L else nrow(ribbon_data) - length(unique(ribbon_kept)),
        if (is.null(gene_data)) 0L else nrow(gene_data) - length(unique(gene_kept))),
      truncated_output_rows = c(sum(seq_out$.truncated_left | seq_out$.truncated_right,
        na.rm = TRUE), if (is.null(ribbon_out)) 0L else sum(ribbon_out$.alignment_clipped),
        if (is.null(gene_out)) 0L else sum(gene_out$.truncated_left | gene_out$.truncated_right)),
      projection = c("source_coordinates", "linear_display", "source_coordinates"),
      boundary = boundary,
      stringsAsFactors = FALSE
    )
  )
}
