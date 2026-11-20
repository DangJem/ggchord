ggchord_link_orientation <- function(data, caller = "link") {
  n <- nrow(data)
  q <- as.character(data$qstrand %||% rep(".", n))
  s <- as.character(data$sstrand %||% rep(".", n))
  explicit <- as.character(data$direction %||% data$relative_strand %||% rep(NA_character_, n))
  known <- q %in% c("+", "-") & s %in% c("+", "-")
  orientation <- ifelse(known, ifelse(q == s, "same", "opposite"), "unknown")
  declared <- unname(c(same = "same", forward = "same", "+" = "same", "1" = "same",
    reverse = "opposite", opposite = "opposite", "-" = "opposite", "-1" = "opposite")[explicit])
  invalid <- !is.na(explicit) & !explicit %in% c("", ".", "?", "*", "0", "unknown") & is.na(declared)
  if (any(invalid)) ggchord_stop(caller, ": invalid relative direction")
  if (any(known & !is.na(declared) & orientation != declared))
    ggchord_stop(caller, ": relative direction conflicts with qstrand/sstrand")
  orientation[!is.na(declared)] <- declared[!is.na(declared)]
  if (!is.null(data$orientation_known)) {
    if (!is.logical(data$orientation_known)) ggchord_stop(caller, ": orientation_known must be logical")
    orientation[data$orientation_known %in% FALSE] <- "unknown"
  }
  orientation
}

ggchord_link_axis_pair <- function(value, caller, argument, positive_integer = FALSE) {
  if (!is.numeric(value) || !(length(value) %in% 1:2) || any(!is.finite(value)) ||
      any(value < 0) || (positive_integer && any(value < 1 | value != floor(value))))
    ggchord_stop(caller, ": invalid ", argument)
  if (length(value) == 1L && is.null(names(value))) return(setNames(rep(value, 2L), c("q", "s")))
  if (!setequal(names(value), c("q", "s"))) ggchord_stop(caller, ": ", argument, " must name q and s")
  value[c("q", "s")]
}

ggchord_link_group_columns <- function(data, group_by, caller) {
  if (is.null(group_by)) return(NULL)
  if (!is.character(group_by) || anyDuplicated(group_by) || any(!group_by %in% names(data)))
    ggchord_stop(caller, ": group_by must name distinct columns")
  data[group_by]
}

ggchord_link_pair_group <- function(data, orientation, group = NULL) {
  a <- list(data$qaccver, data$saccver, data$qstrand %||% rep(".", nrow(data)),
    data$sstrand %||% rep(".", nrow(data)), orientation)
  if (!is.null(group)) a <- c(a, if (is.list(group)) unname(as.list(group)) else list(group))
  if (any(lengths(a) != nrow(data))) ggchord_stop("link grouping: groups must align with rows")
  # Length prefixing keeps user-supplied IDs and missing group keys distinct.
  keys <- lapply(a, function(x) {
    x <- as.character(x)
    ifelse(is.na(x), "-1:", paste0(nchar(x, type = "bytes"), ":", x))
  })
  do.call(paste0, keys)
}

ggchord_link_bins <- function(data, sequences, bins, weight = NULL, group = NULL,
                              caller = "link bins") {
  bins <- ggchord_link_axis_pair(bins, caller, "bins", positive_integer = TRUE)
  weight <- weight %||% rep(1, nrow(data))
  if (!is.numeric(weight) || length(weight) != nrow(data) || any(!is.finite(weight) | weight < 0))
    ggchord_stop(caller, ": weight must provide a finite nonnegative number per record")
  sequences <- as_ggchord_seq(sequences)
  validate_ggchord_data(data, role = "link", sequences = sequences, action = "error",
    checks = c("schema", "coordinates", "units", "references"))
  ggchord_require_columns(data, ggchord_ribbon_required_columns(), caller)
  direction <- ggchord_link_orientation(data, caller)
  axis_bins <- function(side) {
    ids <- data[[paste0(side, "accver")]]
    source <- match(ids, sequences$accver)
    length <- sequences$length[source]
    midpoint <- (data[[paste0(side, "start")]] + data[[paste0(side, "end")]]) / 2
    bin_span <- length / bins[[side]]
    index <- pmin(bins[[side]], floor((midpoint - 1) / bin_span) + 1)
    list(index = index, span = bin_span, unit = sequences$sequence_unit[source])
  }
  q <- axis_bins("q")
  s <- axis_bins("s")
  pair <- ggchord_link_pair_group(data, direction, group)
  data.frame(qbin = q$index, sbin = s$index, area = q$span * s$span,
    q_unit = q$unit, s_unit = s$unit, weight = weight, direction = direction,
    pair_group = pair, bin_group = if (nrow(data)) paste(pair, q$index, s$index, sep = "/") else character())
}

ggchord_link_source_members <- function(data) {
  data$source_members %||% lapply(data$.source_row %||% seq_len(nrow(data)), function(x) x)
}

ggchord_summarize_link_groups <- function(data, groups, weight = NULL,
    extra_columns = "consistent", identity_summary = "none") {
  member <- ggchord_link_source_members(data)
  n <- lengths(groups)
  result <- data[vapply(groups, `[`, integer(1), 1L), , drop = FALSE]
  result$source_members <- I(lapply(groups, function(i) unique(unlist(member[i], recursive = FALSE))))
  result$n_members <- n
  result$is_summary <- n > 1L
  if (!is.null(weight)) result$sum_weight <- vapply(groups, function(i) sum(weight[i]), numeric(1))
  coordinates <- c("qstart", "qend", "sstart", "send")
  for (column in coordinates) {
    bound <- if (endsWith(column, "start")) min else max
    result[[column]] <- vapply(groups, function(i) bound(data[[column]][i]), numeric(1))
  }
  geometry <- c("qaccver", "saccver", coordinates, "qstrand", "sstrand", "direction", "orientation_known",
    "sequence_unit")
  original <- setdiff(names(data), c(geometry, "source_members", "n_members", "sum_weight", "is_summary"))
  provenance <- grepl("^(\\.source|source_)", original)
  per_record <- original %in% c("length", "link_id", "pident", "qsegments", "ssegments")
  for (column in original) {
    if (extra_columns == "list") {
      result[[column]] <- I(lapply(groups, function(i) data[[column]][i]))
    } else {
      for (j in which(result$is_summary)) {
        entries <- groups[[j]]
        different <- !all(vapply(entries, function(i)
          identical(data[[column]][i], data[[column]][entries[1L]]), logical(1)))
        if (different || column %in% original[provenance | per_record])
          result[[column]][j] <- ggchord_typed_na(data[[column]])
      }
    }
  }
  if (identity_summary == "weighted_mean") {
    result$summary_pident <- vapply(groups, function(i) {
      stats::weighted.mean(data$pident[i], data$length[i])
    }, numeric(1))
  }
  rownames(result) <- NULL
  ggchord_attach_data_report(result, "summarize_links",
    inputs = data_report(data), lineage = data.frame(output_row = seq_len(nrow(result)),
      source_members = result$source_members),
    metadata = list(source_records = data,
      grouped_rows = groups, identity_summary = identity_summary))
}

ggchord_bundle_binned_links <- function(data, grid, min_bundle) {
  if (!is.numeric(min_bundle) || length(min_bundle) != 1L ||
      !is.finite(min_bundle) || min_bundle < 2 || min_bundle != floor(min_bundle))
    ggchord_stop("link bundle: min_bundle must be an integer of at least two")
  groups <- split(seq_len(nrow(data)), grid$bin_group)
  groups <- unlist(lapply(groups, function(idx) {
    if (length(idx) < min_bundle || grid$direction[idx[1L]] == "unknown") as.list(idx) else list(idx)
  }), recursive = FALSE)
  groups <- groups[order(vapply(groups, min, integer(1)))]
  ggchord_summarize_link_groups(data, groups, grid$weight)
}

#' Aggregate interval links in a complete-reference grid
#' @param data Interval links.
#' @param sequences Complete reference dictionary.
#' @param bins Bin count for both axes or named counts for q and s.
#' @param min_bundle Minimum source records in a bundle.
#' @param weight A nonnegative numeric weight column, or NULL for unit weights.
#' @param group_by Optional extra grouping columns.
#' @return A link table with source members and a data report.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' bundled <- bundle_ggchord_links(links, seqs, bins=c(q=20, s=20), min_bundle=2)
#' ggchord(seqs) + geom_link_ribbon(data=bundled)
bundle_ggchord_links <- function(data, sequences, bins = 80, min_bundle = 2,
                                 weight = NULL, group_by = NULL) {
  caller <- "bundle_ggchord_links()"
  if (!is.null(weight) && (!is.character(weight) || length(weight) != 1L ||
      !weight %in% names(data))) ggchord_stop(caller, ": weight must name an input column")
  groups <- ggchord_link_group_columns(data, group_by, caller)
  links <- as_ggchord_link(data)
  row_weight <- if (!is.null(weight)) data[[weight]]
  grid <- ggchord_link_bins(links, sequences, bins, row_weight, groups, caller)
  ggchord_bundle_binned_links(links, grid, min_bundle)
}

#' Select original link rows
#' @param data A link table.
#' @param seq_ids,keep_pairs Source reference or pair selections.
#' @param min_pident,max_pident,min_length,max_evalue,min_bitscore Numeric filters.
#' @param min_query_coverage,min_subject_coverage Minimum source coverage.
#' @param drop_self_links Whether to remove same-reference links.
#' @param sort_by Column names, optionally prefixing descending columns with "-".
#' @return A link table with a selection report.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' strong <- filter_ggchord_links(links, min_pident=95, drop_self_links=TRUE)
#' ggchord(seqs) + geom_link_ribbon(data=strong)
#' data_report(strong)
filter_ggchord_links <- function(data, seq_ids = NULL, min_pident = NULL,
    max_pident = NULL, min_length = NULL, max_evalue = NULL, min_bitscore = NULL,
    min_query_coverage = NULL, min_subject_coverage = NULL, keep_pairs = NULL,
    drop_self_links = FALSE, sort_by = NULL) {
  result <- ggchord_filter_link_rows(data, seq_ids = seq_ids, min_pident = min_pident,
    max_pident = max_pident, min_length = min_length, max_evalue = max_evalue,
    min_bitscore = min_bitscore, min_query_coverage = min_query_coverage,
    min_subject_coverage = min_subject_coverage, keep_pairs = keep_pairs,
    drop_self_links = drop_self_links, sort_by = sort_by)
  rows <- attr(result$data, "source_rows")
  ggchord_attach_data_report(result$data, "filter_ggchord_links",
    inputs = data_report(data), lineage = data.frame(output_row = seq_along(rows),
      source_members = I(ggchord_link_source_members(data)[rows])),
    changes = result$report$removed_reasons, metadata = result$report)
}

#' Select one representative of duplicate links
#' @param data Interval links.
#' @param by Compare exact records, coordinates, or reciprocal overlap.
#' @param tolerance Endpoint tolerance in source units.
#' @param keep Representative selection: first, longest or best_pident.
#' @param min_reciprocal_overlap Minimum overlap fraction for overlap matching.
#' @return A link table with duplicate-to-representative lineage.
#' @export
#' @examples
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' duplicated_input <- rbind(links, links[1, ])
#' unique_links <- deduplicate_ggchord_links(duplicated_input, by="exact")
#' data_report(unique_links)$changes
deduplicate_ggchord_links <- function(data, by = "exact", tolerance = 0,
                                       keep = "first", min_reciprocal_overlap = 0.9) {
  by <- match.arg(by, c("exact", "coordinates", "overlap"))
  keep <- match.arg(keep, c("first", "longest", "best_pident"))
  links <- as_ggchord_link(data)
  if (by == "exact" && tolerance != 0) ggchord_stop("deduplicate_ggchord_links(): exact requires zero tolerance")
  # The overlap engine reads original direction evidence independently of
  # endpoint order. Run it within each explicit pair and strand group.
  orientation <- ggchord_link_orientation(links)
  grouped <- split(seq_len(nrow(links)), ggchord_link_pair_group(links, orientation))
  selected <- integer()
  duplicates <- list()
  for (rows in grouped) {
    part <- links[rows, , drop = FALSE]
    if (by == "exact") {
      compared <- part[setdiff(names(part), c(".source_row", ".source_line", ".source_file", "link_id",
        "source_qstart", "source_qend", "source_sstart", "source_send", "source_members"))]
      rank <- switch(keep, first = seq_along(rows), longest = order(-part$length), best_pident = order(-part$pident))
      if (!length(rank)) ggchord_stop("deduplicate_ggchord_links(): representative criterion requires its source column")
      kept <- rank[!duplicated(compared[rank, , drop = FALSE])]
      removed <- setdiff(seq_along(rows), kept)
      representative <- vapply(removed, function(i) {
        candidate <- compared[c(kept, i), , drop = FALSE]
        hit <- which(duplicated(candidate, fromLast = TRUE))[1L]
        kept[hit]
      }, integer(1))
      duplicates[[length(duplicates) + 1L]] <- data.frame(row = rows[removed],
        duplicate_of = rows[representative], reason = rep("exact", length(removed)))
    } else {
      reduced <- ggchord_deduplicate_link_rows(part, tolerance = tolerance, by = by, keep = keep,
        min_reciprocal_overlap = min_reciprocal_overlap)
      kept <- attr(reduced$data, "source_rows")
      removed <- reduced$report$removed
      if (nrow(removed)) {
        duplicates[[length(duplicates) + 1L]] <- data.frame(row = rows[removed$row],
          duplicate_of = rows[removed$duplicate_of], reason = rep(by, nrow(removed)))
      }
    }
    selected <- c(selected, rows[kept])
  }
  selected <- sort(selected)
  ggchord_attach_data_report(data[selected, , drop = FALSE], "deduplicate_ggchord_links",
    inputs = data_report(data),
    lineage = data.frame(output_row = seq_along(selected),
      source_members = I(ggchord_link_source_members(data)[selected])),
    changes = ggchord_rbind_fill(duplicates) %||%
      data.frame(row = integer(), duplicate_of = integer(), reason = character()))
}

#' Summarize adjacent compatible link intervals
#' @param data Interval links.
#' @param max_gap Maximum gap per endpoint, optionally named q and s.
#' @param max_pident_difference Optional maximum source identity difference.
#' @param group_by Additional grouping columns.
#' @param extra_columns Keep consistent values or a list of source values.
#' @param identity_summary Optionally report a length-weighted source mean.
#' @return A table with envelope spans, source members, and summary indicators.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' summaries <- merge_ggchord_links(links, max_gap=20, extra_columns="consistent")
#' ggchord(seqs) + geom_link_ribbon(data=summaries)
#' data_report(summaries)$lineage
merge_ggchord_links <- function(data, max_gap = 0, max_pident_difference = NULL,
    group_by = NULL, extra_columns = "consistent", identity_summary = "none") {
  caller <- "merge_ggchord_links()"
  max_gap <- ggchord_link_axis_pair(max_gap, caller, "max_gap")
  identity_summary <- match.arg(identity_summary, c("none", "weighted_mean"))
  extra_columns <- match.arg(extra_columns, c("consistent", "list"))
  grouping <- ggchord_link_group_columns(data, group_by, caller)
  links <- as_ggchord_link(data)
  if (!is.null(max_pident_difference) && (!is.numeric(max_pident_difference) ||
      length(max_pident_difference) != 1L || !is.finite(max_pident_difference) ||
      max_pident_difference < 0)) ggchord_stop(caller, ": invalid max_pident_difference")
  fields <- c(if (!is.null(max_pident_difference) || identity_summary == "weighted_mean") "pident",
    if (identity_summary == "weighted_mean") "length")
  ggchord_require_columns(links, fields, caller)
  ggchord_ribbon_numeric_input(links, fields, caller)
  if (identity_summary == "weighted_mean" && any(links$length <= 0))
    ggchord_stop(caller, ": source alignment lengths must be positive")
  direction <- ggchord_link_orientation(links, caller)
  groups <- split(seq_len(nrow(links)), ggchord_link_pair_group(links, direction, grouping))
  compatible <- lapply(groups, function(rows) {
    if (direction[rows[1L]] == "unknown") return(as.list(rows))
    reverse <- direction[rows[1L]] == "opposite"
    qlo <- links$qstart; qhi <- links$qend
    slo <- if (reverse) -links$send else links$sstart
    shi <- if (reverse) -links$sstart else links$send
    rows <- rows[order(qlo[rows], slo[rows], qhi[rows], shi[rows])]
    result <- list()
    active <- rows[1L]
    for (j in rows[-1L]) {
      last <- tail(active, 1L)
      progresses <- slo[j] >= slo[last] && shi[j] >= shi[last] && qhi[j] >= qhi[last]
      gaps_ok <- qlo[j] <= qhi[last] + max_gap[["q"]] + 1 &&
        slo[j] <= shi[last] + max_gap[["s"]] + 1
      ids_ok <- is.null(max_pident_difference) ||
        diff(range(links$pident[c(active, j)])) <= max_pident_difference
      if (progresses && gaps_ok && ids_ok) {
        active <- c(active, j)
      } else {
        result[[length(result) + 1L]] <- active
        active <- j
      }
    }
    c(result, list(active))
  })
  compatible <- unlist(compatible, recursive = FALSE)
  compatible <- compatible[order(vapply(compatible, min, integer(1)))]
  out <- ggchord_summarize_link_groups(links, compatible, extra_columns = extra_columns,
    identity_summary = identity_summary)
  out$qspan <- out$qend - out$qstart + 1
  out$sspan <- out$send - out$sstart + 1
  out
}
