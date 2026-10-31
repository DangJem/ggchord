# Independent focused windows are the units of a piecewise display. They never
# become one artificial contiguous genomic sequence.

#' Prepare independent windows for a segmented linear display
#'
#' Crops through [focus_ggchord_data()] while retaining one sequence per
#' window. Pass `segments` to [coord_collinear()] to place those sequences in
#' shared display rows with explicit gaps and local reversals. Coordinates in
#' the returned tables remain window-local and are never concatenated.
#'
#' `cross_gap = "split"` retains only the parts of each ribbon's linear
#' endpoint correspondence that lie in both selected windows. Omitted regions
#' are never bridged by one ribbon polygon. `"drop"` retains only links fully
#' contained in a query and subject window. This is a display projection, not
#' realignment. Point links survive only when both endpoints are retained.
#' Overlapping windows can deliberately duplicate a source row; provenance
#' identifies each copy, and `max_links` bounds the candidate expansion before
#' allocation. Order in the window table is display order within each row.
#'
#' @param seq_data Full sequence table with unique `accver` and integer lengths.
#' @param ribbon_data Optional paired interval table.
#' @param gene_data Optional gene or generic feature table, including segments.
#' @param windows Table with `accver`, `start`, `end`, unique `locus_id`,
#'   optional `row_id` (defaults to source `accver`) and `orientation` (1 or -1,
#'   defaults to 1). Bounds must be ordered, integral and within the sequence.
#' @param point_data Optional point links with `qaccver`, `saccver`, `qpos`, `spos`.
#' @param cross_gap How to handle interval links clipped by selected windows.
#' @param max_links Positive integer upper bound on candidate ribbon/point
#'   copies. Increase explicitly for larger, intentionally overlapping views.
#' @return A list containing focused `seq_data`, `ribbon_data`, `gene_data`,
#'   `point_data`, `segments` for the coordinate system, reversible source to
#'   window `transforms`, and a `report` with input row numbers for every
#'   output connection and omitted input rows. Each window retains its own
#'   source accession, origin, full length and truncation flags.
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' s <- data.frame(accver = "A", length = 1000)
#' w <- data.frame(accver = "A", start = c(100, 700), end = c(300, 900),
#'   locus_id = c("left", "right"), row_id = "isolate", orientation = c(1, -1))
#' view <- segment_ggchord_data(s, windows = w)
#' ggchord(view$seq_data) + geom_seq() +
#'   coord_collinear(segments = view$segments, segment_gap = 50)
segment_ggchord_data <- function(seq_data, ribbon_data = NULL, gene_data = NULL,
    windows, point_data = NULL, cross_gap = c("split", "drop"), max_links = 10000) {
  caller <- "segment_ggchord_data()"
  cross_gap <- match.arg(cross_gap)
  seq_data <- as_ggchord_seq(seq_data)
  if (!is.data.frame(windows) || !nrow(windows)) ggchord_stop(caller, ": windows must be a nonempty table")
  ggchord_require_columns(windows, c("accver", "start", "end", "locus_id"), caller)
  windows <- as.data.frame(windows)
  windows$accver <- ggchord_conversion_ids(windows$accver, caller)
  windows$locus_id <- ggchord_conversion_ids(windows$locus_id, caller)
  if (anyDuplicated(windows$locus_id)) ggchord_stop(caller, ": locus_id must be unique")
  ggchord_conversion_coordinates(windows, c("start", "end"), caller)
  index <- match(windows$accver, seq_data$accver)
  if (anyNA(index) || any(windows$start > windows$end) ||
      any(windows$end > seq_data$length[index])) ggchord_stop(caller,
    ": windows must name known sequences and ordered bounds within sequence length")
  if (!"row_id" %in% names(windows)) windows$row_id <- windows$accver
  windows$row_id <- ggchord_conversion_ids(windows$row_id, caller)
  if (!"orientation" %in% names(windows)) windows$orientation <- 1
  if (!is.numeric(windows$orientation) || anyNA(windows$orientation) ||
      any(!windows$orientation %in% c(-1, 1))) ggchord_stop(caller, ": window orientation must be 1 or -1")
  if (!is.numeric(max_links) || length(max_links) != 1L || !is.finite(max_links) ||
      max_links < 1 || max_links != floor(max_links)) ggchord_stop(caller, ": max_links must be a positive integer")
  if (!is.null(gene_data)) {
    gene_data <- as_ggchord_feature(gene_data)
    if (any(gene_data$start > gene_data$end)) ggchord_stop(caller,
      ": origin-spanning features need explicit ordered segments and a genomic envelope before window selection")
    at <- match(gene_data$accver, seq_data$accver)
    if (anyNA(at) || any(gene_data$end > seq_data$length[at])) ggchord_stop(caller,
      ": features must lie within known source sequences")
  }
  if (!is.null(ribbon_data)) ribbon_data <- as_ggchord_link(ribbon_data)
  if (!is.null(point_data)) {
    if (!is.data.frame(point_data)) ggchord_stop(caller, ": point_data must be a data.frame")
    ggchord_require_columns(point_data, c("qaccver", "saccver", "qpos", "spos"), caller)
    ggchord_conversion_coordinates(point_data, c("qpos", "spos"), caller)
    for (side in c("q", "s")) {
      ids <- ggchord_conversion_ids(point_data[[paste0(side, "accver")]], caller)
      point_data[[paste0(side, "accver")]] <- ids
      at <- match(ids, seq_data$accver)
      if (anyNA(at) || any(point_data[[paste0(side, "pos")]] > seq_data$length[at]))
        ggchord_stop(caller, ": point endpoints must be within known sequences")
    }
  }
  hits <- function(id, a, b, whole = FALSE) {
    which(windows$accver == id & if (whole) {
      windows$start <= min(a, b) & windows$end >= max(a, b)
    } else windows$start <= max(a, b) & windows$end >= min(a, b))
  }
  candidate_count <- 0
  retained <- if (is.null(ribbon_data)) integer() else seq_len(nrow(ribbon_data))
  if (!is.null(ribbon_data)) {
    ggchord_check_ribbon_tables(seq_data, ribbon_data, caller)
    for (side in c("q", "s")) {
      at <- match(ribbon_data[[paste0(side, "accver")]], seq_data$accver)
      if (any(pmax(ribbon_data[[paste0(side, "start")]], ribbon_data[[paste0(side, "end")]]) >
          seq_data$length[at])) ggchord_stop(caller, ": ribbon endpoints must lie within known source sequences")
    }
    whole <- cross_gap == "drop"
    for (i in seq_len(nrow(ribbon_data))) {
      q <- hits(ribbon_data$qaccver[i], ribbon_data$qstart[i], ribbon_data$qend[i], whole)
      s <- hits(ribbon_data$saccver[i], ribbon_data$sstart[i], ribbon_data$send[i], whole)
      count <- length(q) * as.double(length(s))
      if (!count) {
        retained <- retained[retained != i]
        next
      }
      if (whole) count <- length(hits(ribbon_data$qaccver[i], ribbon_data$qstart[i], ribbon_data$qend[i])) *
        as.double(length(hits(ribbon_data$saccver[i], ribbon_data$sstart[i], ribbon_data$send[i])))
      candidate_count <- candidate_count + count
      if (candidate_count > max_links) ggchord_stop(caller,
        ": candidate ribbon expansion exceeds max_links before allocation; use fewer windows or raise the limit explicitly")
    }
  }
  point_hits <- list()
  if (!is.null(point_data)) for (i in seq_len(nrow(point_data))) {
    q <- hits(point_data$qaccver[i], point_data$qpos[i], point_data$qpos[i])
    s <- hits(point_data$saccver[i], point_data$spos[i], point_data$spos[i])
    candidate_count <- candidate_count + length(q) * as.double(length(s))
    if (candidate_count > max_links) ggchord_stop(caller, ": candidate connection expansion exceeds max_links before allocation")
    point_hits[[i]] <- list(q = q, s = s)
  }
  # Add a private input identity independently of any retained source lineage.
  tag <- ".segment_input_row"
  if (!is.null(ribbon_data) && tag %in% names(ribbon_data)) ggchord_stop(caller, ": reserved input column: ", tag)
  if (!is.null(ribbon_data)) ribbon_data[[tag]] <- seq_len(nrow(ribbon_data))
  focused <- focus_ggchord_data(seq_data,
    ribbon_data = if (is.null(ribbon_data)) NULL else ribbon_data[retained, , drop = FALSE],
    gene_data = gene_data, loci = windows)
  ribbon_rows <- if (is.null(focused$ribbon_data)) integer() else focused$ribbon_data[[tag]]
  if (!is.null(focused$ribbon_data)) focused$ribbon_data[[tag]] <- NULL
  # Whole-link policy must be evaluated for each window pair, not merely for
  # an input row which happens to fit some other overlapping selected window.
  if (cross_gap == "drop" && length(ribbon_rows)) {
    q <- focused$ribbon_data$.query_locus_row; s <- focused$ribbon_data$.subject_locus_row
    original <- ribbon_data[ribbon_rows, , drop = FALSE]
    keep <- pmin(original$qstart, original$qend) >= windows$start[q] &
      pmax(original$qstart, original$qend) <= windows$end[q] &
      pmin(original$sstart, original$send) >= windows$start[s] &
      pmax(original$sstart, original$send) <= windows$end[s]
    focused$ribbon_data <- focused$ribbon_data[keep, , drop = FALSE]
    ribbon_rows <- ribbon_rows[keep]
  }
  points <- list(); point_rows <- integer()
  source <- ggchord_focus_source(seq_data)
  if (!is.null(point_data)) for (i in seq_len(nrow(point_data))) {
    for (q in point_hits[[i]]$q) for (s in point_hits[[i]]$s) {
      row <- point_data[i, , drop = FALSE]
      for (side in c("q", "s")) {
        at <- match(row[[paste0(side, "accver")]], seq_data$accver)
        id <- paste0(".source_", side, "accver")
        pos <- paste0(".source_", side, "pos")
        row[[id]] <- ggchord_focus_original(row, id, source$accver[at])
        row[[pos]] <- ggchord_focus_original(row, pos,
          row[[paste0(side, "pos")]] + source$origin[at] - 1)
      }
      row$qaccver <- windows$locus_id[q]; row$saccver <- windows$locus_id[s]
      row$qpos <- row$qpos - windows$start[q] + 1
      row$spos <- row$spos - windows$start[s] + 1
      points[[length(points) + 1L]] <- row; point_rows <- c(point_rows, i)
    }
  }
  segments <- data.frame(accver = focused$seq_data$accver, row_id = windows$row_id,
    orientation = windows$orientation, .truncated_left = focused$seq_data$.truncated_left,
    .truncated_right = focused$seq_data$.truncated_right, stringsAsFactors = FALSE)
  transforms <- data.frame(accver = focused$seq_data$accver, row_id = windows$row_id,
    source_accver = focused$seq_data$.source_accver,
    source_start = focused$seq_data$.source_origin,
    source_end = focused$seq_data$.source_end,
    local_start = 1, local_end = focused$seq_data$length,
    orientation = windows$orientation, stringsAsFactors = FALSE)
  list(seq_data = focused$seq_data, ribbon_data = focused$ribbon_data,
    gene_data = focused$gene_data,
    point_data = if (length(points)) ggchord_rbind_fill(points) else
      if (is.null(point_data)) NULL else point_data[0, , drop = FALSE],
    segments = segments, transforms = transforms, report = list(
      cross_gap = cross_gap, focus = focused$report,
      ribbon_input_rows = ribbon_rows, point_input_rows = point_rows,
      omitted_ribbon_rows = setdiff(if (is.null(ribbon_data)) integer() else seq_len(nrow(ribbon_data)), ribbon_rows),
      omitted_point_rows = setdiff(if (is.null(point_data)) integer() else seq_len(nrow(point_data)), point_rows)))
}
