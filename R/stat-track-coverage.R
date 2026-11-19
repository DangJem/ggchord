#' Summarise interval coverage in sequence windows
#'
#' Each window uses inclusive, one-based endpoints in the reference sequence's
#' bp or aa units. Omitted interval/window units inherit that reference. A
#' discrete input style can form groups; use `aes(group = ...)` to choose them
#' explicitly or `aes(fill = after_stat(coverage))` to colour the pooled result.
#' Fixing a layer's fill also removes its inherited fill mapping before grouping.
#' Fixed-scale sample facets compute each observed series independently.
#' Give sample columns names distinct from statistical and geometry fields.
#' `after_stat()` expressions then see all statistical rows in the layer.
#' @inheritParams geom_track
#' @param geom Track geometry, e.g. `"track_col"`, `"track_line"`, or
#'   `"track_tile"` for a full-band interval coloured by a window statistic.
#' @param windows A data frame of `accver`, `start`, `end` and optional unique
#'   `window_id` and `sequence_unit`. Windows may overlap. Each independently
#'   records its contributing interval rows. Window coordinates refer to the
#'   input reference, including the local coordinates of a focused reference.
#' @param binwidth A positive whole number of sequence units, starting at 1.
#'   Supply exactly one of `binwidth` and `windows`.
#' @param ... Track geometry parameters, including fixed styles.
#' @return A track layer. Exported statistics include `window_id`, local
#'   `window_start` and `window_end`, original `source_accver`, `source_start`,
#'   `source_end`, `sequence_unit`, and the actual window length `denominator`.
#'   `covered` counts the interval union, `depth_sum` counts overlapping
#'   lengths with multiplicity, and `source_rows` indexes the layer source
#'   table. Dividing the two counts by `denominator` gives `coverage` and
#'   `mean_depth`. Windows without a hit in an observed group have zero counts.
#' @export
#' @examples
#' reference <- data.frame(accver = "A", length = 10L, sequence_unit = "bp")
#' intervals <- data.frame(accver = "A", start = c(3, 4, 10), end = c(6, 5, 10))
#' windows <- data.frame(accver = "A", start = c(1, 3), end = c(4, 6),
#'                       window_id = c("first", "overlap"))
#' ggchord(reference) + geom_seq() +
#'   track_band("coverage", width = grid::unit(10, "mm")) +
#'   stat_coverage(data = intervals, windows = windows, geom = "track_point",
#'                 track = "coverage", size = 3) +
#'   coord_collinear()
#' ggchord(reference) + geom_seq() +
#'   geom_track_line(aes(y = after_stat(mean_depth)),
#'                   data = intervals, stat = "coverage", binwidth = 4) +
#'   coord_collinear()
stat_coverage <- function(mapping = NULL, data = NULL, geom = "track_col",
    position = "identity", ..., binwidth = NULL, windows = NULL, track = NULL,
    id = NULL, show.legend = NA, inherit.aes = TRUE) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_track(mapping = mapping, data = data, geom = sub("^track_", "", geom),
    stat = "coverage", position = position, ..., binwidth = binwidth,
    windows = windows, track = track, id = id, show.legend = show.legend, inherit.aes = inherit.aes)
}

ggchord_coverage_grid <- function(sequences, binwidth, windows, observed) {
  unit <- ggchord_sequence_unit(sequences)
  if (is.null(windows)) {
    parts <- lapply(observed, function(i) {
      start <- as.double(seq(1, sequences$length[i], by = binwidth))
      data.frame(accver = as.character(sequences$accver[i]), start = start,
        end = pmin(start + binwidth - 1, sequences$length[i]))
    })
    windows <- ggchord_rbind_fill(parts)
    if (is.null(windows)) windows <- data.frame(accver = character(), start = numeric(), end = numeric())
  }
  ggchord_require_columns(windows, c("accver", "start", "end"), "Coverage windows")
  ggchord_check_interval_unit(windows, unit, "Coverage windows")
  index <- match(windows$accver, sequences$accver)
  if (anyNA(index)) ggchord_stop("Coverage windows contain an unknown reference accession in accver",
    ggchord_input_rows(is.na(index)))
  for (key in c("start", "end")) {
    value <- windows[[key]]
    bad <- if (is.numeric(value)) !is.finite(value) | value != floor(value) else rep(TRUE, nrow(windows))
    if (!is.numeric(value) || any(bad))
      ggchord_stop("Coverage windows require finite whole ", key, " coordinates", ggchord_input_rows(bad))
  }
  bad <- windows$start < 1 | windows$start > windows$end | windows$end > sequences$length[index]
  if (any(bad)) ggchord_stop("Coverage windows require 1 <= start <= end <= reference length", ggchord_input_rows(bad))
  start <- as.double(windows$start)
  end <- as.double(windows$end)
  source <- ggchord_focus_source(sequences)
  id <- if (!is.null(windows$window_id)) as.character(windows$window_id) else
    if (length(start)) paste(windows$accver, start, end, sep = ":") else character()
  bad <- is.na(id) | !nzchar(id) | duplicated(id) | duplicated(id, fromLast = TRUE)
  if (any(bad)) ggchord_stop("Coverage windows need unique non-missing window_id values", ggchord_input_rows(bad))
  grid <- data.frame(window_id = id, window_start = start, window_end = end,
    source_accver = source$accver[index], source_start = source$origin[index] + start - 1,
    source_end = source$origin[index] + end - 1, sequence_unit = rep(unit, length(start)),
    denominator = end - start + 1)
  split(grid, factor(index, levels = seq_len(nrow(sequences))))
}

# Each intersecting pair contains the later of its two starts. The first
# pass includes ties and the second excludes them, for one entry per pair.
ggchord_coverage_windows <- function(data, grid) {
  ord <- order(data$xmin)
  left <- data$xmin[ord]
  right <- data$xmax[ord]
  win <- order(grid$window_start)
  window_starts <- grid$window_start[win]
  before <- findInterval(grid$window_start - 1, left)
  n_intervals <- findInterval(grid$window_end, left) - before
  through <- findInterval(left, window_starts)
  n_windows <- findInterval(right, window_starts) - through
  consecutive <- function(offset, count) rep(offset, count) + sequence(count)
  interval <- c(consecutive(before, n_intervals), rep(seq_along(left), n_windows))
  window <- c(rep(seq_len(nrow(grid)), n_intervals), win[consecutive(through, n_windows)])
  intersections <- split(interval, factor(window, levels = seq_len(nrow(grid))))
  hits <- vector("list", nrow(grid))
  covered <- depth_sum <- numeric(nrow(grid))
  for (j in seq_along(intersections)) {
    i <- sort(intersections[[j]])
    hits[[j]] <- sort(data$.ggchord_source_index[ord[i]])
    if (!length(i)) next
    a <- pmax(left[i], grid$window_start[j])
    b <- pmin(right[i], grid$window_end[j])
    depth_sum[j] <- sum(b - a + 1)
    previous_end <- c(a[1L] - 1, utils::head(cummax(b), -1L))
    covered[j] <- sum(pmax(0, b - pmax(a - 1, previous_end)))
  }
  grid$xmin <- grid$window_start - 1
  grid$xmax <- grid$window_end
  grid$x <- (grid$xmin + grid$xmax) / 2
  grid$width <- grid$denominator
  grid$covered <- covered
  grid$depth_sum <- depth_sum
  grid$coverage <- covered / grid$denominator
  grid$mean_depth <- depth_sum / grid$denominator
  grid$n_interval <- lengths(hits)
  grid$.ggchord_group_rows <- I(hits)
  grid
}

StatGgchordCoverage <- ggplot2::ggproto("StatGgchordCoverage", ggplot2::Stat,
  required_aes = c("xmin", "xmax", ".ggchord_source_index", ".ggchord_reference"),
  dropped_aes = c("xmin", "xmax", ".ggchord_source_index", ".ggchord_reference"),
  default_aes = ggplot2::aes(y = after_stat(coverage)),
  compute_group = function(data, scales, window_grid) {
    ggchord_coverage_windows(data, window_grid[[data$.ggchord_reference[1L]]])
  })

ggchord_coverage_input <- function(input, mapping, sequences, bounds, index) {
  raw <- input
  bad <- !is.finite(sequences$length) | sequences$length < 1 | sequences$length != floor(sequences$length)
  if (any(bad)) ggchord_stop("Coverage requires positive whole sequence lengths in the reference table",
    ggchord_input_rows(bad))
  for (name in c("start", "end")) {
    q <- mapping[[name]]
    if (is.null(q)) ggchord_require_columns(raw, name, "geom_track(stat='coverage')")
    if (ggchord_has_stage(q)) ggchord_stop("Coverage endpoints must be available before statistics")
    value <- if (is.null(q)) input[[name]] else ggchord_eval_input_mapping(name, q, input, "Coverage intervals")
    if (is.character(value) || is.factor(value)) value <- suppressWarnings(as.numeric(as.character(value)))
    if (!length(value) %in% c(1L, nrow(raw)))
      ggchord_stop("Coverage ", name, " must supply one coordinate per interval (got ",
        length(value), " for ", nrow(raw), " rows)")
    value <- rep_len(value, nrow(raw))
    bad <- if (is.numeric(value)) !is.finite(value) | value != floor(value) else rep(TRUE, nrow(raw))
    if (!is.numeric(value) || any(bad)) ggchord_stop("Coverage ", name,
      " must supply finite whole coordinates for each interval", ggchord_input_rows(bad))
    raw[[paste0(".coverage_", name)]] <- as.double(value)
  }
  bad <- raw$.coverage_start > raw$.coverage_end
  if (any(bad)) ggchord_stop("Coverage requires start <= end; split origin-crossing intervals before aggregation",
    ggchord_input_rows(bad))
  parts <- raw[["segments"]]
  if (!is.null(parts)) {
    if (!is.list(parts) || is.data.frame(parts) || length(parts) != nrow(raw))
      ggchord_stop("Coverage expects segments to be a list-column of intervals")
    # Readers such as read_genbank() keep a segment table even for a single
    # ordinary interval. A compound envelope, in contrast, can include gaps.
    simple <- vapply(seq_along(parts), function(i) {
      s <- parts[[i]]
      if (is.null(s) || (is.data.frame(s) && !nrow(s))) return(TRUE)
      is.data.frame(s) && nrow(s) == 1L && all(c("start", "end") %in% names(s)) &&
        isTRUE(s$start == raw$.coverage_start[i]) && isTRUE(s$end == raw$.coverage_end[i]) &&
        (is.null(s$segment_type) || isTRUE(s$segment_type != "gap"))
    }, logical(1))
    if (any(!simple)) ggchord_stop("Coverage needs one interval per row; expand the segments column before aggregation",
      ggchord_input_rows(!simple))
  }
  bad <- raw$.coverage_start < 1 | raw$.coverage_end > sequences$length[index]
  if (bounds == "error" && any(bad))
    ggchord_stop("Coverage interval outside sequence length in start/end; use bounds='clip' to intersect it with the sequence",
      ggchord_input_rows(bad))
  raw
}
