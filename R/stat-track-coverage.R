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
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' intervals <- genes[c("accver", "start", "end")]
#'
#' ggchord(data = seqs) + geom_seq() +
#'   track_band("coverage", thickness = 5) +
#'   stat_coverage(data = intervals, binwidth = 100, track = "coverage") +
#'   scale_track_value("coverage", limits = c(0, 1),
#'                     labels = scales::label_percent()) +
#'   coord_collinear()
#'
#' ggchord(data = seqs) +
#'   geom_track_line(data = intervals, aes(y = after_stat(mean_depth)),
#'                   stat = "coverage", binwidth = 100) +
#'   coord_collinear()
stat_coverage <- function(mapping = NULL, data = NULL, geom = "track_col",
    position = "identity", binwidth = NULL, windows = NULL, bounds = "error",
    track = "default", ..., na.rm = FALSE, show.legend = NA,
    inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!is.character(geom) || length(geom) != 1L || !startsWith(geom, "track_"))
    ggchord_stop("stat_coverage(): choose a compatible track geometry")
  geom_track(mapping, data, geom = sub("^track_", "", geom), stat = "coverage",
    position = position, track = track, ..., binwidth = binwidth, windows = windows,
    bounds = bounds, na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
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
    member <- data$.ggchord_source_index[ord[i]]
    hits[[j]] <- sort(unique(member[!is.na(member)]))
    if (!length(i)) next
    a <- pmax(left[i], grid$window_start[j])
    b <- pmin(right[i], grid$window_end[j])
    part_groups <- split(seq_along(a), member)
    depth_sum[j] <- sum(vapply(part_groups, function(at) {
      aa <- a[at]; bb <- b[at]
      sum(pmax(0, bb - pmax(aa - 1, c(aa[1L] - 1, head(cummax(bb), -1L)))))
    }, numeric(1)))
    previous_end <- c(a[1L] - 1, head(cummax(b), -1L))
    covered[j] <- sum(pmax(0, b - pmax(a - 1, previous_end)))
  }
  grid$xmin <- grid$window_start - 1
  grid$xmax <- grid$window_end
  grid$x <- (grid$window_start + grid$window_end) / 2
  grid$start <- grid$window_start
  grid$end <- grid$window_end
  grid$width <- grid$denominator
  grid$covered <- covered
  grid$depth_sum <- depth_sum
  grid$coverage <- covered / grid$denominator
  grid$mean_depth <- depth_sum / grid$denominator
  grid$n_interval <- lengths(hits)
  grid$.ggchord_group_rows <- I(hits)
  grid
}

StatGgchordCoverage <- ggproto("StatGgchordCoverage", ggplot2::Stat,
  required_aes = c("xmin", "xmax", ".ggchord_source_index", ".ggchord_reference"),
  dropped_aes = c("xmin", "xmax", ".ggchord_source_index", ".ggchord_reference"),
  default_aes = aes(y = after_stat(coverage)),
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
  bad <- raw$.coverage_start < 1 | raw$.coverage_end > sequences$length[index]
  if (bounds == "error" && any(bad))
    ggchord_stop("Coverage interval outside sequence length in start/end; use bounds='clip' to intersect it with the sequence",
      ggchord_input_rows(bad))
  raw
}
