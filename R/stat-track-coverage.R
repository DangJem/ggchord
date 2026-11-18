# Interval statistics in the sequence's one-based, inclusive source units.
# Work proportional to windows and interval/window intersections, never bases.
ggchord_coverage_windows <- function(start, end, rows, sequence_length, binwidth) {
  window_start <- as.double(seq(1, sequence_length, by = binwidth))
  window_end <- pmin(window_start + binwidth - 1, sequence_length)
  n <- length(window_start)
  start <- pmax(1, start)
  end <- pmin(sequence_length, end)
  keep <- which(rows > 0 & start <= end)
  left <- start[keep]
  right <- end[keep]
  first <- floor((left - 1) / binwidth) + 1
  last <- floor((right - 1) / binwidth) + 1
  repeats <- last - first + 1
  idx <- rep(seq_along(keep), repeats)
  window <- rep(first, repeats) + sequence(repeats) - 1
  hits <- split(idx, factor(window, levels = seq_len(n)))
  covered <- depth_sum <- numeric(n)
  for (j in which(lengths(hits) > 0)) {
    i <- hits[[j]]
    a <- pmax(left[i], window_start[j])
    b <- pmin(right[i], window_end[j])
    depth_sum[j] <- sum(b - a + 1)
    by_start <- order(a)
    a <- a[by_start]
    b <- b[by_start]
    # Each interval adds only its extension beyond the preceding union.
    previous_end <- c(a[1L] - 1, utils::head(cummax(b), -1L))
    covered[j] <- sum(pmax(0, b - pmax(a - 1, previous_end)))
  }
  width <- window_end - window_start + 1
  data.frame(x = (window_start - 1 + window_end) / 2, width = width,
    xmin = window_start - 1, xmax = window_end, window_start, window_end,
    covered, coverage = covered / width, mean_depth = depth_sum / width,
    n_interval = lengths(hits, use.names = FALSE),
    .ggchord_group_rows = I(unname(lapply(hits, function(i) rows[keep[i]]))))
}

StatGgchordCoverage <- ggplot2::ggproto("StatGgchordCoverage", ggplot2::Stat,
  required_aes = c("xmin", "xmax", ".ggchord_source_index"),
  dropped_aes = c("xmin", "xmax", ".ggchord_source_index"),
  default_aes = ggplot2::aes(y = after_stat(coverage)),
  compute_group = function(data, scales, sequence_length, binwidth) {
    ggchord_coverage_windows(data$xmin, data$xmax, data$.ggchord_source_index,
      sequence_length, binwidth)
  })

ggchord_coverage_input <- function(raw, mapping, sequences, bounds) {
  if (is.list(raw[["segments"]]) && any(vapply(raw[["segments"]], function(x)
      is.data.frame(x) && nrow(x) > 0L, logical(1))))
    ggchord_stop("Coverage needs one interval per row; expand the segments column before aggregation")
  ids <- as.character(raw$accver)
  index <- match(ids, sequences$accver)
  if (anyNA(index)) ggchord_stop("Coverage needs a known accession for each interval")
  if (any(!is.finite(sequences$length) | sequences$length < 1 |
      sequences$length != floor(sequences$length)))
    ggchord_stop("Coverage requires positive whole sequence lengths")
  for (name in c("start", "end")) {
    q <- mapping[[name]]
    if (is.null(q)) ggchord_require_columns(raw, name, "geom_track(stat='coverage')")
    if (ggchord_has_stage(q)) ggchord_stop("Coverage endpoints must be available before statistics")
    value <- if (is.null(q)) raw[[name]] else rlang::eval_tidy(q, raw)
    if (is.character(value) || is.factor(value)) value <- suppressWarnings(as.numeric(as.character(value)))
    if (!is.numeric(value) || !length(value) %in% c(1L, nrow(raw)) ||
        any(!is.finite(value) | value != floor(value)))
      ggchord_stop("Coverage ", name, " must supply finite whole coordinates for each interval")
    raw[[paste0(".coverage_", name)]] <- rep_len(as.double(value), nrow(raw))
  }
  if (any(raw$.coverage_start > raw$.coverage_end))
    ggchord_stop("Coverage requires start <= end; split origin-crossing intervals before aggregation")
  if (bounds == "error" && any(raw$.coverage_start < 1 | raw$.coverage_end > sequences$length[index]))
    ggchord_stop("Coverage interval outside sequence length; use bounds='clip' to intersect it with the sequence")
  raw
}
