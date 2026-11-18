# From the repository root:
# Rscript --vanilla tools/long-term-audit/v0.50.0-candidates.R [results.rds]
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
  library(ggplot2)
  library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
dot <- flagship_dotplot_example
sequences <- dot$seq_data
links <- dot$links
ids <- sequences$accver
bin_width <- 5000L

# Dotplot source positions are one-based point centres. Compare indexing by
# whole coordinate origin with indexing by the half-integer cell edge.
source_bin <- function(position, origin = 1) floor((position - origin) / bin_width)
midpoints <- links |>
  transmute(qaccver, saccver, source_hsp_id,
    qmid = (qstart + qend) / 2, smid = (sstart + send) / 2,
    q_bin = source_bin(qmid), s_bin = source_bin(smid),
    q_edge_bin = source_bin(qmid, 0.5), s_edge_bin = source_bin(smid, 0.5),
    q_first_bin = source_bin(qstart), s_first_bin = source_bin(sstart))

make_cells <- function(q_bin, s_bin) {
  midpoints |>
    transmute(qaccver, saccver, source_hsp_id, q_bin = {{ q_bin }}, s_bin = {{ s_bin }}) |>
    group_by(qaccver, saccver, q_bin, s_bin) |>
    summarise(count = n(), source_hsp_id = list(sort(source_hsp_id)), .groups = "drop") |>
    mutate(q_start = 1 + q_bin * bin_width, s_start = 1 + s_bin * bin_width,
      q_end = pmin(q_start + bin_width - 1, sequences$length[match(qaccver, ids)]),
      s_end = pmin(s_start + bin_width - 1, sequences$length[match(saccver, ids)]))
}
cells <- make_cells(q_bin, s_bin)
edge_cells <- make_cells(q_edge_bin, s_edge_bin)
first_cells <- make_cells(q_first_bin, s_first_bin)
stopifnot(nrow(links) == 32L, nrow(cells) == 23L, identical(cells, edge_cells),
  identical(sort(unlist(cells$source_hsp_id)), sort(links$source_hsp_id)))
peak <- cells |>
  filter(qaccver == "MT108731.1", saccver == "MT118296.1", q_bin == 4, s_bin == 0)
stopifnot(nrow(peak) == 1L, peak$count == 4L)

base_dot <- ggchord(sequences) + geom_alignment_segment(data = links)
auto <- export_ggchord_layout(base_dot + coord_dotplot(gap = dot$axis_gap))
stopifnot(all(table(auto$axis_layout$axis) == 3L), nrow(auto$alignment) == 32L)
source_fields <- c("qaccver", "saccver", "qstart", "qend", "sstart", "send", "source_hsp_id")
orientations <- expand.grid(x = c(1, -1), y = c(1, -1))
projected <- lapply(seq_len(nrow(orientations)), function(i) {
  directions <- orientations[i, ]
  p <- base_dot + coord_dotplot(x_seqs = ids, y_seqs = ids,
    x_orientation = directions$x, y_orientation = directions$y, gap = dot$axis_gap)
  export <- export_ggchord_layout(p)
  observed <- export$alignment
  stopifnot(all(table(export$axis_layout$axis) == 4L),
    identical(observed[source_fields], links[source_fields]))
  for (axis in c("x", "y")) {
    query <- axis == "x"
    accver <- if (query) observed$qaccver else observed$saccver
    centre <- if (query) midpoints$qmid else midpoints$smid
    midpoint <- (observed[[axis]] + observed[[paste0(axis, "end")]]) / 2
    forward <- map_ggchord_position(p, accver, centre, axis = axis)[[axis]]
    inverse <- map_ggchord_position(p, accver, midpoint, inverse = TRUE, axis = axis)
    frame <- filter(export$axis_layout, .data$axis == .env$axis)
    lengths <- frame$length[match(accver, frame$accver)]
    start <- 1 + source_bin(centre, 0.5) * bin_width
    end <- pmin(start + bin_width - 1, lengths)
    direction <- directions[[axis]]
    left <- map_ggchord_position(p, accver, start, axis = axis)[[axis]] - direction / 2
    right <- map_ggchord_position(p, accver, end, axis = axis)[[axis]] + direction / 2
    stopifnot(all(midpoint >= pmin(left, right)), all(midpoint <= pmax(left, right)))
    stopifnot(all(midpoint == forward), all(inverse$position == centre))
    edge <- map_ggchord_position(p, ids[1L], c(5000, 5000.5, 5001), axis = axis)[[axis]]
    stopifnot(all(diff(edge) == direction * 0.5))
  }
  export$axis_layout
})
point_check <- data.frame(position = c(5000, 5000.5, 5001)) |>
  mutate(first_base_origin = 1 + source_bin(position), cell_edge_origin = 1 + source_bin(position, 0.5))
stopifnot(identical(point_check$first_base_origin, c(1, 1, 2)),
  identical(point_check$cell_edge_origin, c(1, 2, 2)))
last_cells <- sequences |>
  transmute(accver, start = 1 + source_bin(length) * bin_width, end = length,
    denominator = end - start + 1)

# One reference per query and one observed comparison per subject. Compute
# the independent answer by splitting each window at interval boundaries.
reference <- filter(sequences, accver == "MT108731.1")
intervals <- links |>
  filter(qaccver == reference$accver) |>
  transmute(accver = qaccver, comparison = saccver, start = pmin(qstart, qend),
    end = pmax(qstart, qend), source_hsp_id)
groups <- unique(intervals$comparison)
windows <- data.frame(window_start = seq(1, reference$length, by = bin_width)) |>
  mutate(window_end = pmin(window_start + bin_width - 1, reference$length))
coverage_rows <- bind_rows(lapply(groups, function(group) {
  rows <- filter(intervals, comparison == group)
  bind_rows(lapply(seq_len(nrow(windows)), function(i) {
    win <- windows[i, , drop = FALSE]
    hits <- filter(rows, start <= win$window_end, end >= win$window_start)
    left <- pmax(hits$start, win$window_start)
    right <- pmin(hits$end + 1, win$window_end + 1)
    edges <- sort(unique(c(win$window_start, left, right, win$window_end + 1)))
    depth <- vapply(head(edges, -1L), function(x) sum(left <= x & right > x), integer(1))
    data.frame(comparison = group, win, denominator = win$window_end - win$window_start + 1,
      covered = sum(diff(edges)[depth > 0L]), depth_sum = sum(diff(edges) * depth),
      source_hsp_id = I(list(sort(hits$source_hsp_id))))
  }))
}))
stopifnot(nrow(coverage_rows) == 39L, all(coverage_rows$denominator[c(13, 26, 39)] == 4323),
  identical(as.integer(table(intervals$comparison)[groups]), c(7L, 15L, 1L)))

base_tracks <- ggchord(reference) + geom_seq() + coord_collinear()
separate_tracks <- function(order, width) {
  Reduce(function(p, i) {
    group <- order[i]
    p + track_band(group, order = i, width = grid::unit(width, "mm")) +
      stat_coverage(data = filter(intervals, comparison == group), binwidth = bin_width,
        track = group, id = group) +
      scale_track_value_continuous(group, limits = c(0, 1))
  }, seq_along(order), init = base_tracks)
}
grouped <- base_tracks + track_band("coverage") +
  stat_coverage(aes(group = comparison, colour = comparison), data = intervals,
    binwidth = bin_width, geom = "track_line", track = "coverage", id = "grouped") +
  scale_track_value_continuous("coverage", limits = c(0, 1))
faceted <- base_tracks + facet_wrap(vars(comparison)) + track_band("coverage") +
  stat_coverage(data = intervals, binwidth = bin_width, id = "faceted", track = "coverage")
fixed <- faceted + scale_track_value_continuous("coverage", limits = c(0, 1))

window_fields <- c("comparison", "window_start", "window_end", "denominator", "covered", "depth_sum")
check_coverage <- function(p, comparison) {
  out <- export_ggchord_layout(p, data_stage = c("source", "statistic"))
  rows <- bind_rows(lapply(names(out$track_statistics), function(id) {
    data <- out$track_statistics[[id]]
    stage <- out$layer_data[[id]]
    stopifnot(nrow(data) == nrow(stage$statistic), all(data$source_accver == reference$accver),
      all(data$source_start == data$window_start), all(data$source_end == data$window_end),
      all(data$sequence_unit == "bp"), all(data$coverage == data$covered / data$denominator),
      all(data$mean_depth == data$depth_sum / data$denominator))
    data$comparison <- comparison(data, id)
    data$source_hsp_id <- lapply(data$source_rows, function(i) sort(stage$source$source_hsp_id[i]))
    data[c(window_fields, "source_hsp_id")]
  })) |>
    arrange(match(comparison, groups), window_start)
  stopifnot(isTRUE(all.equal(rows[window_fields], coverage_rows[window_fields], check.attributes = FALSE)),
    identical(unclass(rows$source_hsp_id), unclass(coverage_rows$source_hsp_id)),
    all(vapply(out$track_registry$limits, function(x) all(x == c(0, 1)), logical(1))))
  invisible(out)
}
check_coverage(separate_tracks(groups, 5), function(data, id) id)
check_coverage(separate_tracks(rev(groups), 8), function(data, id) id)
check_coverage(grouped, function(data, id) data$colour)
fixed_result <- check_coverage(fixed, function(data, id) data$comparison)
local_result <- export_ggchord_layout(faceted +
  scale_track_value_continuous("coverage", facet_domain = "panel"))
panel_domains <- local_result$track_registry |>
  transmute(comparison, lower = vapply(limits, min, numeric(1)),
    upper = vapply(limits, max, numeric(1)))
expected_upper <- tapply(coverage_rows$covered / coverage_rows$denominator, coverage_rows$comparison, max)
stopifnot(identical(local_result$track_statistics$faceted, fixed_result$track_statistics$faceted),
  all(panel_domains$lower == 0), all(panel_domains$upper == expected_upper[panel_domains$comparison]))

first_last <- coverage_rows |>
  filter(window_start == 1 | window_end == reference$length) |>
  mutate(coverage = covered / denominator, mean_depth = depth_sum / denominator)
stopifnot(identical(first_last$covered, c(0, 0, 637, 3033, 0, 0)))

free <- tryCatch({
  export_ggchord_layout(fixed + facet_wrap(vars(comparison), scales = "free_y"))
  NA_character_
}, error = conditionMessage)
stopifnot(grepl("fixed scales", free, fixed = TRUE))
single_window <- ggplot(coverage_rows,
  aes(factor(window_start), covered / denominator,
      group = interaction(comparison, window_start))) + geom_boxplot()
boxes <- ggplot_build(single_window)$data[[1L]]
stopifnot(nrow(boxes) == 39L, all(boxes$ymin == boxes$lower & boxes$lower == boxes$middle &
  boxes$middle == boxes$upper & boxes$upper == boxes$ymax))

source_manifest <- tools::md5sum("data/flagship_dotplot_example.rda")
result <- list(source_manifest = source_manifest, midpoints = midpoints, cells = cells,
  first_endpoint_cells = first_cells, point_check = point_check, last_cells = last_cells,
  axes = projected, first_last = first_last, coverage = coverage_rows, panel_domains = panel_domains,
  free_error = free,
  source_hsp = links[c(source_fields, ".source_file", ".source_line")])
if (length(args)) saveRDS(result, args[[1L]])
cat(nrow(links), "HSPs;", nrow(cells), "midpoint cells;", nrow(first_cells), "first-endpoint cells.\n")
print(point_check, row.names = FALSE)
print(last_cells, row.names = FALSE)
print(peak)
print(first_last)
print(panel_domains)
cat(free, "\n")
cat(nrow(boxes), "singleton boxplots with coincident whiskers and quartiles.\n")
