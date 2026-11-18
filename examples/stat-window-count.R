# A third-party ggplot2 Stat for ggchord source-coordinate tracks.
# Each supplied closed window counts the intervals that touch it.
library(ggchord)
library(ggplot2)

StatWindowCount <- ggproto("StatWindowCount", Stat,
  required_aes = c("x", "xmin", "xmax"),
  dropped_aes = c("xmin", "xmax", "source_row"),
  compute_group = function(data, scales, windows) {
    members <- lapply(seq_len(nrow(windows)), function(i) {
      selected <- data$xmin <= windows$end[i] & data$xmax >= windows$start[i]
      data$source_row[selected]
    })
    data.frame(
      # The display edges of [start, end] are start - 1 and end.
      x = (windows$start - 1 + windows$end) / 2,
      n_interval = lengths(members),
      source_rows = I(members)
    )
  }
)

reference <- data.frame(accver = "A", length = 12L)
intervals <- data.frame(
  accver = "A", sample = c("first", "second", "first"),
  start = c(2L, 5L, 4L), end = c(5L, 6L, 4L)
)
windows <- data.frame(start = c(1L, 5L, 9L), end = c(4L, 8L, 12L))
counts <- aes(x = start, xmin = start, xmax = end, y = after_stat(n_interval))

window_plot <- ggchord(reference) +
  geom_track_col(counts, data = intervals, stat = StatWindowCount,
    windows = windows, width = 4, position = "identity", id = "bars") +
  geom_track_line(counts, data = intervals, stat = StatWindowCount,
    windows = windows, id = "line") +
  facet_wrap(vars(sample)) +
  coord_collinear()
window_result <- export_ggchord_layout(window_plot, data_stage = "statistic")
fields <- c("sample", "pos", "n_interval", "track_value", "source_rows")
bar_counts <- window_result$layer_data$bars$statistic[fields]
line_counts <- window_result$layer_data$line$statistic[fields]
stopifnot(identical(bar_counts, line_counts))
# The two samples have counts (2, 1, 0) and (0, 1, 0) across the three windows.
bar_counts
