# Direct checks for the source-space track grammar; no package/device matrix.
pkgload::load_all(".", quiet = TRUE)
s <- data.frame(accver = "A", length = 1000)
d <- data.frame(accver = "A", pos = c(100, 200, 400), value = c(2, 3, 4),
  low = c(1, 2, 3), high = c(3, 4, 5), label = letters[1:3])
results <- list()
for (geom in c("point", "line", "path", "step", "col", "bar", "tile", "area",
    "ribbon", "segment", "linerange", "errorbar", "text", "rect")) {
  mapping <- aes(x = pos, y = value)
  if (geom == "bar") mapping <- aes(x = pos)
  if (geom %in% c("ribbon", "linerange", "errorbar")) mapping <- aes(x = pos, ymin = low, ymax = high)
  if (geom == "segment") mapping <- aes(x = pos, y = value, xend = pos + 40, yend = high)
  if (geom == "text") mapping <- aes(x = pos, y = value, label = label)
  if (geom == "rect") mapping <- aes(xmin = pos - 20, xmax = pos + 20, ymin = low, ymax = high)
  p <- ggchord(s) + geom_seq() + geom_track(mapping, data = d, geom = geom,
    track = position_track("signal"), position = "identity") + coord_collinear()
  ggplot2::ggplotGrob(p)
  results[[geom]] <- "built and rendered to grob in collinear source frame"
}
events <- data.frame(accver = c("A", "A", "A", "B"), pos = c(100, 100, 300, 100), weight = c(2, 3, 1, 7))
seqs <- data.frame(accver = c("A", "B"), length = 1000)
p <- ggchord(seqs) + geom_seq() + geom_track(aes(x = pos, weight = weight,
  y = after_stat(count), fill = after_stat(count)), data = events,
  geom = "bar", track = position_track("count")) + coord_collinear()
e <- export_ggchord_layout(p, include = "track", original_data = TRUE)
computed <- e$track_statistics[[1]]
stopifnot(identical(computed$count, c(5, 1, 7)),
  identical(computed$source_rows[[1]], 1:2), !is.null(e$track_diagnostics))
invisible(ggplot2::ggplotGrob(p))
results$count <- "weighted counts match supplied source groups without crossing accessions"

overlay <- ggchord(s) + geom_seq() +
  geom_track(aes(x = pos, y = value), data = d, geom = "line", track = position_track("signal")) +
  geom_track(aes(x = pos, y = value), data = d, geom = "point", track = position_track("signal")) +
  coord_circular()
e <- export_ggchord_layout(overlay, include = "track")
stopifnot(nrow(e$track_registry) == 1L)
invisible(ggplot2::ggplotGrob(overlay))
results$overlay <- "line and point share one circular band/domain"

columns <- data.frame(accver = "A", start = c(1, 101), end = c(100, 200), value = c(.3, .7))
p <- ggchord(s) + geom_seq() + geom_track_col(aes(track_value = value), data = columns,
  position = position_track("signal")) + coord_chord()
e <- export_ggchord_layout(p, include = "track", original_data = TRUE)
stopifnot(identical(e$track_statistics[[1]]$track_value, columns$value),
  identical(e$track_statistics[[1]]$source_rows, I(list(1L, 2L))))
results$col <- "supplied values and full interval widths retain identity semantics"

p <- ggchord(s) + geom_seq() + geom_track(aes(x = pos, y = after_stat(count)),
  data = events[1:3, ], geom = "col", stat = "bin", binwidth = 100,
  track = position_track("bins"), bounds = "clip") + coord_collinear()
e <- export_ggchord_layout(p, include = "track")
stopifnot(identical(e$track_statistics[[1]]$count, c(2, 0, 1)))
results$bin <- "explicit binwidth; empty bin retained and membership recorded"
p <- ggchord(s) + geom_seq() + geom_track(aes(x = pos, y = weight),
  data = events[1:3, ], geom = "point", stat = "summary", fun = mean,
  track = position_track("mean")) + coord_collinear()
e <- export_ggchord_layout(p, include = "track")
stopifnot(identical(e$track_statistics[[1]]$track_value, c(2.5, 1)))
results$summary <- "native summary operates in source positions and value units"
jsonlite::write_json(results, "tools/long-term-audit/track-grammar-results-v0.17.0.json",
  auto_unbox = TRUE, pretty = TRUE)
cat("Track grammar checks passed\n")
