# Synthetic numeric teaching inputs; no biological measurements are invented.
# For existing coverage, GC-content or GC-skew tables, map the supplied value
# column in the same way as signal below. Choose a meaningful unit and limits
# with scale_track_value_continuous(); signed GC-skew may use a zero baseline.
# The plotting layer does not calculate coverage or sequence composition.
library(ggchord)
sequences <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
bins <- data.frame(accver = rep(c("A", "B"), each = 20),
  start = rep(seq(1, 951, 50), 2), end = rep(seq(50, 1000, 50), 2),
  signal = c(sin(seq(0, 2 * pi, length.out = 20)), cos(seq(0, 2 * pi, length.out = 20))))
sites <- data.frame(accver = c("A", "B"), start = c(150, 500), end = c(350, 750),
  category = c("site", "reference"))
track_plot <- ggchord(sequences) + geom_seq() +
  geom_track_bar(aes(track_value = signal, track_fill = signal), data = bins,
    position = position_track("signal", width = unit(5, "mm"))) +
  scale_track_value_continuous("signal", limits = c(-1, 1), unit = "arbitrary units") +
  scale_track_fill_gradientn("signal", colours = c("#31537D", "#F3F4EE", "#D26B38")) +
  geom_track_tile(aes(track_fill = category), data = sites,
    position = position_track("sites", side = "inside", width = unit(2.5, "mm"))) +
  scale_track_fill_manual("sites", values = c(site = "#468C86", reference = "#D1AB59")) +
  scale_seq_position_continuous(sequences = "B") + geom_seq_scale(length = 100) +
  coord_collinear() + labs(title = "Explicit tracks with shared numerical domains")
track_plot
# Saving remains standard ggsave(). Diagnostics are explicitly requested:
# diagnose_ggchord(track_plot, width = 8, height = 5, units = "in")
# ggsave("tracks.pdf", track_plot, width = 8, height = 5)
