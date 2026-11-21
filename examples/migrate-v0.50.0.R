# Historical example for the first v0.50.0 API design at commit d6faf312.
# The short intervals below are teaching data in bp.
library(ggchord)
library(ggplot2)
library(dplyr)

reference <- tibble(accver = "A", length = 10L, sequence_unit = "bp")
intervals <- tibble(
  accver = "A", sample = c("S1", "S1", "S2"),
  feature_id = c("f1", "f2", "f3"), strand = "+",
  start = c(3L, 4L, 10L), end = c(6L, 5L, 10L)
)

# Each sample uses the same physical band, source windows and numerical scale.
coverage_plot <- ggchord(reference) +
  geom_seq(id = "reference") +
  track_band("coverage", width = grid::unit(8, "mm")) +
  stat_coverage(data = intervals, binwidth = 4, track = "coverage", id = "windows") +
  geom_track_line(
    aes(y = after_stat(coverage)), data = intervals, stat = "coverage",
    binwidth = 4, track = "coverage", colour = "#216579"
  ) +
  scale_track_value_continuous("coverage", limits = c(0, 1)) +
  facet_wrap(vars(sample)) +
  coord_collinear()

window_data <- export_ggchord_layout(coverage_plot)$track_statistics$windows |>
  select(sample, window_start, window_end, denominator, covered, coverage, source_rows)

# A filtered label layer follows the lanes computed from all the feature rows.
captions <- intervals |>
  filter(feature_id %in% c("f2", "f3"))
feature_plot <- ggchord(reference) +
  geom_seq() +
  geom_feature(data = intervals, id = "features", position = position_feature_stack()) +
  geom_feature_label(
    aes(label = feature_id), data = captions,
    anchor_layer = "features", placement = "repel"
  ) +
  facet_wrap(vars(sample)) +
  coord_collinear()

print(window_data)
print(coverage_plot)
print(feature_plot)

# A site occupies one base. A threshold belongs to the same numerical scale
# as the observed coverage and spans only the supplied source positions.
samples <- distinct(intervals, accver, sample)
sites <- mutate(samples, start = 4L, end = 4L)
threshold <- bind_rows(mutate(samples, pos = 1L), mutate(samples, pos = 10L)) |>
  mutate(fraction = 0.5)
marked_plot <- coverage_plot +
  geom_seq_region(data = sites) +
  geom_track_line(aes(x = pos, y = fraction), data = threshold,
    measure = "coverage", track = "coverage", linetype = "dashed")
print(marked_plot)

# Independent samples can also be compared by the colour of each whole window.
tile_plot <- ggchord(reference) +
  geom_seq() +
  stat_coverage(aes(fill = after_stat(coverage)), data = intervals,
    binwidth = 4, geom = "track_tile", id = "tiles") +
  facet_wrap(vars(sample)) +
  coord_collinear()
print(tile_plot)
