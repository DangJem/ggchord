# Visualize the explicit subset produced by prototype-graph-paths.py.
pkgload::load_all(quiet = TRUE)
library(dplyr)
seqs <- read.csv('tools/long-term-audit/v0.27.0-graph-path-sources.csv',
  colClasses = c('character', 'numeric', 'character'))
steps <- read.csv('tools/long-term-audit/v0.27.0-graph-path-steps.csv',
  colClasses = c('character', 'numeric', 'numeric', 'character', 'character', 'integer', 'integer'))
stopifnot(nrow(steps) == 200L, nrow(seqs) == 13L, !anyDuplicated(seqs$accver),
  all(steps$end >= steps$start), all(steps$end <= seqs$length[match(steps$accver, seqs$accver)]))
for (id in unique(steps$accver)) {
  d <- steps |> filter(accver == id)
  stopifnot(identical(d$source_step, 1:100), d$start[1L] == 1,
    all(tail(d$start, -1L) == head(d$end, -1L)+1))
}
windows <- steps |> group_by(accver) |> summarise(start = 1, end = max(end), .groups = 'drop')
focused <- focus_ggchord_data(seqs, feature_data = steps, loci = windows)
plot <- ggchord(focused$seq_data, feature_data = focused$feature_data) +
  geom_seq(show.legend = FALSE) + geom_seq_label(size = 2.8) +
  geom_feature(aes(feature_fill = strand), shape = 'block') + coord_collinear() +
  labs(title = 'Existing human LPA graph paths: explicit first 100 steps',
    subtitle = 'Path-local coordinates; node traversal orientation; no inferred pair links') +
  theme(plot.title.position = 'plot', plot.subtitle = element_text(hjust = .5))
out <- export_ggchord_layout(plot, original_data = TRUE)
ggplot2::ggsave('tools/long-term-audit/v0.27.0-graph-paths.png', plot,
  width = 10, height = 5, dpi = 110, bg = 'white')
cat('Graph path subset build: passed; original steps and complete supplied path lengths retained\n')
