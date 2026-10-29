# A small synthetic coordinate example, not biological source records.
library(ggchord)
sequences <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
windows <- data.frame(accver = c("A", "A", "B", "B"),
  start = c(100, 700, 100, 700), end = c(300, 900, 300, 900),
  locus_id = c("A_left", "A_right", "B_left", "B_right"),
  row_id = c("A", "A", "B", "B"), orientation = c(1, -1, 1, 1))
links <- data.frame(qaccver = "A", saccver = "B", qstart = 100,
  qend = 900, sstart = 100, send = 900, pident = 95)
features <- data.frame(accver = c("A", "A", "B", "B"),
  start = c(130, 730, 130, 730), end = c(260, 860, 260, 860),
  strand = c("+", "+", "+", "+"), anno = c("A left", "A reversed", "B left", "B right"))
view <- segment_ggchord_data(sequences, links, features, windows)
segmented_plot <- ggchord(view$seq_data, view$ribbon_data, view$gene_data) +
  geom_seq() + geom_link_ribbon(fill = "#83AEB8") + geom_gene() +
  geom_gene_label_repel() +
  scale_seq_position_continuous(coordinates = "source", breaks = c(100, 200, 300, 700, 800, 900)) +
  coord_collinear(segments = view$segments, segment_gap = 80) +
  labs(title = "Independent windows with a local reversal")
segmented_plot
# Each window has its own reversible source mapping. Gaps have no mapping.
map_ggchord_position(segmented_plot, "A_right", 800, source = TRUE)
