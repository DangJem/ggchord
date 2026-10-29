# Optional issue #2 workflow. These are declared teaching intervals, not a
# measurement of a named host/virus integration event.
library(ggchord)

if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
    !requireNamespace("InteractionSet", quietly = TRUE)) {
  message("This example requires the optional GenomicRanges and InteractionSet packages.")
} else {
  intervals <- GenomicRanges::GRanges(
    seqnames = c("host", "virus"),
    ranges = IRanges::IRanges(start = c(100, 10), end = c(180, 50)),
    strand = c("-", "*"),
    seqlengths = c(host = 1000L, virus = 300L)
  )
  S4Vectors::mcols(intervals)$label <- c("host annotation", "viral interval")
  interactions <- InteractionSet::GInteractions(intervals[1L], intervals[2L])
  S4Vectors::mcols(interactions)$support <- 12L

  # Explicit conversion is also available for downstream packages that need
  # plain tables. Source row, anchor identities and metadata are retained.
  seqs <- as_ggchord_seq(intervals)
  features <- as_ggchord_feature(intervals)
  links <- as_ggchord_link(interactions)

  # Direct S4 inputs follow the same conversion contract.
  # The supplied negative-strand anchor is intentional; its decreasing link
  # endpoints preserve direction, so avoid the generic reversed-interval hint.
  bioc_example <- ggchord(GenomicRanges::seqinfo(intervals), interactions,
    validate = "none") +
    geom_link_ribbon(aes(ribbon_alpha = support)) +
    scale_ribbon_alpha_continuous(range = c(.2, .6)) +
    geom_seq() + geom_seq_label() +
    geom_feature(data = intervals) +
    geom_feature_label_repel(data = intervals) +
    coord_collinear()
  if (interactive()) print(bioc_example)

  # Usual ggplot2 export, when requested:
  # ggsave("integration.png", bioc_example, width = 9, height = 5)
}
