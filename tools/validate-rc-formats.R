# One representative flexible-parameter composition, not equivalent-form enumeration.
pkgload::load_all(quiet = TRUE)
local({
  grDevices::pdf(NULL, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  seqs <- seq_data_example[c("accver", "length")]
  genes <- transform(gene_data_example, strand = source_strand)[c("accver", "start", "end", "strand")]
  radii <- setNames(c(.9, 1.1, 1, 1.2), seqs$accver)
  widths <- setNames(list(c("+" = .04, "-" = .06)), seqs$accver[1L])
  offsets <- setNames(list(c("+" = .1, "-" = -.2)), seqs$accver[1L])
  plot <- ggchord(seqs) + geom_seq(seq_order = rev(seqs$accver), seq_radius = radii) +
    geom_gene(data = genes, width = widths, position = position_strand(offsets))
  layout <- get_chord_layout(plot)
  stopifnot(identical(layout$seqs, rev(seqs$accver)),
    identical(unname(layout$seqRadius[seqs$accver]), unname(radii)))
  out <- export_ggchord_layout(plot, include = "gene")
  stopifnot(setequal(out$gene$source_row, seq_len(nrow(genes))))
  write.csv(data.frame(case = "named_reordering_and_nested_strands", status = "PASS_FINITE",
    detail = "Named radius follows accver across reversal; one nested strand width/offset composes with all true gene rows; other forms source-reviewed, not enumerated",
    source = "R/helpers_params.R; R/position-feature.R; tools/validate-rc-formats.R",
    R = R.version.string, ggplot2 = as.character(utils::packageVersion("ggplot2"))),
    "design-roadmap/audit-generated/rc-format-validation-v0.16.6.csv", row.names = FALSE)
  cat("One reordered named-vector/nested-strand composition passed.\n")
})
