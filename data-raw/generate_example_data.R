# Rebuild the packaged example data from the unchanged files in examples/.
# Run from the package root with:
#   Rscript data-raw/generate_example_data.R

pkgload::load_all(".", quiet = TRUE)

seq_data_example <- utils::read.delim(
  "examples/seq_track.tsv", stringsAsFactors = FALSE, check.names = FALSE
)
raw_genes <- utils::read.delim(
  "examples/gene_track.tsv", stringsAsFactors = FALSE, check.names = FALSE
)

# Choose a small, readable set that spans each complete sequence. Informative
# annotations receive a modest preference over repeated "hypothetical
# protein" entries, but genomic coverage remains the dominant criterion.
select_demo_genes <- function(data, sequence_data, n_per_sequence = 8L) {
  selected <- vector("list", nrow(sequence_data))
  for (i in seq_len(nrow(sequence_data))) {
    sid <- sequence_data$seq_id[i]
    sequence_length <- sequence_data$length[i]
    candidates <- data[data$seq_id == sid, , drop = FALSE]
    candidates$.source_order <- seq_len(nrow(candidates))
    candidates$.midpoint <- (candidates$start + candidates$end) / 2
    candidates$.feature_length <- abs(candidates$end - candidates$start) + 1

    # Very short features disappear at normal example sizes. Keep them in the
    # raw file, but prefer features that can be inspected in the package demo.
    visible <- candidates$.feature_length >= max(700, sequence_length * 0.01)
    if (sum(visible) >= n_per_sequence) candidates <- candidates[visible, ]

    targets <- seq(0.08, 0.92, length.out = n_per_sequence) * sequence_length
    chosen <- integer(0)
    for (target in targets) {
      available <- setdiff(seq_len(nrow(candidates)), chosen)
      distance <- abs(candidates$.midpoint[available] - target) / sequence_length
      informative <- candidates$anno[available] != "hypothetical protein"
      visibility <- pmin(candidates$.feature_length[available] / sequence_length, 0.03)
      score <- distance - 0.035 * informative - 0.15 * visibility
      chosen <- c(chosen, available[order(
        score, candidates$.source_order[available]
      )[1L]])
    }

    demo <- candidates[chosen, c(
      "seq_id", "start", "end", "strand", "anno"
    ), drop = FALSE]
    demo <- demo[order(demo$start, demo$end), , drop = FALSE]
    demo$source_strand <- demo$strand

    # The source annotations are strongly strand-biased (two sequences contain
    # no reverse-strand records). The packaged object is a plotting fixture, so
    # alternate its display strand while retaining source_strand explicitly.
    pattern <- if (i %% 2L) c("+", "-") else c("-", "+")
    demo$strand <- rep(pattern, length.out = nrow(demo))
    selected[[i]] <- demo
  }

  out <- do.call(rbind, selected)
  row.names(out) <- NULL
  out
}

gene_data_example <- select_demo_genes(raw_genes, seq_data_example)

blast_files <- sort(list.files(
  "examples/blastn", pattern = "\\.o7$", full.names = TRUE
))
raw_ribbons <- read_blast(files = blast_files, format = "outfmt7")

# Keep sparse sequence pairs intact while preventing one dense comparison from
# dominating the teaching figure. Within a dense pair, deterministic
# farthest-point sampling retains long, high-identity alignments distributed
# across both sequences instead of several nearly coincident short ribbons.
select_demo_ribbons <- function(data, max_per_pair = 3L) {
  data <- data[data$length >= 300, , drop = FALSE]
  data$.source_order <- seq_len(nrow(data))
  data$.pair <- paste(
    pmin(data$qaccver, data$saccver),
    pmax(data$qaccver, data$saccver),
    sep = "\r"
  )
  groups <- split(seq_len(nrow(data)), data$.pair, drop = TRUE)
  selected <- integer(0)

  for (rows in groups) {
    if (length(rows) <= max_per_pair) {
      selected <- c(selected, rows)
      next
    }
    qmid <- (data$qstart[rows] + data$qend[rows]) / 2 / data$qlen[rows]
    smid <- (data$sstart[rows] + data$send[rows]) / 2 / data$slen[rows]
    length_score <- log1p(data$length[rows])
    length_score <- (length_score - min(length_score)) /
      max(diff(range(length_score)), 1)
    identity_score <- data$pident[rows] / 100
    utility <- 0.12 * length_score + 0.03 * identity_score
    chosen <- which.max(utility)

    while (length(chosen) < max_per_pair) {
      available <- setdiff(seq_along(rows), chosen)
      separation <- vapply(available, function(i) {
        min(sqrt(
          (qmid[i] - qmid[chosen])^2 +
            (smid[i] - smid[chosen])^2
        ))
      }, numeric(1))
      score <- separation + utility[available]
      chosen <- c(chosen, available[order(
        -score, data$.source_order[rows[available]]
      )[1L]])
    }
    selected <- c(selected, rows[chosen])
  }

  out <- data[sort(selected), setdiff(names(data), c(
    ".source_order", ".pair"
  )), drop = FALSE]
  row.names(out) <- NULL
  out
}

ribbon_data_example <- select_demo_ribbons(raw_ribbons)
row.names(ribbon_data_example) <- NULL

save(seq_data_example, file = "data/seq_data_example.rda", compress = "xz")
save(ribbon_data_example, file = "data/ribbon_data_example.rda", compress = "xz")
save(gene_data_example, file = "data/gene_data_example.rda", compress = "xz")
