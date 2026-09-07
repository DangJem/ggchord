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

names(seq_data_example)[names(seq_data_example) == "seq_id"] <- "accver"
names(raw_genes)[names(raw_genes) == "seq_id"] <- "accver"

# Choose a small, readable set that spans each complete sequence. Informative
# annotations receive a modest preference over repeated "hypothetical
# protein" entries, but genomic coverage remains the dominant criterion.
select_demo_genes <- function(data, sequence_data, n_per_sequence = 8L) {
  selected <- vector("list", nrow(sequence_data))
  for (i in seq_len(nrow(sequence_data))) {
    sid <- sequence_data$accver[i]
    sequence_length <- sequence_data$length[i]
    candidates <- data[data$accver == sid, , drop = FALSE]
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
      "accver", "start", "end", "strand", "anno"
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

# Minimal single-genome fixture. The source FASTA deliberately contains a few
# common restriction motifs and is kept intact under examples/single-genome/.
single_lines <- readLines("examples/single-genome/minimal.fasta", warn = FALSE)
single_sequence <- paste0(single_lines[!grepl("^>", single_lines)], collapse = "")
single_genome_example <- data.frame(
  accver = "minimal_genome", length = nchar(single_sequence),
  stringsAsFactors = FALSE
)
single_gene_example <- utils::read.delim(
  "examples/single-genome/features.tsv", stringsAsFactors = FALSE
)
restriction_site_example <- find_restriction_sites(
  stats::setNames(single_sequence, "minimal_genome"),
  enzymes = c("EcoRI", "BamHI", "HindIII", "PstI", "SmaI")
)

save(single_genome_example, file = "data/single_genome_example.rda", compress = "xz")
save(single_gene_example, file = "data/single_gene_example.rda", compress = "xz")
save(restriction_site_example, file = "data/restriction_site_example.rda", compress = "xz")

# Circular-plasmid fixtures. FASTA records are read without rewriting them.
read_one_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  list(
    header = sub("^>", "", lines[1L]),
    sequence = paste0(lines[!grepl("^>", lines)], collapse = "")
  )
}
plasmid_files <- c(
  "pUC19c" = "examples/plasmid/pUC19c.fna",
  "pBR322" = "examples/plasmid/pBR322.fna",
  "pBluescript II SK(+)" = "examples/plasmid/pBluescript II SK(+).fna"
)
plasmids <- lapply(plasmid_files, read_one_fasta)
plasmid_sequence_example <- data.frame(
  accver = sub(" .*$", "", vapply(plasmids, `[[`, character(1), "header")),
  label = names(plasmids),
  length = nchar(vapply(plasmids, `[[`, character(1), "sequence")),
  sequence = vapply(plasmids, `[[`, character(1), "sequence"),
  stringsAsFactors = FALSE
)

# A compact, independently specified teaching panel of common motifs keeps the
# installed fixture reproducible without redistributing the complete REBASE
# database. Full REBASE parsing is audited separately by
# generate_rebase_database.R and remains license-gated.
common_motifs <- c(
  EcoRI = "GAATTC", BamHI = "GGATCC", HindIII = "AAGCTT",
  PstI = "CTGCAG", SmaI = "CCCGGG", KpnI = "GGTACC",
  SacI = "GAGCTC", SalI = "GTCGAC", XbaI = "TCTAGA",
  XhoI = "CTCGAG", SpeI = "ACTAGT", NotI = "GCGGCCGC",
  EagI = "CGGCCG", ApaI = "GGGCCC", EcoRV = "GATATC"
)
plasmid_restriction_example <- find_restriction_sites(
  stats::setNames(
    plasmid_sequence_example$sequence, plasmid_sequence_example$accver
  ),
  patterns = common_motifs
)

# Core annotations for examples. pBR322 coordinates follow GenBank J01749.1;
# pUC19c coordinates follow the feature summary embedded in GenBank L09137.2.
# The fixture is deliberately compact rather than a replacement for either
# source record's complete feature table.
plasmid_feature_example <- data.frame(
  accver = c(rep("J01749.1", 11), rep("L09137.2", 3)),
  start = c(
    27, 43, 86, 1515, 1788, 1905, 1915, 2011, 2351, 2535, 3293,
    238, 396, 1629
  ),
  end = c(
    33, 49, 1276, 1519, 1792, 1910, 2106, 2167, 2414, 2540, 4153,
    682, 452, 2417
  ),
  strand = c(
    "-", "+", "+", "-", "-", "+", "+", "+", "-", "+", "-",
    "-", "-", "-"
  ),
  type = c(
    "promoter", "promoter", "CDS", "repeat_region", "repeat_region",
    "RBS", "CDS", "misc_feature", "misc_feature", "rep_origin", "CDS",
    "CDS", "misc_feature", "CDS"
  ),
  anno = c(
    "P1", "P2", "tet", "direct repeat", "direct repeat", "RBS",
    "rop", "H-strand effector", "L-strand effector", "ori", "bla",
    "lacZ alpha", "MCS", "bla"
  ),
  source = c(rep("GenBank J01749.1", 11), rep("GenBank L09137.2", 3)),
  stringsAsFactors = FALSE
)

save(plasmid_sequence_example,
  file = "data/plasmid_sequence_example.rda", compress = "xz")
save(plasmid_feature_example,
  file = "data/plasmid_feature_example.rda", compress = "xz")
save(plasmid_restriction_example,
  file = "data/plasmid_restriction_example.rda", compress = "xz")
