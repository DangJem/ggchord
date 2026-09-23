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

# Circular-plasmid fixtures. Every checked-in FASTA record under
# examples/plasmid becomes one plasmid_example_* object. The file stem is the
# stable sequence ID and display label; FASTA headers remain source metadata
# and are not treated as identifiers because several contain spaces.
read_one_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines) || !grepl("^>", lines[1L])) {
    stop("Invalid FASTA record: ", path, call. = FALSE)
  }
  sequence <- toupper(paste0(lines[!grepl("^>", lines)], collapse = ""))
  if (!nzchar(sequence) || !grepl("^[ACGTRYSWKMBDHVN]+$", sequence)) {
    stop("Invalid DNA sequence in ", path, call. = FALSE)
  }
  list(
    header = sub("^>", "", lines[1L]),
    sequence = sequence
  )
}

plasmid_files <- sort(list.files(
  "examples/plasmid", pattern = "\\.fna$", full.names = TRUE
))
if (!length(plasmid_files)) {
  stop("No plasmid FASTA records found under examples/plasmid", call. = FALSE)
}
plasmid_labels <- tools::file_path_sans_ext(basename(plasmid_files))
names(plasmid_files) <- plasmid_labels
plasmids <- lapply(plasmid_files, read_one_fasta)

plasmid_object_name <- function(label) {
  suffix <- gsub("\\(\\+\\)", "_plus", label)
  suffix <- gsub("[^[:alnum:]]+", "_", suffix)
  suffix <- gsub("^_+|_+$", "", suffix)
  paste0("plasmid_example_", suffix)
}
plasmid_object_names <- vapply(
  plasmid_labels, plasmid_object_name, character(1L)
)
if (anyDuplicated(plasmid_object_names)) {
  stop("Plasmid file names do not produce unique example object names",
    call. = FALSE)
}

make_plasmid_object <- function(key) {
  record <- plasmids[[key]]
  data.frame(
    accver = key,
    label = key,
    length = nchar(record$sequence),
    sequence = record$sequence,
    stringsAsFactors = FALSE
  )
}

plasmid_examples <- stats::setNames(
  lapply(plasmid_labels, make_plasmid_object), plasmid_object_names
)
invisible(list2env(plasmid_examples, envir = environment()))

# Retain the pre-release fixture spelling as a compatibility alias while the
# visual benchmark moves to the official pUC19 reference.
plasmid_example_pUC19c <- plasmid_example_pUC19
save(plasmid_example_pUC19c,
  file = "data/plasmid_example_pUC19c.rda", compress = "xz")
for (object_name in names(plasmid_examples)) {
  save(list = object_name,
    file = file.path("data", paste0(object_name, ".rda")), compress = "xz")
}
