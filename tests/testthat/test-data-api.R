test_that("validation, cleaning and one malformed input have public outcomes", {
  seq <- smoke_sequences()
  ribbon <- smoke_ribbon()
  genes <- smoke_genes()
  expect_true(validate_ggchord_data(seq, ribbon, genes)$valid)
  cleaned <- clean_ggchord_data(seq, ribbon, transform(genes, accver = "unknown"),
    unknown_id = "keep")
  expect_s3_class(cleaned, "ggchord_clean")
  expect_equal(nrow(cleaned$ribbon_data), 1L)
  expect_equal(nrow(cleaned$gene_data), nrow(genes))
  expect_equal(cleaned$ribbon_data$direction, "reverse")
  expect_error(ggchord(transform(seq, length = "bad")), "numeric")
})

test_that("import helpers parse minimal files with explicit coordinate conventions", {
  fasta <- tempfile(fileext = ".fna")
  blast <- tempfile(fileext = ".tsv")
  gff <- tempfile(fileext = ".gff3")
  bed <- tempfile(fileext = ".bed")
  paf <- tempfile(fileext = ".paf")
  on.exit(unlink(c(fasta, blast, gff, bed, paf)), add = TRUE)
  writeLines(c(">A", "ACGTACGT"), fasta)
  writeLines(paste("A", "B", 98, 100, 0, 0, 1, 100, 1, 100, "1e-20", 200,
    sep = "\t"), blast)
  writeLines(c("##gff-version 3",
    "A\tsource\tCDS\t1\t8\t.\t+\t0\tID=cds1;product=test%20protein",
    "##FASTA", ">A", "ACGTACGT"), gff)
  expect_equal(read_fasta_lengths(fasta)$length, 8L)
  expect_equal(nrow(read_blast(blast, format = "outfmt6")), 1L)
  expect_equal(read_gff3(gff)$anno, "test protein")
  writeLines("A\t0\t1\tpoint\t0\t.", bed)
  expect_equal(unname(unlist(read_bed(bed)[c("start", "end")])), c(1, 1))
  writeLines("A\t1000\t99\t200\t-\tB\t1000\t299\t400\t90\t101\t255", paf)
  mapping <- read_paf(paf)
  expect_equal(unname(unlist(mapping[c("qstart", "qend", "sstart", "send")])),
    c(100, 200, 400, 300))
})

test_that("ribbon filter, deduplication and merge each have one success path", {
  ribbon <- transform(smoke_ribbon(), pident = 95, length = 101)
  expect_equal(nrow(filter_ggchord_ribbons(ribbon, min_pident = 90)$data), 1L)
  expect_equal(nrow(deduplicate_ggchord_ribbons(rbind(ribbon, ribbon))$data), 1L)
  blocks <- rbind(ribbon, transform(ribbon,
    qstart = 201, qend = 301, sstart = 299, send = 199))
  expect_equal(nrow(merge_ggchord_ribbons(blocks)$data), 1L)
})

test_that("explicit bundling and layout suggestions preserve sequence membership", {
  ribbon <- transform(smoke_ribbon(), pident = 95, length = 101)
  bundled <- bundle_ggchord_ribbons(rbind(ribbon, ribbon),
    smoke_sequences(), bins = 10)
  expect_equal(nrow(bundled$data), 1L)
  expect_equal(bundled$data$.bundle_n, 2L)
  suggestion <- optimize_ggchord_layout(smoke_sequences(), ribbon)
  expect_setequal(suggestion$seq_order, smoke_sequences()$accver)
})

test_that("focus composes nested segment coordinates and keeps reverse links", {
  features <- smoke_genes()[1, ]
  features$segments <- I(list(data.frame(start = c(100, 160), end = c(140, 200),
    segment_index = 1:2, segment_type = "standard", length_bp = 41)))
  focused <- focus_ggchord_data(smoke_sequences(), smoke_ribbon(), features,
    loci = data.frame(accver = c("A", "B"), start = c(120, 320),
      end = c(180, 380)))
  again <- focus_ggchord_data(focused$seq_data, focused$ribbon_data,
    focused$gene_data, loci = data.frame(accver = c("A", "B"), start = 2, end = 60))
  expect_true(all(again$gene_data$segments[[1]]$start >= 1))
  expect_true(all(again$gene_data$segments[[1]]$end <= again$seq_data$length[1]))
  expect_equal(again$seq_data$.source_origin, c(121, 321))
  expect_gt(again$ribbon_data$sstart, again$ribbon_data$send)
  expect_equal(again$gene_data$.source_start, features$start)
})


test_that("preparation and target loci compose without changing source identities", {
  sequences <- transform(smoke_sequences(), genome = c("isolate1", "isolate2"))
  genes <- transform(smoke_genes(), feature_id = c("gA", "gB"))
  hierarchy <- validate_ggchord_hierarchy(sequences,
    data.frame(genome = c("isolate1", "isolate2")))
  expect_true(hierarchy$valid)
  homologs <- link_ggchord_clusters(genes,
    data.frame(feature_id = c("gA", "gB"), cluster_id = "family"),
    data.frame(qaccver = "A", saccver = "B"))
  expect_equal(nrow(homologs$data), 1L)
  expect_false(any(c("pident", "length") %in% names(homologs$data)))
  prepared <- prepare_ggchord_data(sequences,
    rbind(smoke_ribbon(), smoke_ribbon()), features = list(genes = genes),
    steps = list(deduplicate = list(keep = "first")))
  expect_true(prepared$report$valid)
  expect_equal(nrow(prepared$ribbon_data), 1L)
  expect_equal(prepared$features$genes$feature_id, genes$feature_id)
  loci <- find_ggchord_loci(prepared$seq_data, prepared$features$genes,
    targets = "gB", upstream = 50, downstream = 20)
  expect_equal(unname(unlist(loci[c("start", "end")])), c(380, 550))
  focused <- focus_ggchord_data(prepared$seq_data,
    gene_data = prepared$features$genes, loci = loci)
  expect_equal(focused$gene_data$feature_id, "gB")
  expect_equal(focused$seq_data$genome, "isolate2")
  opposite <- rbind(smoke_ribbon(),
    transform(smoke_ribbon(), sstart = 300, send = 400))
  bundled <- prepare_ggchord_data(sequences, opposite,
    steps = list(clean = list(), bundle = list(bins = 1, weight = "count")))
  expect_equal(nrow(bundled$ribbon_data), 2L)
  expect_setequal(bundled$ribbon_data$direction, c("same", "reverse"))
})

test_that("spliced CDS domains project and retain their parent identity", {
  coding <- data.frame(feature_id = "cds", accver = "A", start = c(100, 200),
    end = c(103, 204), strand = "+", phase = c(0, 2))
  cds <- build_ggchord_cds(coding)
  domain <- data.frame(feature_id = "domain", parent_id = "cds", start = 2, end = 3)
  genomic <- project_ggchord_features(domain, cds)
  expect_equal(genomic$parent_id, "cds")
  expect_equal(nrow(genomic$segments[[1L]]), 2L)
  protein <- project_ggchord_features(genomic, cds, from = "genome")
  expect_equal(protein[c("feature_id", "start", "end")], domain[c("feature_id", "start", "end")])
})

test_that("cleaning retains one-base inclusive gene and ribbon intervals", {
  ribbon <- transform(smoke_ribbon(), qend = qstart, send = sstart)
  gene <- transform(smoke_genes()[1, ], end = start)
  cleaned <- clean_ggchord_data(smoke_sequences(), ribbon, gene)
  expect_equal(nrow(cleaned$ribbon_data), 1L)
  expect_equal(nrow(cleaned$gene_data), 1L)
  expect_equal(cleaned$ribbon_data$qstart, cleaned$ribbon_data$qend)
  expect_equal(cleaned$gene_data$start, cleaned$gene_data$end)
})
