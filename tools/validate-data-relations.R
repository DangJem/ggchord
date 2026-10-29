# Run explicitly from repository root; no testthat or package check.
pkgload::load_all(quiet = TRUE)

sequences <- data.frame(accver = c("A", "B", "C"), length = c(1000, 1000, 1000),
  genome = c("isolate1", "isolate1", "isolate2"))
genomes <- data.frame(genome = c("isolate1", "isolate2", "unused"), condition = c("a", "b", "c"))
hierarchy <- validate_ggchord_hierarchy(sequences, genomes)
stopifnot(hierarchy$valid, identical(hierarchy$seq_data, sequences),
  identical(hierarchy$genomes, genomes), identical(hierarchy$contigs$n_contigs, c(2L, 1L, 0L)),
  hierarchy$issues$category == "orphan_genome",
  identical(hierarchy$contigs$sequence_rows[[1L]], 1:2))
broken <- sequences
broken$accver[2L] <- "A"
broken$genome[1L] <- "missing"
broken$genome[3L] <- NA_character_
invalid <- validate_ggchord_hierarchy(broken, genomes)
stopifnot(!invalid$valid, all(c("duplicate_key", "unknown_genome", "missing_key") %in% invalid$issues$category))

features <- data.frame(feature_id = c("a1", "a2", "b1", "c1"), accver = c("A", "A", "B", "C"),
  start = c(10, 100, 300, 500), end = c(30, 130, 340, 540), strand = c("+", "+", "-", "."))
membership <- data.frame(feature_id = c("a1", "a2", "b1", "c1", "missing"),
  cluster_id = c("group1", "group1", "group1", "group1", "unknown"), evidence = "input assertion")
pairs <- data.frame(qaccver = c("A", "B"), saccver = c("B", "C"))
strict <- link_ggchord_clusters(features, membership, pairs)
stopifnot(nrow(strict$data) == 1L, !strict$report$valid,
  "multiple_members" %in% strict$report$comparisons$reason,
  "unknown_feature" %in% strict$report$members$.status,
  !strict$data$orientation_known, !any(c("length", "pident") %in% names(strict$data)))
all <- link_ggchord_clusters(features, membership, pairs, paralogs = "all", max_links = 3)
stopifnot(nrow(all$data) == 3L, all$data$qstart[1L] == 10,
  all$data$sstart[1L] == 340, all$data$send[1L] == 300,
  all$data$direction[1L] == "reverse", is.na(tail(all$data$direction, 1L)),
  identical(all$report$members$evidence, membership$evidence),
  identical(all$data$.query_feature_rows[[1L]], 1L))
limit <- tryCatch(link_ggchord_clusters(features, membership, pairs,
  paralogs = "all", max_links = 2), error = identity)
stopifnot(inherits(limit, "error"), grepl("max_links", conditionMessage(limit)))

# Explicit discontinuous feature spans retain actual segments and source rows.
discontinuous <- rbind(features[1L, ], transform(features[1L, ], start = 50, end = 70), features[3L, ])
span <- link_ggchord_clusters(discontinuous, membership[c(1L, 3L), ], pairs[1L, ], segment_policy = "span")
stopifnot(span$data$qstart == 10, span$data$qend == 70,
  identical(span$data$.query_feature_rows[[1L]], 1:2),
  identical(span$data$qsegments[[1L]]$start, c(10, 50)))
prepared <- prepare_ggchord_data(sequences, span$data,
  features = list(genes = discontinuous))
stopifnot(prepared$report$valid, nrow(prepared$ribbon_data) == 1L,
  !any(c("length", "pident") %in% names(prepared$ribbon_data)))
cat("Explicit hierarchy, cluster links, ambiguity and expansion-limit checks passed.\n")
