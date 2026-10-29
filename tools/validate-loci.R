# Run explicitly from the repository root; no package test/check is performed.
pkgload::load_all(quiet = TRUE)

sequences <- data.frame(accver = c("A", "B"), length = c(1000, 600),
  genome = c("isolate1", "isolate2"))
features <- data.frame(feature_id = c("forward", "reverse", "split", "split", "unknown"),
  accver = c("A", "A", "B", "B", "B"), start = c(10, 800, 100, 200, 500),
  end = c(50, 950, 150, 250, 550), strand = c("+", "-", "+", "+", "."))
loci <- find_ggchord_loci(sequences, features, c("reverse", "forward", "split"),
  upstream = 100, downstream = 20)
stopifnot(identical(loci$target_id, c("reverse", "forward", "split")),
  identical(loci$start, c(780, 1, 1)), identical(loci$end, c(1000, 70, 270)),
  loci$clipped_right[1], loci$clipped_left[2],
  identical(loci$target_source_rows[[3]], c(3L, 4L)),
  identical(loci$genome, c("isolate1", "isolate1", "isolate2")))
focused <- focus_ggchord_data(sequences, gene_data = features, loci = loci)
stopifnot(identical(focused$seq_data$accver, loci$locus_id),
  identical(focused$seq_data$length, loci$end - loci$start + 1),
  sum(focused$seq_data$.source_accver == "A") == 2L)
unknown <- find_ggchord_loci(sequences, features, "unknown", 10, 10)
stopifnot(unknown$start == 490, unknown$end == 560, unknown$strand == ".")

expect_failure <- function(expr, message) {
  error <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(error, "error"), grepl(message, conditionMessage(error), fixed = TRUE))
}
expect_failure(find_ggchord_loci(sequences, features, "absent"), "not found")
expect_failure(find_ggchord_loci(sequences, features, "unknown", 10, 20), "equal upstream")
ambiguous <- features
ambiguous$feature_id[2] <- "forward"
expect_failure(find_ggchord_loci(sequences, ambiguous, "forward"), "ambiguous")
expect_failure(find_ggchord_loci(sequences, features, "forward", upstream = 0.5), "integer")
invalid <- features
invalid$end[1] <- 1001
expect_failure(find_ggchord_loci(sequences, invalid, "forward"), "within")
cat("Target locus development checks passed.\n")
