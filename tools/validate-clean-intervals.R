# Explicit development check for inclusive cleaning, without package tests.
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver = c("q", "s"), length = c(10, 20))
cases <- c("single", "before", "after", "clip_left", "clip_right")
genes <- data.frame(accver = "q", start = c(5, -3, 11, -4, 10),
  end = c(5, 0, 12, 1, 13), strand = "+", case = cases)
links <- data.frame(qaccver = "q", saccver = "s",
  qstart = c(5, -3, 2, -4, 10), qend = c(5, 0, 4, 1, 13),
  sstart = c(7, 2, 21, 6, 15), send = c(7, 5, 25, 2, 18),
  orientation_known = c(FALSE, TRUE, TRUE, TRUE, TRUE), case = cases)
out <- clean_ggchord_data(seqs, links, genes)
stopifnot(identical(out$gene_data$case, cases[c(1, 4, 5)]),
  identical(out$gene_data$start, c(5, 1, 10)),
  identical(out$gene_data$end, c(5, 1, 10)),
  identical(out$ribbon_data$case, cases[c(1, 4, 5)]),
  identical(out$ribbon_data$qstart, c(5, 1, 10)),
  identical(out$ribbon_data$qend, c(5, 1, 10)),
  !out$ribbon_data$orientation_known[1L])
removed <- out$report[out$report$action == "drop", ]
stopifnot(nrow(removed) == 4L,
  all(grepl("no overlap", removed$reason)),
  identical(sort(unique(removed$row)), c(2L, 3L)))
cat("Inclusive single-base and clipping checks passed.\n")
