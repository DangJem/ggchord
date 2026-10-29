# Small reusable inputs keep default package tests independent of teaching assets.
smoke_sequences <- function() data.frame(accver = c("A", "B"), length = 1000)
smoke_ribbon <- function() data.frame(qaccver = "A", saccver = "B",
  qstart = 100, qend = 200, sstart = 400, send = 300)
smoke_genes <- function() data.frame(accver = c("A", "B"),
  start = c(100, 400), end = c(200, 500), strand = c("+", "-"),
  anno = c("alpha", "beta"))

build_ggchord_smoke <- function(plot) {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  ggplot2::ggplot_build(plot)
}
