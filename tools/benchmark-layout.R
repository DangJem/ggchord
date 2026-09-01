# Lightweight manual benchmark for the v0.10.0 layout pipeline.
# Run from the package root with: Rscript tools/benchmark-layout.R

devtools::load_all(quiet = TRUE)

make_ribbons <- function(n, seq_length = 100000) {
  starts <- seq(1, seq_length - 200, length.out = n)
  data.frame(
    qaccver = rep("A", n),
    saccver = rep("B", n),
    length = rep(200, n),
    pident = seq(50, 100, length.out = n),
    qstart = starts,
    qend = starts + 199,
    sstart = rev(starts),
    send = rev(starts) + 199
  )
}

seq_data <- data.frame(
  seq_id = c("A", "B"), length = c(100000, 100000)
)

for (n in c(1000L, 5000L)) {
  ribbons <- make_ribbons(n)
  plot <- ggchord(seq_data, ribbons, validate = "none") +
    geom_seq() + geom_ribbon()
  elapsed <- system.time(invisible(ggplot2::ggplot_build(plot)))[["elapsed"]]
  message(sprintf("%d ribbons: %.3f seconds", n, elapsed))
}
