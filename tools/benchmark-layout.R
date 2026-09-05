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
    geom_seq() + geom_link_ribbon()
  raw_elapsed <- system.time(
    invisible(ggplot2::ggplot_build(plot))
  )[["elapsed"]]

  bundle_elapsed <- system.time(
    bundled <- bundle_ggchord_ribbons(ribbons, seq_data)
  )[["elapsed"]]
  optimize_elapsed <- system.time(
    optimized <- optimize_ggchord_layout(seq_data, ribbons)
  )[["elapsed"]]

  bundled_plot <- ggchord(seq_data, bundled$data, validate = "none") +
    geom_seq(
      seq_order = optimized$seq_order,
      seq_orientation = optimized$seq_orientation
    ) +
    geom_link_ribbon(aes(ribbon_alpha = .bundle_density)) +
    scale_ribbon_alpha_continuous()
  bundled_build_elapsed <- system.time(
    invisible(ggplot2::ggplot_build(bundled_plot))
  )[["elapsed"]]

  message(sprintf(
    paste(
      "%d ribbons: raw build %.3fs; bundle %.3fs (%d rows);",
      "optimize %.3fs (%d scored, approximate=%s); bundled build %.3fs"
    ),
    n, raw_elapsed, bundle_elapsed, nrow(bundled$data),
    optimize_elapsed, optimized$report$n_scored,
    optimized$report$approximate, bundled_build_elapsed
  ))
}
