# Explicit finite RC benchmark, not a release or stress run.
# Rscript tools/benchmark-rc.R [1..5 repetitions]
# Reports observations, never a fixed-seconds machine-independent gate.
pkgload::load_all(quiet = TRUE)
source("tools/rc-fixtures.R")
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args)) suppressWarnings(as.integer(args[1L])) else 1L
stopifnot(length(args) <= 1L, !is.na(n), n >= 1L, n <= 5L)
local({
  grDevices::pdf(NULL, width = 8, height = 8, family = "Helvetica")
  on.exit(grDevices::dev.off(), add = TRUE)
  f <- rc_fixtures()
  factories <- list(
    chord = function() ggchord(f$seqs, f$links, validate = "none") + geom_seq() + geom_link_ribbon(),
    circular = function() f$circular(),
    collinear = function() f$linear,
    track = function() f$track)
  results <- list()
  for (case in names(factories)) for (i in seq_len(n)) {
    p <- factories[[case]]()
    # A fresh combination owns a fresh layout cache, including for stored fixtures.
    p <- p + theme(legend.position = "bottom")
    elapsed <- system.time(value <- ggplot2::ggplot_build(p))[["elapsed"]]
    results[[length(results) + 1L]] <- data.frame(case = case, repetition = i,
      elapsed_seconds = elapsed, built_object_bytes = as.numeric(object.size(value)),
      measure = "single build, fresh plot cache; object footprint is not peak allocation",
      sequence_rows = switch(case, chord = nrow(f$seqs), collinear = nrow(f$view$seq_data), 1L),
      primary_rows = switch(case, chord = nrow(f$links), circular = nrow(f$features) + nrow(f$primers) + nrow(f$sites),
        collinear = nrow(f$view$gene_data) + nrow(f$view$ribbon_data), track = nrow(f$bins)),
      R = R.version.string, ggplot2 = as.character(utils::packageVersion("ggplot2")),
      platform = R.version$platform, device = "PDF metrics 8x8 in; Helvetica; no output file",
      revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
    cat(case, i, sprintf("%.3fs", elapsed), "\n")
  }
  write.csv(do.call(rbind, results),
    "design-roadmap/audit-generated/rc-benchmark-v0.16.6.csv", row.names = FALSE)
})
