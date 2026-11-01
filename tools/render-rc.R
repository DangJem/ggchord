# One explicit visual scene at a time, PNG only. Not a release device matrix.
# Rscript tools/render-rc.R [circular|inside|dense|chord|collinear|track] [output-dir]
# Images/session/manifest go only to .codex-tmp or an explicitly supplied directory.
pkgload::load_all(quiet = TRUE)
source("tools/rc-fixtures.R")
args <- commandArgs(trailingOnly = TRUE)
case <- if (length(args)) args[1L] else "circular"
stopifnot(case %in% c("circular", "inside", "dense", "chord", "collinear", "track"), length(args) <= 2L)
out <- if (length(args) == 2L) args[2L] else ".codex-tmp/rc-v0.16.6/visual"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
f <- rc_fixtures()
p <- switch(case,
  circular = f$circular(), inside = f$circular(side = "inside"),
  dense = f$circular(dense = TRUE), collinear = f$linear, track = f$track,
  chord = ggchord(f$seqs, f$links, validate = "none") + geom_seq() + geom_link_ribbon() + geom_gene(data = f$genes))
file <- file.path(out, paste0("rc-", case, "-v0.16.6.png"))
ggplot2::ggsave(file, p + theme(text = element_text(family = "sans")),
  width = 9, height = 8, units = "in", dpi = 120, device = "png", bg = "white")
writeLines(capture.output(sessionInfo()), file.path(out, paste0(case, "-session.txt")))
write.csv(data.frame(case = case, file = normalizePath(file), md5 = unname(tools::md5sum(file)),
  width = 9, height = 8, units = "in", dpi = 120, device = "png", family = "sans",
  source = "tools/rc-fixtures.R; bundled references, source_strand and exact CDS DNA",
  inspection = "RENDERED; human/model visual review recorded separately"),
  file.path(out, paste0(case, "-manifest.csv")), row.names = FALSE)
cat(normalizePath(file), "\n")
