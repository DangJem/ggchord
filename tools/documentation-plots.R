# One PNG documentation path, not a release-device matrix.
# Run from the repository root. All consumers share these same files.
pkgload::load_all(quiet = TRUE)
text <- readLines("tools/documentation-examples.R", warn = FALSE)
markers <- grep("^# ---- [a-z]+ ----$", text)
chunks <- setNames(lapply(seq_len(length(markers) - 1L), function(i)
  text[seq.int(markers[i] + 1L, markers[i + 1L] - 1L)]),
  sub("^# ---- ([a-z]+) ----$", "\\1", text[markers[-length(markers)]]))
env <- new.env(parent = globalenv())
for (name in setdiff(names(chunks), c("save", "preview", "imports")))
  eval(parse(text = chunks[[name]]), envir = env)
stopifnot(nrow(env$genomic) > 0L, nrow(env$primers) > 0L, nrow(env$sites) > 0L,
  nrow(env$rows) > 0L)
assets <- data.frame(id = c("chord", "circular", "collinear", "tracks", "focus", "labels"),
  object = c("p_chord", "p_circular", "p_collinear", "p_tracks", "p_focus", "p_labels"),
  width = c(7, 7, 8, 7, 8, 8), height = c(6, 6, 5, 6, 5, 5),
  units = "in", dpi = 144, device = "png", background = "white",
  input = "tools/documentation-examples.R (synthetic teaching inputs)",
  generator = "tools/documentation-plots.R")
dir.create("man/figures", showWarnings = FALSE)
for (i in seq_len(nrow(assets))) {
  path <- paste0("man/figures/", assets$id[i], ".png")
  ggplot2::ggsave(path, get(assets$object[i], env), width = assets$width[i],
    height = assets$height[i], units = assets$units[i], dpi = assets$dpi[i],
    device = "png", bg = assets$background[i])
}
# Remaining documented geometry combinations are built once, without another device.
for (name in c("p_mapped", "p_sites", "p_segments", "p_migrated", "p_offsets"))
  invisible(ggplot2::ggplot_build(get(name, env)))
assets$file <- paste0("man/figures/", assets$id, ".png")
assets$md5 <- unname(tools::md5sum(assets$file))
assets$code_md5 <- unname(tools::md5sum("tools/documentation-examples.R"))
write.csv(assets, "website/assets/plots.csv", row.names = FALSE)
cat("Six shared PNG figures and documented geometry snippets completed.\n")
