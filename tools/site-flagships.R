# Each download is the canonical, standalone source for one complete website figure.
# Keep its intermediate plots available to the homepage preview at a smaller size.
flagship_cases <- lapply(c("chord", "collinear", "radial", "sankey", "circular",
  "dotplot", "fusion"), function(view) {
  env <- new.env(parent = globalenv())
  env$complete <- source(file.path("website/examples", paste0("case-", view, ".R")),
    local = env)$value
  env
})
names(flagship_cases) <- c("chord", "collinear", "radial", "sankey", "circular", "dotplot", "fusion")
