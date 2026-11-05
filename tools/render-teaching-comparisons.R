# Helpers for the registered independent comparison family.
ggchord_comparison_assets <- function(registry) {
  do.call(rbind, lapply(registry, function(row) data.frame(
    id = row$id, object = row$object, width = row$width, height = row$height,
    units = "in", dpi = 144, device = "png", background = "white",
    generator = "tools/documentation-plots.R", input = row$input,
    coordinate = row$coordinate, stringsAsFactors = FALSE
  )))
}

ggchord_render_comparisons <- function(registry, requested, chunks, fingerprints) {
  previous <- read.csv("website/assets/plots.csv", stringsAsFactors = FALSE)
  retained <- previous[!previous$id %in% requested, , drop = FALSE]
  current <- unname(unlist(fingerprints$assets)[retained$id])
  if (anyNA(current) || any(current != retained$dependency_md5) ||
      any(unname(tools::md5sum(retained$file)) != retained$md5)) {
    stop("Unrequested figure dependencies or bytes changed; use the full generator")
  }
  selected <- Filter(function(row) row$id %in% requested, registry)
  environment <- new.env(parent = globalenv())
  evaluated <- character()
  assets <- ggchord_comparison_assets(selected)
  for (i in seq_along(selected)) {
    row <- selected[[i]]
    for (key in setdiff(unlist(row$snippets), evaluated)) {
      eval(parse(text = chunks[[key]]), environment)
      evaluated <- c(evaluated, key)
    }
    plot <- get(row$object, environment)
    actual <- if (inherits(plot$coordinates, "CoordCollinear")) "collinear" else
      if (inherits(plot$coordinates, "CoordCircular")) "circular" else "chord"
    stopifnot(identical(actual, row$coordinate))
    preview <- view_ggchord(plot, width = row$width, height = row$height,
      units = "in", dpi = 144, device = "png", bg = "white", viewer = "none")
    destination <- paste0("man/figures/", row$id, ".png")
    stopifnot(file.copy(preview, destination, overwrite = TRUE))
    cat("Rendered comparison:", row$id, "\n")
  }
  assets$file <- paste0("man/figures/", assets$id, ".png")
  assets$md5 <- unname(tools::md5sum(assets$file))
  assets$code_md5 <- fingerprints$code
  assets$dependency_md5 <- unname(unlist(fingerprints$assets)[assets$id])
  write.csv(rbind(retained, assets), "website/assets/plots.csv", row.names = FALSE)
}
