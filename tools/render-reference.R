# Render help text only; never execute examples or build vignettes.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the output directory")
out <- args[[1L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(quiet = TRUE)
Sys.setenv("_R_HELP_LINKS_TO_TOPICS_" = "false")
files <- sort(list.files("man", "\\.Rd$", full.names = TRUE))
rds <- lapply(files, tools::parse_Rd)
tag_text <- function(rd, tag) paste(vapply(Filter(function(x)
  identical(attr(x, "Rd_tag"), tag), rd), function(x) paste(unlist(x), collapse = ""), ""), collapse = " ")
aliases <- lapply(rds, function(rd) vapply(Filter(function(x)
  identical(attr(x, "Rd_tag"), "\\alias"), rd), function(x) paste(unlist(x), collapse = ""), ""))
stems <- sub("\\.Rd$", "", basename(files))
links <- setNames(rep(paste0(stems, ".html"), lengths(aliases)), unlist(aliases))
metadata <- vector("list", length(files))
for (i in seq_along(files)) {
  tools::Rd2HTML(rds[[i]], out = file.path(out, paste0(stems[i], ".html")),
    package = "ggchord", Links = links, Links2 = links, texmath = "ascii")
  metadata[[i]] <- list(slug = stems[i], title = tag_text(rds[[i]], "\\title"),
    aliases = unname(aliases[[i]]), source = files[i])
}
writeLines(jsonlite::toJSON(metadata, auto_unbox = TRUE, pretty = TRUE), file.path(out, "index.json"))
