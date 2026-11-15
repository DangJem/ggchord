# Website-only dependencies. Keep cache keys tied to this explicit list.
cran <- c("ggplot2", "rlang", "scales", "pkgload", "rmarkdown", "jsonlite",
  "pkgdown", "dplyr")
bioconductor <- c("Biostrings", "VariantAnnotation")
missing <- function(packages) packages[!vapply(packages, requireNamespace,
  logical(1), quietly = TRUE)]
needed <- missing(cran)
if (length(needed)) install.packages(needed, Ncpus = 2L)
needed <- missing(bioconductor)
if (length(needed)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(needed, ask = FALSE, update = FALSE, Ncpus = 2L)
}
needed <- missing(c(cran, bioconductor))
if (length(needed)) stop("Missing website dependencies: ", paste(needed, collapse = ", "))
cat("Website dependencies available; installed packages reused.\n")
