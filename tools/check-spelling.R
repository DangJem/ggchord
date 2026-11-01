# Explicit documentation-only spelling review. No tests, knitting or package build.
# Rscript tools/check-spelling.R
# Requires the optional development package spelling; never edits WORDLIST.
main <- function() {
  if (!requireNamespace("spelling", quietly = TRUE)) {
    stop("Install the optional development package 'spelling' to review documentation")
  }
  words <- readLines("inst/WORDLIST", encoding = "UTF-8", warn = FALSE)
  stopifnot(length(words) > 0L, !anyDuplicated(words),
    all(nzchar(words)), !any(grepl("[[:space:]]", words)))
  # Current DESCRIPTION has no Language field: spelling uses en-US for Rd and
  # metadata. The wordlist deliberately retains legitimate British spellings.
  rd <- spelling::spell_check_package(".", vignettes = FALSE, use_wordlist = TRUE)
  prose <- spelling::spell_check_files(
    c("README.md", "MIGRATION.md", "vignettes/ggchord_vignette.Rmd"),
    ignore = words, lang = "en_US")
  # An English dictionary cannot judge the Han navigation labels that appear
  # in paired documents; keep every finding containing an English letter.
  rd <- rd[grepl("[A-Za-z]", rd$word), , drop = FALSE]
  prose <- prose[grepl("[A-Za-z]", prose$word), , drop = FALSE]
  cat("Rd/DESCRIPTION remaining candidates (manual review, not confirmed errors):\n")
  print(rd)
  cat("English README/migration/tutorial remaining candidates:\n")
  print(prose)
  cat("WORDLIST:", length(words), "reviewed entries; no files modified.\n")
  invisible(list(rd = rd, english = prose))
}
main()
