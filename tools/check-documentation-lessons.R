# Verify that task pages prepare their own R objects, without building the package,
# knitting vignettes, saving user files, opening previews or running a test suite.
pkgload::load_all(quiet = TRUE)
lines <- readLines('tools/documentation-examples.R', warn = FALSE)
markers <- grep('^# ---- [a-z]+ ----$', lines)
chunks <- setNames(lapply(seq_len(length(markers)-1L), function(i)
  lines[seq.int(markers[i]+1L, markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$', '\\1', lines[markers[-length(markers)]]))
pages <- jsonlite::fromJSON("website/pages.json")$slug
pages <- Filter(function(page) {
  text <- paste(readLines(paste0("website/en/", page, ".md"), warn = FALSE), collapse = "\n")
  grepl("\\{\\{(code|prepare):", text)
}, pages)
for(page in pages) {
  text <- paste(readLines(paste0('website/en/', page, '.md'), warn=FALSE), collapse='\n')
  matches <- regmatches(text, gregexpr('\\{\\{(code|prepare):[a-z,]+\\}\\}', text))[[1]]
  environment <- new.env(parent=globalenv())
  for(marker in matches) {
    keys <- strsplit(sub('\\}\\}$','',sub('^\\{\\{(code|prepare):','',marker)),',',fixed=TRUE)[[1]]
    for(key in setdiff(keys,c('save','preview','firstsave'))) {
      tryCatch(eval(parse(text=chunks[[key]]),environment), error=function(e)
        stop(page, ' / ', key, ': ', conditionMessage(e), call.=FALSE))
    }
  }
  cat('PASS:', page, 'in a separate documentation environment\n')
}
cat(length(pages),'page preparations completed; no package checks or rendered exports.\n')
