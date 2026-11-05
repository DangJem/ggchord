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
  grepl("\\{\\{(code|prepare|recipe):", text)
}, pages)
recipes <- jsonlite::fromJSON("website/showcase-recipes.json")
outputs <- list()
for(page in pages) {
  text <- paste(readLines(paste0('website/en/', page, '.md'), warn=FALSE), collapse='\n')
  for (name in names(recipes)) {
    text <- sub(paste0("{{recipe:", name, "}}"),
      paste0("{{prepare:", paste(recipes[[name]], collapse = ","), "}}"),
      text, fixed = TRUE)
  }
  matches <- regmatches(text, gregexpr('\\{\\{(code|prepare|output|console):[a-z,]+\\}\\}', text))[[1]]
  environment <- new.env(parent=globalenv())
  evaluated <- list()
  for(marker in matches) {
    if (startsWith(marker, "{{output:") || startsWith(marker, "{{console:")) {
      key <- sub("^.*:([a-z]+)[}][}]$", "\\1", marker)
      if (!key %in% names(evaluated)) {
        stop(page, ' / ', key, ': output requires its code first', call. = FALSE)
      }
      value <- evaluated[[key]]
      if (startsWith(marker, "{{console:")) {
        output <- list(console = paste(if (is.data.frame(value)) {
          capture.output(print(utils::head(value, 4L), width = 10000))
        } else {
          capture.output(str(value, max.level = 1L, vec.len = 3L, list.len = 6L,
                             give.attr = FALSE, width = 76))
        }, collapse = "\n"))
      } else output <- list(columns = as.list(names(value)), rows = lapply(seq_len(nrow(value)), function(i)
        as.list(unname(unlist(value[i, , drop = FALSE])))))
      if (!is.null(outputs[[key]]) && !identical(outputs[[key]], output)) {
        stop(page, ' / ', key, ': output differs between lessons', call. = FALSE)
      }
      outputs[[key]] <- output
      next
    }
    keys <- strsplit(sub('\\}\\}$','',sub('^\\{\\{(code|prepare):','',marker)),',',fixed=TRUE)[[1]]
    for(key in setdiff(keys,c('save','preview','frontsave','frontpreview','firstsave','firstsavebasic'))) {
      evaluated[key] <- list(tryCatch(eval(parse(text=chunks[[key]]),environment), error=function(e)
        stop(page, ' / ', key, ': ', conditionMessage(e), call.=FALSE)))
    }
  }
  cat('PASS:', page, 'in a separate documentation environment\n')
}
output_file <- 'website/assets/lesson-outputs.json'
serialized <- jsonlite::toJSON(outputs, auto_unbox = TRUE, pretty = TRUE)
if ('--write-outputs' %in% commandArgs(trailingOnly = TRUE)) {
  writeLines(serialized, output_file)
} else if (!identical(paste(readLines(output_file, warn = FALSE), collapse = '\n'), as.character(serialized))) {
  stop('Lesson outputs are stale. Run tools/check-documentation-lessons.R --write-outputs.')
}
cat(length(pages),'page preparations completed; no package checks or rendered exports.\n')
