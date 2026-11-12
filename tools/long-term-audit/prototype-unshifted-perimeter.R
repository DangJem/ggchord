# Isolated packing experiment; never edits runtime sources or source records.
script <- readLines('tools/long-term-audit/check-flagship-label-quality.R', warn = FALSE)
cut <- grep('^results <- list', script)[1L]
eval(parse(text = script[seq_len(cut)]))
namespace <- asNamespace('ggchord')
candidate <- get('ggchord_place_dense_perimeter', namespace)
replacements <- 0L
events <- new.env(parent = emptyenv()); events$rows <- list()
rewrite <- function(expr) {
  if (!is.call(expr)) return(expr)
  if (identical(expr[[1L]], as.name('if')) &&
      grepl('^ggchord_perimeter_rail_crossings', paste(deparse(expr[[2L]]), collapse = ' '))) {
    replacements <<- replacements+1L
    return(quote(NULL))
  }
  as.call(lapply(as.list(expr), rewrite))
}
body(candidate) <- rewrite(body(candidate)); stopifnot(replacements == 1L)
private <- new.env(parent = namespace); private$events <- events
environment(candidate) <- private
unlockBinding('ggchord_place_dense_perimeter', namespace)
assign('ggchord_place_dense_perimeter', candidate, namespace)
lockBinding('ggchord_place_dense_perimeter', namespace)
prefix <- 'v0.27.0-unshifted-perimeter-prototype-'
eval(parse(text = script[-seq_len(cut)]))
file <- paste0('tools/long-term-audit/', prefix, paste(names(profiles), collapse = '-'), '-results.json')
report <- jsonlite::fromJSON(file, simplifyVector = FALSE)
report$baseline_runtime_sha256 <- report$runtime_sha256; report$runtime_sha256 <- NULL
report$production_override <- 'Keep preferred compact columns; defer crossings to the measured route solver; no production certification'
report$candidate_body_sha256 <- digest::digest(deparse(body(candidate)), algo = 'sha256')
report$candidate_events <- events$rows
jsonlite::write_json(report, file, pretty = TRUE, auto_unbox = TRUE, na = 'null')
