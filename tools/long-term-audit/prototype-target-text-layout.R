# Isolated text-metric candidate; keep production source unchanged.
pkgload::load_all(quiet = TRUE)
native_candidate <- ggchord_text_boxes
replacements <- 0L
rewrite <- function(expr) {
  if (!is.call(expr)) return(expr)
  if (identical(expr, quote(ggchord_measurement_device(text = texts[valid])))) {
    replacements <<- replacements+1L
    return(quote(if (grDevices::dev.cur() > 1L) function() invisible(NULL) else
      ggchord_measurement_device(text = texts[valid])))
  }
  as.call(lapply(as.list(expr), rewrite))
}
body(native_candidate) <- rewrite(body(native_candidate))
stopifnot(replacements == 1L)
namespace <- asNamespace('ggchord')
unlockBinding('ggchord_text_boxes', namespace)
assign('ggchord_text_boxes', native_candidate, envir = namespace)
lockBinding('ggchord_text_boxes', namespace)
script <- readLines('tools/long-term-audit/check-label-size-quality.R', warn = FALSE)
script <- script[-grep('^pkgload::load_all', script)]
at <- grep('^cases <-', script)[1L]
script <- append(script, "prefix <- 'v0.27.0-label-size-target-prototype-'", after = at-1L)
eval(parse(text = script))
file <- 'tools/long-term-audit/v0.27.0-label-size-target-prototype-results.json'
report <- jsonlite::fromJSON(file, simplifyVector = FALSE)
report$baseline_runtime_sha256 <- report$runtime_sha256
report$runtime_sha256 <- NULL
report$production_override <- 'Only ggchord_text_boxes uses the target device when it is active'
report$candidate_body_sha256 <- digest::digest(deparse(body(native_candidate)), algo = 'sha256')
report$scope <- paste('In-process text-metric candidate only.', report$scope)
jsonlite::write_json(report, file, pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Target text layout prototype recorded; no production-source certification\n')
