# Highlight and replay the same R source that the reference offers for download.
args <- commandArgs(trailingOnly = TRUE)
inputs <- jsonlite::fromJSON(args[[1]], simplifyVector = FALSE)
output_dir <- args[[2]]
plot_width <- 1000L
plot_height <- 625L
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(".", export_all = FALSE, quiet = TRUE)

highlight <- function(source) {
  markup <- downlit::highlight(source, classes = downlit::classes_pandoc())
  if (is.na(markup)) stop("Unable to highlight reference code")
  # Downlit adds documentation links by default. Keep its token spans and the
  # terminal newline, which is omitted by its line splitter.
  markup <- gsub("<a\\b[^>]*>|</a>", "", markup, perl = TRUE)
  if (endsWith(source, "\n")) markup <- paste0(markup, "\n")
  doc <- xml2::read_html(paste0("<pre><code>", markup, "</code></pre>"), encoding = "UTF-8")
  if (!identical(xml2::xml_text(xml2::xml_find_first(doc, "//code")), source))
    stop("Syntax highlighting changed the reference code")
  markup
}

render <- function(topic, source) {
  env <- new.env(parent = globalenv())
  device <- NULL
  open_png <- function(...) {
    ragg::agg_png(file.path(output_dir, paste0(topic, "-%02d.png")),
      width = plot_width, height = plot_height, res = 120, background = "white")
    device <<- grDevices::dev.cur()
  }
  old_options <- options(device = open_png, width = 98, max.print = 240)
  on.exit({
    if (!is.null(device) && device %in% grDevices::dev.list())
      grDevices::dev.off(device)
    options(old_options)
  })
  capture.output({
    for (expr in parse(text = source, keep.source = FALSE)) {
      result <- withVisible(eval(expr, env))
      if (result$visible) print(result$value)
    }
  })
}

rendered_figures <- function(topic) {
  paths <- sort(Sys.glob(file.path(output_dir, paste0(topic, "-[0-9][0-9].png"))))
  # Device-based measurements can create a blank canvas before any plot is drawn.
  drawn <- Filter(function(path) {
    pixels <- png::readPNG(path, native = TRUE)
    any(pixels != pixels[1L])
  }, paths)
  basename(drawn)
}

results <- lapply(names(inputs$topics), function(topic) {
  entry <- inputs$topics[[topic]]
  tryCatch({
    signature <- highlight(entry$signature)
    example <- highlight(entry$example)
    console <- render(topic, entry$example)
    figures <- rendered_figures(topic)
    # Discard leading and trailing empty lines from the printed R transcript.
    console <- trimws(paste(console, collapse = "\n"), whitespace = "[\r\n]")
    if (!length(figures) && !nzchar(console)) stop("Example has no visible result")
    list(signature = signature, example = example, console = console, figures = as.list(figures))
  }, error = function(error) stop(topic, ": ", conditionMessage(error), call. = FALSE))
})
names(results) <- names(inputs$topics)
jsonlite::write_json(list(topics = results,
  width = plot_width, height = plot_height), file.path(output_dir, "results.json"), auto_unbox = TRUE)
