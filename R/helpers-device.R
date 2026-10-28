# Open a private text-measurement device and return its cleanup function.
# dev.off() selects the next open device, not necessarily the device that was
# active before this one. Nested measurement inside a preview must restore the
# preview explicitly, otherwise layout continues on the IDE's smaller device.
ggchord_measurement_device <- function(text = NULL, ...) {
  previous <- grDevices::dev.cur()
  unicode <- is.null(text) || any(grepl("[^\\x01-\\x7F]", text,
    perl = TRUE), na.rm = TRUE)
  path <- NULL
  opened <- FALSE
  on.exit(if (!opened && !is.null(path)) unlink(path))
  if (!unicode) {
    # Keep the established metrics for ordinary ASCII-only layout requests.
    grDevices::pdf(NULL, ...)
  } else {
    # Base PDF uses single-byte font encodings on some platforms. Its text
    # bounds can fail before drawing a Unicode label such as lacZ-alpha.
    # These devices measure the original string without changing the caller's
    # chosen output device or transliterating its biological annotation.
    path <- tempfile("ggchord-measure-", fileext = ".pdf")
    if (isTRUE(capabilities("aqua"))) {
      grDevices::quartz(type = "pdf", file = path, ...)
    } else if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf(filename = path, ...)
    } else if (.Platform$OS.type == "windows") {
      path <- sub("[.]pdf$", ".png", path)
      args <- list(...)
      args$filename <- path
      args$width <- args$width %||% 7
      args$height <- args$height %||% 7
      args$units <- "in"
      args$res <- 144
      do.call(grDevices::png, args)
    } else if (requireNamespace("svglite", quietly = TRUE)) {
      path <- sub("[.]pdf$", ".svg", path)
      svglite::svglite(file = path, ...)
    } else if (!is.null(text)) ggchord_stop(
      "Unicode text measurement needs a Unicode-capable graphics device; ",
      "use an R build with Cairo support or install the optional svglite package"
    ) else {
      # Unknown preview contents do not by themselves require an optional
      # backend on a platform where ordinary ASCII PDF measurement works.
      path <- NULL
      grDevices::pdf(NULL, ...)
    }
  }
  owned <- grDevices::dev.cur()
  opened <- TRUE
  closed <- FALSE
  function() {
    if (closed) return(invisible(NULL))
    closed <<- TRUE
    on.exit({
      if (!is.null(path)) unlink(path)
      if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
    })
    if (owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
    invisible(NULL)
  }
}
