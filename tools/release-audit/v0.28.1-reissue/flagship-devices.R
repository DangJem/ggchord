pkgload::load_all(quiet = TRUE)
source("tools/site-flagships.R")
out <- "/tmp/ggchord-v0281-republish-flagships"
dir.create(out, showWarnings = FALSE)
for (key in c("showcaseplasmiddata", "showcaseplasmid")) invisible(eval(parse(text = shared_blocks[[key]])))
metadata <- jsonlite::fromJSON("website/manual-figures.json")
cases <- metadata[metadata$id %in% c("flagchord", "flaglinear", "flagradial", "flagsankey", "flagdotplot", "flagfusion", "flagcircular"), ]
cases$name <- cases$object
stopifnot(nrow(cases) == 7L)
records <- list()
for (i in seq_len(nrow(cases))) for (format in c("png", "pdf", "svg")) {
  path <- file.path(out, paste0(cases$name[i], ".", format))
  warnings <- character()
  pdf <- if (cases$name[i] %in% c("flag_dotplot", "p_show_plasmid")) {
    function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height, ...)
  } else grDevices::pdf
  withCallingHandlers(ggplot2::ggsave(path, get(cases$name[i]), width = cases$width[i], height = cases$height[i], dpi = cases$dpi[i], bg = "#FFFCF5", device = switch(format, png = ragg::agg_png, pdf = pdf, svg = svglite::svglite)), warning = function(w) warnings <<- c(warnings, conditionMessage(w)))
  if (length(warnings)) stop(cases$name[i], "/", format, ": ", paste(warnings, collapse = "; "))
  stopifnot(file.exists(path), file.info(path)$size > 0)
  if (format == "pdf") stopifnot(qpdf::pdf_length(path) == 1L)
  if (format == "svg") stopifnot(xml2::xml_name(xml2::read_xml(path)) == "svg")
  records[[length(records) + 1L]] <- data.frame(case = cases$name[i], format, width = cases$width[i], height = cases$height[i], bytes = file.info(path)$size, sha256 = digest::digest(file = path, algo = "sha256"), result = "PASS")
  cat("PASS:", cases$name[i], format, "\n"); flush.console()
}
write.csv(do.call(rbind, records), file.path(out, "devices.csv"), row.names = FALSE)
