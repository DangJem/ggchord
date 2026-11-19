# Read the exact clipboard expressions in a clean working directory. A file
# device exercises the normal print method without opening an interactive window.
options(device = function(...) grDevices::png(filename = "plot-%03d.png", width = 900, height = 700))
script <- readLines(file("stdin"), warn = FALSE)
plot_seen <- FALSE
for (expr in parse(text = script)) {
  result <- withVisible(eval(expr, envir = .GlobalEnv))
  if (result$visible && inherits(result$value, c("ggplot", "ggchord_fusion", "grob", "gtable"))) {
    print(result$value)
    plot_seen <- TRUE
  }
  if (is.character(result$value) && length(result$value) == 1L &&
      !is.na(result$value) && grepl("[.](png|pdf|svg)$", result$value, ignore.case = TRUE) && file.exists(result$value))
    plot_seen <- TRUE
}
grDevices::graphics.off()
if (commandArgs(TRUE)[1] == "plot" && !plot_seen && !length(list.files(pattern = "[.](png|pdf|svg)$", ignore.case = TRUE)))
  stop("The copied plotting example returned no visible plot or saved figure")
