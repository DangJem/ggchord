#' Preview a plot using explicit device dimensions
#' @inheritParams export_ggchord_layout
#' @return Invisibly, the plot, measurement metadata and diagnostics.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' view_ggchord(p_basic, width = 180, height = 100, units = "mm")
view_ggchord <- function(plot, width = 8, height = 6, units = "in",
    device = "png", dpi = 150) {
  caller <- "view_ggchord()"
  ggchord_output_plot(plot, caller)
  request <- ggchord_output_request(width, height, units, device, dpi, caller)
  if (is.null(request$size) || is.null(request$device))
    ggchord_stop(caller, ": specify a width, height and preview device")
  shown <- interactive()
  warnings <- character()
  preview <- withCallingHandlers(ggchord_with_output_device(request, function() {
    if (shown) grDevices::dev.control(displaylist = "enable")
    result <- ggchord_output_build_geometry(ggchord_query_clone(plot))
    list(issues = ggchord_diagnose_built(result, request$size * 25.4),
      recording = if (shown) grDevices::recordPlot())
  }), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if (!is.null(preview$recording)) grDevices::replayPlot(preview$recording)
  measurement <- ggchord_output_measurement(request, TRUE,
    actual = attr(preview, "measurement_device"))
  warning_rows <- lapply(unique(warnings), function(message)
    ggchord_output_issue("warning", "plot_warning", message))
  issues <- ggchord_rbind_fill(c(list(preview$issues), warning_rows))
  invisible(list(plot = plot, measurement = measurement,
    diagnostics = structure(list(issues = issues, measurement = measurement),
      class = "ggchord_diagnostics")))
}
