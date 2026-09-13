# Lightweight export-size preview for ggchord/ggplot objects.

ggchord_preview_cleanup <- function(path, keep = 20L) {
  html <- list.files(path, pattern = "^ggchord-preview-.*\\.html$",
                     full.names = TRUE)
  if (length(html) <= keep) return(invisible(NULL))
  info <- file.info(html)
  html <- html[order(info$mtime, decreasing = TRUE)]
  old_html <- html[seq.int(keep + 1L, length(html))]
  old_stems <- sub("\\.html$", "", old_html)
  old_files <- c(old_html, paste0(old_stems, ".png"), paste0(old_stems, ".svg"))
  unlink(old_files[file.exists(old_files)])
  invisible(NULL)
}

ggchord_preview_pixels <- function(value, units, dpi) {
  switch(units,
    "in" = value * dpi,
    cm = value / 2.54 * dpi,
    mm = value / 25.4 * dpi,
    px = value
  )
}

# Measure the fixed gtable decorations and its equal-unit panel. The null
# panel dimensions carry its aspect ratio; guide boxes contribute their real
# physical widths/heights. This fits whitespace without stretching geometry.
ggchord_preview_layout <- function(plot, width_inches, height_inches = 8) {
  close_device <- ggchord_measurement_device(width = width_inches, height = height_inches)
  on.exit(close_device())
  table <- ggplot2::ggplotGrob(plot)
  fixed_width <- grid::convertWidth(sum(table$widths), "inches", valueOnly = TRUE)
  fixed_height <- grid::convertHeight(sum(table$heights), "inches", valueOnly = TRUE)
  horizontal_guides <- which(table$layout$name %in% c("guide-box-top", "guide-box-bottom"))
  guide_width <- max(c(0, vapply(horizontal_guides, function(i) {
    grid::convertWidth(grid::grobWidth(table$grobs[[i]]), "inches", valueOnly = TRUE)
  }, numeric(1))))
  if (fixed_width >= width_inches || guide_width > width_inches) {
    warning("view_ggchord(): legends exceed the requested width; increase width ",
      "or set nrow in guide_ggchord_legend() via guides()",
      call. = FALSE)
  }
  nw <- grid::unitType(table$widths) == "null"
  nh <- grid::unitType(table$heights) == "null"
  if (!any(nw) || !any(nh)) return(list(plot = table, height = width_inches * 0.7))
  ratio <- sum(as.numeric(table$heights[nh])) / sum(as.numeric(table$widths[nw]))
  panel_height <- max(width_inches - fixed_width, 0.5) * ratio
  guide_rows <- which(table$layout$name %in% c("guide-box-left", "guide-box-right"))
  guide_height <- max(c(0, vapply(guide_rows, function(i) {
    grid::convertHeight(grid::grobHeight(table$grobs[[i]]), "inches", valueOnly = TRUE)
  }, numeric(1))))
  fitted_height <- max(1, fixed_height + max(panel_height, guide_height))
  # Freeze the measured gtable for export: rebuilding at a new height can
  # select different label tracks and oscillate between two aspect ratios.
  list(plot = table, height = fitted_height)
}

ggchord_preview_height <- function(plot, width_inches) {
  ggchord_preview_layout(plot, width_inches)$height
}

#' Preview a ggchord plot at its intended export size
#'
#' Renders a plot with \code{ggsave()} into a temporary PNG or SVG and opens a
#' small static preview page in the configured IDE viewer or web browser. This
#' is an export-size preview, not an interactive plot conversion; the original
#' plot and its standard \code{ggsave()} workflow remain unchanged.
#' Automatic height fitting does not resize an overwide legend. If decorations
#' exceed the requested width, increase \code{width} or use, for example,
#' \code{guides(seq_colour = guide_ggchord_legend(nrow = 2))}.
#'
#' @section RStudio plot pane and export-size rendering:
#' Label packing, guide sizing, and fitted limits use the active graphics
#' device's physical dimensions. Printing \code{p} in a small or narrow
#' RStudio Plots pane can therefore produce a substantially different layout
#' from an export, and very small panes may not have enough room for all radial
#' labels. Resize or Zoom the pane and print again for a larger direct preview.
#' For size-sensitive chord plots, use \code{view_ggchord()} while composing,
#' then save the final figure with \code{ggsave()} at an explicit width and
#' height matching the intended output. Preview measurements restore the
#' previously active graphics device.
#'
#' @param plot A ggchord or ggplot object, default \code{last_plot()}.
#' @param width Positive output width, default 12.39 inches when \code{units = "in"}.
#' @param height Positive output height, default 9.71 inches. Use \code{NULL} to fit the
#'   sequence, labels and legends while preserving equal coordinate units.
#'   Explicit width/height values are always respected.
#' @param units Output units: \code{"in"}, \code{"cm"}, \code{"mm"}, or
#'   \code{"px"}.
#' @param device Preview device, \code{"png"} or \code{"svg"}. SVG output
#'   requires the optional svglite package or a working base Cairo SVG device.
#' @param dpi Positive raster resolution, default 144.
#' @param bg Optional background colour passed to \code{ggsave()}.
#' @param viewer How to open the preview: \code{"auto"} prefers
#'   \code{getOption("viewer")} and falls back to a browser in interactive
#'   sessions; \code{"ide"}, \code{"browser"}, and \code{"none"} select an
#'   explicit destination. \code{"none"} is useful in scripts and tests.
#'
#' @return The normalized temporary image path, invisibly.
#' @export
#' @examples
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq()
#' if (interactive()) view_ggchord(p)
view_ggchord <- function(
    plot = ggplot2::last_plot(),
    width = 12.39,
    height = 9.71,
    units = c("in", "cm", "mm", "px"),
    device = c("png", "svg"),
    dpi = 144,
    bg = NULL,
    viewer = c("auto", "ide", "browser", "none")) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  caller <- "view_ggchord()"
  if (!inherits(plot, "ggplot")) {
    ggchord_stop(caller, ": plot must be a ggchord or ggplot object")
  }
  units <- match.arg(units)
  device <- match.arg(device)
  viewer <- match.arg(viewer)
  valid_dimension <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  }
  if (!valid_dimension(width) || (!is.null(height) && !valid_dimension(height))) {
    ggchord_stop(caller, ": width and height must be positive finite numbers (or height = NULL)")
  }
  if (!valid_dimension(dpi)) {
    ggchord_stop(caller, ": dpi must be one positive finite number")
  }
  if (!is.null(bg) &&
      (!is.character(bg) || length(bg) != 1L || is.na(bg))) {
    ggchord_stop(caller, ": bg must be NULL or one colour string")
  }

  preview_plot <- plot
  if (is.null(height)) {
    width_inches <- ggchord_preview_pixels(width, units, dpi) / dpi
    fitted <- ggchord_preview_layout(plot, width_inches)
    height_inches <- fitted$height
    preview_plot <- fitted$plot
    height <- switch(units, "in" = height_inches, cm = height_inches * 2.54,
                     mm = height_inches * 25.4, px = height_inches * dpi)
  }

  preview_dir <- file.path(tempdir(), "ggchord-preview")
  if (!dir.exists(preview_dir)) {
    dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
  }
  stem <- tempfile("ggchord-preview-", tmpdir = preview_dir)
  image_file <- paste0(stem, ".", device)
  html_file <- paste0(stem, ".html")

  save_device <- if (identical(device, "svg") &&
      requireNamespace("svglite", quietly = TRUE)) {
    "svg"
  } else if (identical(device, "svg")) {
    grDevices::svg
  } else {
    "png"
  }
  tryCatch(
    ggplot2::ggsave(
        filename = image_file,
        plot = preview_plot,
        device = save_device,
        width = width,
        height = height,
        units = units,
        dpi = dpi,
        bg = bg
      ),
    error = function(e) {
      if (identical(device, "svg")) {
        ggchord_stop(
          caller, ": SVG preview requires the optional svglite package or ",
          "a working Cairo SVG device; ", conditionMessage(e)
        )
      }
      stop(e)
    }
  )
  if (!file.exists(image_file)) {
    if (identical(device, "svg")) {
      ggchord_stop(
        caller, ": SVG preview requires the optional svglite package or ",
        "a working Cairo SVG device"
      )
    }
    ggchord_stop(caller, ": preview device did not create an output file")
  }

  width_px <- max(1L, round(ggchord_preview_pixels(width, units, dpi)))
  height_px <- max(1L, round(ggchord_preview_pixels(height, units, dpi)))
  image_name <- utils::URLencode(basename(image_file), reserved = TRUE)
  metadata <- paste0(
    format(width, trim = TRUE), " ", units, " x ",
    format(height, trim = TRUE), " ", units,
    "  |  ", toupper(device),
    if (identical(device, "png")) paste0("  |  ", dpi, " dpi") else ""
  )
  html <- c(
    "<!doctype html>",
    "<html><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>ggchord preview</title>",
    "<style>",
    "html,body{margin:0;min-height:100%;background:#eef1f4;color:#30343b;",
    "font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif}",
    ".bar{position:sticky;top:0;z-index:1;padding:8px 12px;background:#fff;",
    "border-bottom:1px solid #d8dde3;font-size:12px}",
    ".stage{padding:18px;overflow:auto;text-align:center}",
    ".frame{display:inline-block;background:#fff;box-shadow:0 2px 14px #0002;",
    "line-height:0}",
    paste0("img{display:block;width:", width_px, "px;height:", height_px,
           "px;object-fit:contain}"),
    "</style></head><body>",
    paste0("<div class=\"bar\">ggchord preview &nbsp; ", metadata, "</div>"),
    paste0("<div class=\"stage\"><div class=\"frame\"><img src=\"",
           image_name, "\" alt=\"ggchord preview\"></div></div>"),
    "</body></html>"
  )
  writeLines(html, html_file, useBytes = TRUE)
  ggchord_preview_cleanup(preview_dir, keep = 20L)

  html_file <- normalizePath(html_file, winslash = "/", mustWork = TRUE)
  image_file <- normalizePath(image_file, winslash = "/", mustWork = TRUE)
  ide_viewer <- getOption("viewer")
  if (identical(viewer, "ide")) {
    if (!is.function(ide_viewer)) {
      ggchord_stop(caller, ": no IDE viewer is configured in getOption('viewer')")
    }
    ide_viewer(html_file)
  } else if (identical(viewer, "browser")) {
    utils::browseURL(html_file)
  } else if (identical(viewer, "auto")) {
    if (is.function(ide_viewer)) {
      ide_viewer(html_file)
    } else if (interactive()) {
      utils::browseURL(html_file)
    }
  }

  invisible(image_file)
}
