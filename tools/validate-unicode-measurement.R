# Direct checks for the Unicode measurement failure reported by Ubuntu CI.
# Run explicitly from the repository root; no package tests/checks are run.
# macOS may not reproduce the original base-PDF conversion error. This script
# checks the original Unicode data, private measurement, and device ownership;
# a Linux/Windows CI rerun is still needed to establish those platform results.
pkgload::load_all(quiet = TRUE)

validate_unicode_measurement <- function() {
  old_options <- options(warn = 2)
  on.exit(options(old_options), add = TRUE)
  original <- grDevices::dev.list()
  previous <- grDevices::dev.cur()
  on.exit({
    for (id in setdiff(grDevices::dev.list(), original)) grDevices::dev.off(id)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
  }, add = TRUE)

  annotation <- data.frame(text = "lacZ\u03b1", text_x = 0, text_y = 0, size = 2.5)
  boxes <- ggchord:::ggchord_text_boxes(annotation)
  stopifnot(all(is.finite(boxes$w)), all(boxes$w > 0), all(boxes$h > 0),
    identical(annotation$text, "lacZ\u03b1"),
    identical(grDevices::dev.cur(), previous))

  # Measuring must preserve an explicitly selected caller device, including
  # base PDF. It does not make that output device able to draw Unicode.
  grDevices::pdf(NULL, width = 6, height = 4)
  caller <- grDevices::dev.cur()
  nested <- function() {
    close_outer <- ggchord:::ggchord_measurement_device(text = "ASCII")
    on.exit(close_outer())
    outer <- grDevices::dev.cur()
    ggchord:::ggchord_text_boxes(annotation)
    stopifnot(identical(grDevices::dev.cur(), outer))
    stop("deliberate measurement-scope error")
  }
  error <- tryCatch(nested(), error = identity)
  stopifnot(inherits(error, "error"),
    grepl("deliberate", conditionMessage(error)),
    identical(grDevices::dev.cur(), caller))

  # A cleanup called twice must not close a new caller device reusing its ID.
  close_private <- ggchord:::ggchord_measurement_device(text = annotation$text)
  close_private()
  grDevices::pdf(NULL)
  replacement <- grDevices::dev.cur()
  close_private()
  stopifnot(identical(grDevices::dev.cur(), replacement))
  grDevices::dev.off(replacement)
  stopifnot(identical(grDevices::dev.cur(), caller))

  # The actual public CI input retains its biological label through layout.
  data(plasmid_example_pBluescript_II_SK_plus)
  features <- find_common_features(plasmid_example_pBluescript_II_SK_plus)
  tracks <- position_feature_stack(spacing = .10,
    base_position = position_plasmid())
  plot <- ggchord(plasmid_example_pBluescript_II_SK_plus) + geom_seq() +
    geom_feature_plasmid(data = features, position = tracks) +
    geom_feature_label_repel(data = features, position = tracks) +
    coord_circular(rotation = 90)
  layout <- export_ggchord_layout(plot, include = c("feature", "labels"))
  stopifnot(any(layout$labels$label == "lacZ\u03b1", na.rm = TRUE),
    all(is.finite(layout$feature$x)),
    identical(grDevices::dev.cur(), caller),
    !length(list.files(tempdir(), pattern = "^ggchord-measure-")))
  cat("Unicode annotation measurement and nested device ownership passed.\n")
}

validate_unicode_measurement()
