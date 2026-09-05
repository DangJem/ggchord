# Open a private text-measurement device and return its cleanup function.
# dev.off() selects the next open device, not necessarily the device that was
# active before this one. Nested measurement inside a preview must restore the
# preview explicitly, otherwise layout continues on the IDE's smaller device.
ggchord_measurement_device <- function(...) {
  previous <- grDevices::dev.cur()
  grDevices::pdf(NULL, ...)
  owned <- grDevices::dev.cur()
  function() {
    if (owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
    invisible(NULL)
  }
}
