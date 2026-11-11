# Research helpers only; not exported or sourced by the package runtime.
# Character display maps individual source positions, not alignment columns.
source_letters <- function(seqs, accver, start, end, alphabet = c("DNA", "RNA", "AA"),
                           reverse = FALSE) {
  alphabet <- match.arg(alphabet)
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) stop("invalid reverse")
  if (!is.data.frame(seqs) || !all(c("accver", "length", "sequence") %in% names(seqs)))
    stop("actual sequence strings are required")
  at <- which(seqs$accver == accver)
  if (length(at) != 1L) stop("select exactly one source sequence")
  text <- seqs$sequence[at]
  if (is.na(text) || nchar(text) != seqs$length[at]) stop("source string/length mismatch")
  if (length(start) != 1L || length(end) != 1L ||
      any(!is.finite(c(start, end))) || any(c(start, end) != floor(c(start, end))) ||
      start < 1 || end < start || end > nchar(text)) stop("invalid source interval")
  letters <- strsplit(substr(text, start, end), "", fixed = TRUE)[[1]]
  allowed <- switch(alphabet, DNA = "ACGTRYSWKMBDHVN", RNA = "ACGURYSWKMBDHVN",
    AA = "ACDEFGHIKLMNPQRSTVWYBXZJUO*")
  if (any(!toupper(letters) %in% strsplit(allowed, "", fixed = TRUE)[[1]]))
    stop("alphabet mismatch or gap: supply an explicit alignment map for gapped input")
  position <- seq.int(start, end)
  original <- letters
  if (reverse) {
    position <- rev(position); original <- rev(original); letters <- rev(letters)
    if (alphabet != "AA") {
      to <- switch(alphabet, DNA = "TGCAYRSWMKVHDBN", RNA = "UGCAYRSWMKVHDBN")
      letters <- chartr(paste0(allowed, tolower(allowed)), paste0(to, tolower(to)), letters)
    }
  }
  origin <- if (".source_origin" %in% names(seqs)) seqs$.source_origin[at] else 1
  id <- if (".source_accver" %in% names(seqs)) seqs$.source_accver[at] else accver
  data.frame(accver = accver, source_accver = id, source_position = position + origin - 1,
    local_position = position, display_position = seq_along(position),
    source_symbol = original, symbol = letters, alphabet = alphabet,
    sequence_unit = if (alphabet == "AA") "aa" else "bp", reverse = reverse)
}

# Measure on the actual active device; size is physical mm, independent of DPI.
letter_budget <- function(letters, width_mm, size_mm = 3) {
  if (!nrow(letters) || length(width_mm) != 1L || !is.finite(width_mm) || width_mm <= 0)
    stop("positive width and nonempty letters required")
  metrics <- grid::textGrob(unique(letters$symbol),
    gp = grid::gpar(fontfamily = "mono", fontsize = size_mm * 72.27 / 25.4))
  glyph <- grid::convertWidth(grid::grobWidth(metrics), "mm", valueOnly = TRUE)
  required <- glyph * 1.25
  list(representation = if (width_mm / nrow(letters) >= required) "letters" else "tiles",
    available_mm_per_position = width_mm / nrow(letters), required_mm_per_position = required,
    size_mm = size_mm, n_positions = nrow(letters),
    reason = "Measured glyph width plus spacing; narrow panels retain tiles without shrinking fonts")
}

# M0 deliberately accepts only untransformed native scales and a single locus.
# Reuse the same continuous source frame and measured Fusion context. This is
# an isolated identity-scale prototype, not a promise for arbitrary Coord/scales.
prototype_front <- function(export, panel_id, start, end) {
  layout <- export$source_layouts[[panel_id]]
  context <- export$diagnostics[[panel_id]]$render_context
  id <- layout$seqs
  if (length(id) != 1L) stop("M0 needs exactly one visible instance")
  origin <- layout$sequence_layout$source_origin[match(id, layout$sequence_layout$accver)]
  if (!length(origin) || is.na(origin)) origin <- 1
  pos <- seq(start - 1, end, length.out = 80) - origin + 1
  if (isTRUE(layout$collinear)) {
    ref <- layout$sequence_reference
    xy <- cbind(ggchord:::ggchord_collinear_x(pos, id, ref), rep(unname(ref$y[id]), length(pos)))
    if (identical(layout$seq_axis, "vertical")) xy <- cbind(-xy[, 2], xy[, 1])
  } else xy <- ggchord:::ggchord_track_frame(pos, id, layout)$base
  bounds <- context$panel_bounds_mm
  data.frame(source_boundary = pos + origin - 1,
    x = bounds[1] + (xy[, 1] - context$x_range[1]) / diff(context$x_range) * diff(bounds[1:2]),
    y = bounds[3] + (xy[, 2] - context$y_range[1]) / diff(context$y_range) * diff(bounds[3:4]))
}

# Render panels and decorative leaders together on the current target device.
# No fake HSP table is created; biologically meaningful connections stay empty.
makeContent.ggchord_m0_prototype <- function(x) {
  result <- ggchord:::ggchord_build_fusion(x$composition)
  e <- result$export
  a <- prototype_front(e, "center", x$start, x$end)
  b <- prototype_front(e, "detail", x$start, x$end)
  ends <- c(1L, nrow(a))
  leaders <- grid::segmentsGrob(a$x[ends], a$y[ends], b$x[ends], b$y[ends],
    default.units = "mm", gp = grid::gpar(col = "#B45309", lty = 2, lwd = .8))
  # Leaders behind independently rebuilt panel foregrounds, with no filled band.
  # Composition panels have transparent backgrounds in this prototype.
  title <- grid::textGrob("pETDuet-1: source 101–1900 bp → independent linear detail",
    x=.5,y=.98,gp=grid::gpar(fontsize=13,fontface="bold"))
  grid::setChildren(x, grid::gList(leaders, result$grob, title))
}
