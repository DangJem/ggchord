# Focused checks for decorative linear sequence-end breaks.
# Run from the repository root: Rscript tools/validate-seq-break.R
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver = c("A", "B"), length = c(1000, 800))
focused <- focus_ggchord_data(seqs,
  loci = data.frame(accver = c("A", "B"), start = c(200, 1), end = c(700, 500)))
ref <- list(lens = stats::setNames(focused$seq_data$length, c("A", "B")),
  orientation = c(A = 1, B = -1), x_scale = 2 / max(focused$seq_data$length),
  x_offset = c(A = -.2, B = .3), y = c(A = 0, B = -.6))
params <- list(break_width = .02, break_height = .04)
geometry <- ggchord:::ggchord_seq_break_geometry(focused$seq_data, params,
  list(collinear = TRUE, sequence_reference = ref))
ends <- unique(geometry[c("accver", "source_end", "display_end")])
stopifnot(nrow(ends) == 3L,
  all(c("left", "right") %in% ends$source_end[ends$accver == "A"]),
  all(ends$source_end[ends$accver == "B"] == "right"),
  all(ends$display_end[ends$accver == "B"] == "left"),
  length(unique(geometry$group)) == 2L * nrow(ends),
  all(is.finite(geometry$x)), all(is.finite(geometry$y)))
b_mark <- geometry[geometry$accver == "B", , drop = FALSE]
stopifnot(abs(mean(range(b_mark$x)) - ref$x_offset["B"]) < .Machine$double.eps^0.5)
complete <- ggchord:::ggchord_seq_break_geometry(seqs, params,
  list(collinear = TRUE, sequence_reference = ref))
stopifnot(nrow(complete) == 0L)
fallback <- focused$seq_data
fallback$.truncated_left <- fallback$.truncated_right <- NULL
fallback_geometry <- ggchord:::ggchord_seq_break_geometry(fallback, params,
  list(collinear = TRUE, sequence_reference = ref))
stopifnot(identical(geometry, fallback_geometry))
error <- tryCatch(ggchord:::ggchord_seq_break_geometry(focused$seq_data, params,
  list(collinear = FALSE)), error = identity)
stopifnot(inherits(error, "error"), grepl("coord_collinear", conditionMessage(error)))
layer <- geom_seq_break(aes(colour = accver))
stopifnot(identical(layer$ggchord_type, "seq_break"), is.null(layer$aes_params$colour))
cat("Linear sequence break decorations passed.\n")
