# Finite v0.16.6 development checks; no testthat, build/check or release matrix.
# Rscript tools/validate-rc.R [layers|circular|linear|labels|output|data|all]
# New evidence only: rc-behavior-v0.16.6.csv. Unselected cases stay NOT_RUN.
pkgload::load_all(quiet = TRUE)
source("tools/rc-fixtures.R")
mode <- commandArgs(trailingOnly = TRUE)
if (!length(mode)) mode <- "all"
cases <- c("layers", "circular", "linear", "labels", "output", "data")
stopifnot(length(mode) == 1L, mode %in% c(cases, "all"))
result <- data.frame(case = cases, status = "NOT_RUN", detail = "Not selected",
  R = R.version.string, ggplot2 = as.character(utils::packageVersion("ggplot2")),
  platform = R.version$platform, device = "private PDF metrics; target PNG only in output case",
  revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE))
path <- "design-roadmap/audit-generated/rc-behavior-v0.16.6.csv"
save <- function(value) write.csv(value, path, row.names = FALSE)
snapshot <- function(p) lapply(p$layers, function(l) list(data = l$data,
  input = l$ggchord_input_data, mapping = l$mapping, params = l$ggchord_params,
  fixed = l$aes_params, id = l$ggchord_layer_id))
fails <- function(expr, text) {
  e <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(e, "error"), grepl(text, conditionMessage(e), fixed = TRUE))
}
local({
  original <- grDevices::dev.cur()
  grDevices::pdf(NULL, width = 8, height = 8, family = "Helvetica")
  owned <- grDevices::dev.cur()
  on.exit({
    if (owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
    if (original %in% grDevices::dev.list()) grDevices::dev.set(original)
  }, add = TRUE)
  f <- rc_fixtures()
  checks <- list(
    layers = function() {
      g <- f$genes[c("accver", "start", "end", "strand")]
      a <- g[g$accver %in% f$seqs$accver[1:2], ]; b <- g[!g$accver %in% f$seqs$accver[1:2], ]
      shared <- geom_gene(data = a, width = .04)
      p <- ggchord(f$seqs, f$links) + geom_seq(seq_radius = c(1, 1.1, 1, .95),
        seq_curvature = c(.8, 1, -.6, 1), seq_gap = .05, seq_orientation = c(1, -1, 1, -1)) +
        shared + geom_gene(data = b, width = .07) +
        geom_feature(data = a, shape = "block", position = position_strand(.15), feature_fill = "#578B74") +
        geom_feature(data = b, shape = "marker", position = position_strand(-.15), feature_fill = "#9774A0") +
        geom_link_ribbon(data = f$links[1:6, ], alpha = .4) +
        geom_link_ribbon(data = f$links[7:13, ], fill = "#CC9966") +
        geom_link_line(data = f$points[1:6, ], colour = "#216754") +
        geom_link_line(data = f$points[7:13, ], colour = "#8B4C72") +
        scale_gene_fill_manual(values = c("+" = "#123456", "-" = "#765432")) + coord_chord(rotation = 12)
      before <- snapshot(p); build <- ggplot2::ggplot_build(p)
      stopifnot(all(build$data[[2]]$gene_fill %in% c("#123456", "#765432")),
        all(build$data[[3]]$gene_fill %in% c("#123456", "#765432")))
      o <- export_ggchord_layout(p, include = c("gene", "feature", "ribbon", "link"), original_data = TRUE)
      for (component in c("gene", "feature", "ribbon", "link")) {
        stopifnot(length(unique(o[[component]]$layer_id)) == 2L)
        for (id in unique(o[[component]]$layer_id)) {
          inputs <- o$original_data[[id]][[unique(o[[component]]$component[o[[component]]$layer_id == id])]]
          stopifnot(setequal(o[[component]]$source_row[o[[component]]$layer_id == id],
            seq_len(nrow(inputs))))
        }
      }
      # One repeated object gets distinct owners; sibling builds do not mutate p.
      sibling <- ggchord(f$seqs) + geom_seq() + shared + shared + coord_collinear()
      out <- export_ggchord_layout(sibling, include = "gene")
      stopifnot(length(unique(out$gene$layer_id)) == 2L, identical(before, snapshot(p)))
      "Two gene/feature/ribbon/point layers retain source sets and styles; four-column genes honour user scale; reused layer IDs isolated; nondefault curve/radius/gap/orientation"
    },
    circular = function() {
      p <- f$circular(); l <- get_chord_layout(p)
      paths <- rc_leaders(l)
      stopifnot(nrow(paths) > 0L, rc_crossings(paths) == 0L)
      out <- export_ggchord_layout(p, include = c("feature", "labels", "restriction"))
      stopifnot(setequal(out$feature$source_row[out$feature$layer_id == unique(out$feature$layer_id)[1L]],
        seq_len(nrow(f$features))))
      # F54 regression, with all true feature/primer/site inputs present.
      inner <- export_ggchord_layout(f$circular(side = "inside"), include = "labels")
      d <- inner$labels
      d <- d[d$.component %in% "text" & d$annotation_class %in% "primer", ]
      stopifnot(setequal(d$source_row, seq_len(nrow(f$primers))), all(d$annotation_side == "inner"),
        all(sqrt(d$text_x^2 + d$text_y^2) < sqrt(d$anchor_x^2 + d$anchor_y^2)),
        any(inner$annotation_registry$kind == "primer_label" & inner$annotation_registry$side == "inner"))
      # Inside restriction placement must not be forced into the exterior fan.
      restricted <- export_ggchord_layout(f$circular(restriction_side = "inside"), include = "restriction")
      r <- restricted$restriction; r <- r[r$restriction_component %in% "label", ]
      rr <- restricted$annotation_registry
      stopifnot(nrow(r) == nrow(f$sites), all(sqrt(r$x^2 + r$y^2) < 1),
        sum(rr$kind == "restriction_label") == nrow(f$sites),
        all(rr$side[rr$kind == "restriction_label"] == "inner"))
      "pUC19 full features + six primers + three true cuts: sparse exterior leaders do not cross; F54 inside primers and inside cuts retain content/side; dense visual case is opt-in"
    },
    linear = function() {
      out <- export_ggchord_layout(f$linear, include = c("seq", "ribbon", "feature"), original_data = TRUE)
      stopifnot(nrow(out$ribbon) > 0L, all(out$connection_report$retained),
        setequal(out$sequence_layout$orientation, c(1, -1)))
      for (i in seq_len(nrow(f$windows))) {
        pos <- mean(c(f$cds$start[i], f$cds$end[i]))
        forward <- map_ggchord_position(f$linear, f$windows$locus_id[i], pos, source = TRUE)
        inverse <- map_ggchord_position(f$linear, f$windows$locus_id[i], forward$x, source = TRUE, inverse = TRUE)
        stopifnot(isTRUE(all.equal(inverse$position, pos)))
      }
      stopifnot(all(f$view$gene_data$.source_accver %in% f$cds$accver))
      "Source-backed pUC19/pETDuet AmpR exact CDS membership; focused windows, local reversal and source-coordinate inverse retain provenance; display links are not a realignment"
    },
    labels = function() {
      # Keep a small source-backed subset; one missing annotation stays missing.
      g <- f$genes[c(1, 4, 9, 17), ]; g$anno[1] <- NA_character_
      s <- f$seqs[f$seqs$accver %in% g$accver, ]
      for (mode in c("auto", "radial", "arc")) {
        p <- ggchord(s) + geom_seq() + geom_gene(data = g) +
          geom_gene_label_repel(data = g, label_layout = mode, label_wrap = NULL)
        l <- get_chord_layout(p)
        normalize <- function(x) gsub("[[:space:]]+", " ", trimws(gsub("\\n", " ", x, fixed = TRUE)))
        text <- normalize(l$gene_labels$label)
        stopifnot(all(normalize(na.omit(g$anno)) %in% text), !any(text %in% "NA"),
          !ggchord:::ggchord_label_box_conflicts(l$gene_labels,
            units_per_inch = l$text_units_per_inch, box_padding = 0),
          rc_crossings(rc_leaders(l)) == 0L)
      }
      "auto/radial/arc finite source-backed subsets retain nonmissing text, no literal NA, no label-box overlap or leader crossings at this size"
    },
    output = function() {
      p <- f$track; before <- snapshot(p)
      out <- export_ggchord_layout(p, include = "track")
      cache <- get_chord_layout(p, build = FALSE)
      caller <- grDevices::dev.cur(); context <- ggchord:::ggchord_render_context()
      for (size in list(c(6, 5), c(9, 7))) {
        d <- diagnose_ggchord(p, width = size[1], height = size[2], device = "png")
        stopifnot(d$measurement$status == "measured", all(is.finite(d$measurement$panel_mm)),
          grDevices::dev.cur() == caller, identical(ggchord:::ggchord_render_context(), context),
          nrow(d$track_registry) > 0L,
          all(vapply(d$track_registry$limits, function(x) identical(as.numeric(x), c(0, 1)), logical(1))))
      }
      preview <- view_ggchord(p, width = 6, height = 5, device = "png", viewer = "none")
      on.exit(unlink(c(preview, sub("[.]png$", ".html", preview))), add = TRUE)
      stopifnot(file.exists(preview), grDevices::dev.cur() == caller,
        identical(cache, get_chord_layout(p, build = FALSE)), identical(before, snapshot(p)),
        identical(out, export_ggchord_layout(p, include = "track")))
      fails(diagnose_ggchord(p + scale_track_value_continuous("gc", limits = c(.2, 1)),
        width = 6, height = 5, device = "png"), "baseline")
      stopifnot(grDevices::dev.cur() == caller, identical(ggchord:::ggchord_render_context(), context))
      "Real 100-bp DNA GC bins; two PNG measurement sizes and one preview; domain, caller device, source layers/cache and error context preserved"
    },
    data = function() {
      # Existing standard BLAST table (not a fabricated alignment).
      # Locate the declared source rather than infer statistics for six-column links.
      standard <- list.files("examples/blastn-standard", pattern = "[.]tsv$", full.names = TRUE)
      stopifnot(length(standard) > 0L)
      read <- read_blast(standard[1L])
      stopifnot(nrow(read) > 0L, all(c("qaccver", "saccver", "pident") %in% names(read)))
      adapted <- as_ggchord_link(data.frame(query = f$links$qaccver,
        subject = f$links$saccver, from = f$links$qstart, to = f$links$qend,
        target_from = f$links$sstart, target_to = f$links$send),
        mapping = c(qaccver = "query", saccver = "subject", qstart = "from", qend = "to",
          sstart = "target_from", send = "target_to"))
      stopifnot(identical(adapted$qstart, f$links$qstart), identical(adapted$send, f$links$send))
      model <- build_ggchord_cds(transform(f$cds, phase = 0), genetic_code = 11)
      residues <- data.frame(feature_id = "display_residues", parent_id = f$cds$feature_id[1], start = 2, end = 20)
      projected <- project_ggchord_features(residues, model)
      back <- project_ggchord_features(projected, model, from = "genome")
      stopifnot(back$start == 2, back$end == 20, projected$strand == f$cds$strand[1])
      fails(geom_gene(gene_offset = .1), "position_strand")
      fails(geom_feature(feature_shape = "block"), "shape")
      fails(geom_gene_label_repel(label_rotation = 90), "geom_gene_label(label_rotation")
      "Existing standard BLAST import, explicitly renamed six-column non-BLAST roles, true CDS residue projection round-trip and representative migration errors; no native aa axis claim"
    }
  )
  for (case in if (mode == "all") cases else mode) {
    i <- match(case, result$case)
    value <- tryCatch(checks[[case]](), error = identity)
    result$status[i] <- if (inherits(value, "error")) "FAIL" else "PASS_FINITE"
    result$detail[i] <- if (inherits(value, "error")) conditionMessage(value) else value
    save(result); cat(case, result$status[i], result$detail[i], "\n")
    if (inherits(value, "error")) stop(value)
  }
})
