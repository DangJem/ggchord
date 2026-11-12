# M1 direct checks; no package test suite or release device matrix.
pkgload::load_all(quiet = TRUE)
library(dplyr)
seqs <- seq_data_example |> select(accver, length)
features <- gene_data_example |> mutate(type = 'gene')
hsp <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hsp, accver = qaccver, start = pmin(qstart, qend), end = pmax(qstart, qend)),
  transmute(hsp, accver = saccver, start = pmin(sstart, send), end = pmax(sstart, send))
)
# Deliberate user metadata collision: retain it only in the raw selection table.
regions$region_row <- c('query window', 'subject window')
regions$instance_note <- regions$detail_instance_note <- c('source note', 'detail note')
local <- focus_ggchord_data(seqs, feature_data = features,
  ribbon_data = hsp, loci = regions)
original_features <- features[local$feature_data$.source_row, ]
stopifnot(nrow(local$ribbon_data) == 1L,
  identical(local$feature_data$strand, original_features$strand),
  identical(local$feature_data$.source_accver, original_features$accver),
  all(local$feature_data$.source_start == original_features$start),
  all(local$feature_data$.source_end == original_features$end))
# A cropped query window must also project the subject end of the real HSP.
# Compare original-coordinate progress, allowing the documented integer rounding.
cropped_regions <- regions
cropped_regions$start[1] <- regions$start[1] + floor((regions$end[1]-regions$start[1])/4)
cropped <- focus_ggchord_data(seqs, ribbon_data = hsp, loci = cropped_regions)$ribbon_data
stopifnot(nrow(cropped) == 1L, cropped$.alignment_clipped,
  cropped$.source_row == 1L, cropped$pident == hsp$pident)
progress <- list()
for (end in c('q', 's')) {
  original <- unlist(hsp[paste0(end, c('start', 'end'))], use.names = FALSE)
  actual <- unlist(cropped[paste0(end, c('start', 'end'))], use.names = FALSE) +
    cropped[[paste0('.', end, 'source_origin')]] - 1
  stopifnot(all(unlist(cropped[paste0('.source_', end, c('start', 'end'))],
    use.names = FALSE) == original), sign(diff(actual)) == sign(diff(original)))
  progress[[end]] <- (actual-original[1])/diff(original)
}
stopifnot(progress$q[1]>0, max(abs(progress$q-progress$s)) <=
  1/min(abs(hsp$qend-hsp$qstart),abs(hsp$send-hsp$sstart)))
source <- ggchord(seqs, feature_data = features, ribbon_data = hsp) +
  geom_seq(show.legend = FALSE) + geom_gene(show.legend = FALSE) +
  geom_link_ribbon(show.legend = FALSE)
detail <- ggchord(local$seq_data, feature_data = local$feature_data,
  ribbon_data = local$ribbon_data) + geom_seq(show.legend = FALSE) +
  geom_gene(show.legend = FALSE) + geom_link_ribbon(show.legend = FALSE) + coord_collinear()
check <- function(x) {
  e <- export_ggchord_fusion(x, 12, 8)
  stopifnot(identical(e$focus_region_data, regions), nrow(e$focus_regions) == 2L,
    nrow(e$leaders) == 4L, nrow(e$connections) == 0L, length(e$relationship_data) == 0L,
    identical(e$focus_regions$region_row, 1:2))
  for (i in 1:2) {
    fronts <- filter(e$focus_fronts, region_row == i)
    for (panel in c('center', 'detail')) {
      stopifnot(isTRUE(all.equal(range(fronts$position_boundary[fronts$panel_id == panel]),
        c(regions$start[i]-1, regions$end[i]))))
      # Independent anchor check: interpolate the actual rendered backbone,
      # using input windows and physical panel bounds, not fusion_front().
      plot <- x$panels[[panel]]$plot
      layout <- e$source_layouts[[panel]]
      id <- regions$accver[i]
      backbone <- layout$seq_arcs[[id]]
      input <- plot$ggchord$data$seq_data
      row <- match(id, input$accver)
      origin <- if ('.source_origin' %in% names(input)) input$.source_origin[row] else 1
      fraction <- (c(regions$start[i]-1,regions$end[i])-origin+1)/input$length[row]
      if (layout$sequence_reference$orientation[id]<0) fraction <- 1-fraction
      frame <- e$diagnostics[[panel]]$render_context
      expected <- matrix(NA_real_,2,2,dimnames=list(NULL,c('x','y')))
      for (axis in c('x','y')) {
        values <- stats::approx(seq(0,1,length.out=nrow(backbone)),backbone[[axis]],
          xout=fraction)$y
        scale <- plot$scales$get_scales(axis)
        if (!is.null(scale)) values <- scale$transform(values)
        limits <- frame[[paste0(axis,'_range')]]
        bounds <- frame$panel_bounds_mm[if(axis=='x') 1:2 else 3:4]
        expected[,axis] <- bounds[1]+(values-limits[1])/diff(limits)*diff(bounds)
      }
      front <- fronts[fronts$panel_id==panel,]
      actual <- as.matrix(front[c(1,nrow(front)),c('x','y')])
      # This is a physical attachment tolerance, not frozen device coordinates.
      stopifnot(all(sqrt(rowSums((actual-expected)^2)) < x$focus$width/10))
    }
  }
  e
}
chord <- magnify_ggchord(source + coord_chord(), regions, detail, side = 'bottom')
e <- check(chord)
linear <- magnify_ggchord(source + coord_collinear(), regions, detail, side = 'bottom')
check(linear) |> invisible()
check(magnify_ggchord(source + coord_collinear(seq_axis = 'vertical') +
  ggplot2::scale_x_reverse(), regions, detail, side = 'bottom')) |> invisible()
# Distinct selectors are resolved independently for each selected source.
repeated <- ggchord(seqs) + geom_seq(seq_order = setNames(
  c(seqs$accver[1], seqs$accver[2], seqs$accver[1]), c('q1', 's', 'q2')),
  show.legend = FALSE) + coord_collinear()
selected <- regions
selected$instance <- c('q2', 's')
selected$detail_instance <- selected$accver
routed <- export_ggchord_fusion(magnify_ggchord(repeated, selected, detail, side = 'bottom'))
stopifnot(identical(routed$focus_regions$source_instance, c('q2', 's')))
ambiguous <- tryCatch(export_ggchord_fusion(magnify_ggchord(repeated, regions, detail)), error = identity)
stopifnot(inherits(ambiguous, 'error'), grepl('ambiguous', conditionMessage(ambiguous)))
outside <- regions; outside$end[2] <- seqs$length[2] + 1
rejected <- tryCatch(export_ggchord_fusion(magnify_ggchord(source, outside, detail)), error = identity)
stopifnot(inherits(rejected, 'error'))
one <- magnify_ggchord(source + coord_collinear(), regions[1, ], detail, side = 'bottom')
stopifnot(nrow(export_ggchord_fusion(one)$leaders) == 2L)
none <- magnify_ggchord(source + coord_collinear(), regions, detail,
  side = 'bottom', connector = 'none')
stopifnot(nrow(export_ggchord_fusion(none)$leaders) == 0L)
bad <- regions; bad$end[2] <- Inf
err <- tryCatch(magnify_ggchord(source, bad, detail), error = identity)
stopifnot(inherits(err, 'error'))
# No false HSPs from selection; source fields refer to the original observation.
stopifnot(identical(local$ribbon_data$pident, hsp$pident),
  local$ribbon_data$.source_row == 1L,
  local$ribbon_data$qstart+local$ribbon_data$.qsource_origin-1 == hsp$qstart,
  local$ribbon_data$send+local$ribbon_data$.ssource_origin-1 == hsp$send)
ggplot2::ggsave('tools/long-term-audit/v0.23.0-magnify-pairs.png', linear,
  width = 12, height = 8, dpi = 120, bg = 'white')
jsonlite::write_json(list(version = '0.23.0.9000', source_sequences = nrow(seqs),
  selected_regions = nrow(regions), local_hsp = nrow(local$ribbon_data),
  decorative_leaders = nrow(e$leaders), exact_source_boundaries = TRUE,
  raw_metadata_preserved = TRUE, chord_and_horizontal_vertical_collinear = TRUE,
  reflected_scale = TRUE, single_region_compatible = TRUE,
  per_record_feature_direction = TRUE, cropped_hsp_source_projection = TRUE,
  independent_backbone_anchors_mm = TRUE,
  no_release_matrix = TRUE), 'tools/long-term-audit/v0.23.0-magnify-pairs.json',
  pretty = TRUE, auto_unbox = TRUE)
cat('PASS: M1 paired source regions, real HSP, native frames and raw metadata\n')
