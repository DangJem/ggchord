# Direct source-coordinate checks, not package/release acceptance.
pkgload::load_all(quiet = TRUE)
source('examples/origin-detail.R')
fail <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
stopifnot(identical(cells$source_position, as.numeric(positions)),
  all(cells$source_verified), all(cells$source_length == source_length))
inverse <- map_alignment_position(cells, 'origin-window', positions, inverse = TRUE)
stopifnot(all(inverse$matched), identical(inverse$alignment_column, seq_along(positions)))
# Negative strand across the same true source origin, with a retained gap.
negative_positions <- rev(positions)
negative_text <- paste0(ggchord_alignment_complement(
  substring(genome$sequence, negative_positions, negative_positions), 'DNA'), collapse = '')
negative_text <- paste0(substr(negative_text, 1, 20), '-', substring(negative_text, 21))
negative_source <- source
negative_source$start <- region$end
negative_source$strand <- '-'
negative <- as_ggchord_alignment(c('origin-window' = negative_text), 'DNA',
  negative_source, circular = TRUE)
stopifnot(identical(negative$source_position[!negative$is_gap], as.numeric(negative_positions)),
  all(negative$source_verified), is.na(negative$source_position[21]))
negative_plot <- ggplot2::ggplot(negative, aes(alignment_column, accver, label = symbol)) +
  geom_residue(show.legend = FALSE) + ggplot2::theme_void()
negative_plot <- register_ggchord_alignment(negative_plot, negative)
negative_detail <- magnify_ggchord(whole, region, negative_plot)
export <- export_ggchord_fusion(origin_detail, width = 16, height = 4)
reverse_export <- export_ggchord_fusion(negative_detail, width = 16, height = 4)
for (d in list(export, reverse_export)) {
  stopifnot(nrow(d$focus_regions) == 2, nrow(d$leaders) == 2,
    all(d$focus_regions$wrap), all(d$focus_regions$region_row == 1),
    identical(d$focus_region_data, as.data.frame(region)) || identical(d$focus_region_data, region),
    nrow(d$connections) == 0,
    identical(d$leaders$boundary, c('start', 'end')),
    identical(d$leaders$source_position_boundary, c(region$start - 1, region$end)))
}
# No implicit wrapping, unknown length, or repeated source-position tours.
fail(as_ggchord_alignment(c(A = 'AC'), 'DNA', circular = TRUE), 'explicit source_data')
fail(as_ggchord_alignment(c(A = 'AC'), 'DNA',
  data.frame(accver = 'A', start = 3, strand = '+'), circular = TRUE), 'full source_length')
fail(as_ggchord_alignment(c(A = 'AAAAA'), 'DNA',
  data.frame(accver = 'A', start = 3, strand = '+', source_length = 4), circular = TRUE), 'revisit')
fail(as_ggchord_alignment(c(A = 'A'), 'DNA',
  data.frame(accver = 'A', start = 1, strand = '+', sequence = 'AC', source_length = 3), circular = TRUE), 'disagrees')
undeclared <- region |> select(-wrap)
fail(magnify_ggchord(whole, undeclared, characters), 'explicit wrap')
unknown <- characters
unknown$ggchord_alignment$source_length <- NULL
fail(export_ggchord_fusion(magnify_ggchord(whole, region, unknown)), 'matching declared')
# Focused genomic copies carry complete source length, not their local widths.
windows <- tibble(accver = genome$accver, start = c(region$start, 1), end = c(source_length, region$end))
focused <- focus_ggchord_data(genome, loci = windows)
# Different focus instances must not be silently spliced into one displayed instance.
local <- ggchord(focused$seq_data) + geom_seq() + coord_collinear()
fail(export_ggchord_fusion(magnify_ggchord(whole, region, local)), 'no fully visible')
# Precision at the maximum exactly representable source length, without full strings.
large <- as_ggchord_alignment(c(A = 'ACGT'), 'DNA', data.frame(accver = 'A',
  start = 2^53 - 2, strand = '+', source_length = 2^53 - 1), circular = TRUE)
stopifnot(identical(large$source_position, c(2^53 - 2, 2^53 - 1, 1, 2)))
ggplot2::ggsave('tools/long-term-audit/v0.27.0-origin-detail.png', origin_detail,
  width = 16, height = 4, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed', source = 'examples/plasmid/pUC19.fna',
  source_length = source_length, source_positions = positions, source_verified = TRUE,
  original_regions = nrow(region), drawn_pieces = nrow(export$focus_regions),
  external_leaders = nrow(export$leaders), negative_strand_with_gap = TRUE,
  implicit_wrap_rejected = TRUE, repeated_tour_rejected = TRUE,
  incomplete_instance_rejected = TRUE, exact_integer_boundary = TRUE),
  'tools/long-term-audit/v0.27.0-origin-results.json', pretty = TRUE, auto_unbox = TRUE)
cat('Explicit origin detail: passed\n')
