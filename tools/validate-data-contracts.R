# Targeted v0.16.2 checks. Synthetic fixtures, no package test/check or device matrix.
pkgload::load_all(quiet = TRUE)
local({
  directory <- tempfile('ggchord-data-')
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  error <- function(expr, pattern) {
    value <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(value, 'error'), grepl(pattern, conditionMessage(value), fixed = TRUE))
  }
  seqs <- data.frame(accver = c('A', 'B'), length = 1000)
  links <- data.frame(qaccver = 'A', saccver = 'B', qstart = c(1, 101),
    qend = c(1000, 200), sstart = c(1, 101), send = c(1000, 200),
    length = c(1100, 100), pident = 95)
  original <- links
  validation <- validate_ggchord_data(seqs, links)
  stopifnot(!any(validation$warnings$category == 'overlapping'))
  stopifnot(nrow(deduplicate_ggchord_ribbons(links, by = 'overlap')$data) == 2L)
  reverse <- transform(links[2L, ], qstart = 200, qend = 101, sstart = 200, send = 101)
  reverse2 <- transform(reverse, qstart = 199, qend = 100, sstart = 199, send = 100)
  stopifnot(any(validate_ggchord_data(seqs, rbind(reverse, reverse2))$warnings$category == 'overlapping'),
    nrow(deduplicate_ggchord_ribbons(rbind(reverse, reverse2), by = 'overlap')$data) == 1L)
  opposite <- transform(links[2L, ], sstart = 200, send = 101)
  stopifnot(nrow(deduplicate_ggchord_ribbons(rbind(links[2L, ], opposite), by = 'overlap')$data) == 2L)
  single <- merge_ggchord_ribbons(links[1L, ])
  stopifnot(single$data$length == 1100, !single$data$.length_estimated)
  merged <- merge_ggchord_ribbons(transform(links, qend = c(100, 200),
    send = c(100, 200), length = c(120, 130)))
  stopifnot(nrow(merged$data) == 1L, merged$data$.length_estimated,
    identical(merged$report$source_length[[1L]], c(120, 130)), identical(links, original))
  separate <- merge_ggchord_ribbons(transform(links, qaccver = c('A', 'C'), tag = 'x'), group_by = 'tag')
  stopifnot(nrow(separate$data) == 2L)
  crossing <- data.frame(qaccver = 'A', saccver = 'B', qstart = c(1, 11),
    qend = c(10, 20), sstart = c(11, 1), send = c(20, 10), length = 10, pident = 95)
  stopifnot(nrow(merge_ggchord_ribbons(crossing)$data) == 2L)
  bundle <- bundle_ggchord_ribbons(rbind(links[2L, ], reverse), seqs, bins = 1, weight = 'count')
  stopifnot(nrow(bundle$data) == 1L, bundle$data$.bundle_n == 2L,
    bundle$report$source_rows == '1,2', bundle$data$qstart < bundle$data$qend)
  points <- data.frame(qaccver = 'A', saccver = 'B', qstart = c(101, 102),
    qend = c(101, 102), sstart = c(900, 899), send = c(900, 899),
    qstrand = '+', sstrand = '-', direction = 'reverse', orientation_known = TRUE,
    length = 1, pident = 100)
  prepared <- prepare_ggchord_data(seqs, points,
    steps = list(clean = list(), merge = list(), bundle = list(bins = 1, weight = 'count')))
  stopifnot(nrow(prepared$ribbon_data) == 2L, all(prepared$ribbon_data$direction == 'reverse'),
    all(prepared$ribbon_data$.bundle_n == 1L))
  distinct_points <- rbind(points[1, ], transform(points[1, ], sstrand = '+', direction = 'same'))
  stopifnot(nrow(deduplicate_ggchord_ribbons(distinct_points)$data) == 2L)
  density_input <- transform(links[c(2, 2), ], orientation_known = c(TRUE, FALSE))
  density_build <- ggplot2::ggplot_build(ggchord(seqs, density_input, validate = 'none') +
    geom_seq() + stat_ribbon_density(aes(ribbon_alpha = after_stat(bundle_n)),
      bins = 1, weight = 'count') + coord_collinear())
  stopifnot(all(density_build$data[[2L]]$bundle_n == 1L),
    all(is.finite(density_build$data[[2L]]$ribbon_alpha)))
  unknown <- transform(links[1, ], qaccver = NA_character_, saccver = '')
  kept <- clean_ggchord_data(seqs, unknown, unknown_id = 'keep')
  stopifnot(identical(kept$ribbon_data, unknown), nrow(kept$report) == 2L)

  path <- file.path(directory, 'input.blast')
  path2 <- file.path(directory, 'other.blast')
  line <- '001\tNA\t95\t10\t0\t0\t1\t10\t10\t1\t0\t1'
  writeLines(c('# source comment', line), path)
  writeLines(line, path2)
  read <- read_blast(file = file.path(directory, '*.blast'), files = path, source_file = TRUE)
  stopifnot(nrow(read) == 2L, identical(read$qaccver, c('001', '001')),
    identical(read$saccver, c('NA', 'NA')), identical(read$.source_line, c(2L, 1L)),
    identical(attr(read, 'ggchord_import')$rows_per_file, c(1L, 1L)))
  bad <- file.path(directory, 'bad.blast')
  writeLines(c(line, '001\tNA\t95'), bad)
  error(read_blast(bad), 'inconsistent field counts')
  error(read_blast(file = NA_character_), 'non-missing')
  error(read_blast(file = directory), 'file not found')
  error(filter_ggchord_ribbons(links, min_pident = c(10, 20)), 'one finite number')
  error(find_common_features(c(A = 'ACGT'), min_protein_identity = c(.8, .9)), 'each be one finite')
  error(find_primer_bindings(c(A = 'ACGTACGT'), primers = c(test = 'ACGTACGT'),
    min_annealed_bases = 8.5), 'integer of at least 8')
  invalid <- file.path(directory, 'encoding.bed')
  writeBin(as.raw(c(0xff, 0x09, 0x30, 0x09, 0x31, 0x0a)), invalid)
  error(read_bed(invalid), 'invalid UTF-8')
  bed <- file.path(directory, 'one.bed')
  writeLines('A\t0\t1', bed)
  stopifnot(read_bed(bed)$start == 1, read_bed(bed)$end == 1)
  writeLines('A\t1\t1', bed)
  error(read_bed(bed), 'zero-width')

  patterns <- data.frame(enzyme = 'custom', motif = 'ACGT', blocked_by = 'Dcm',
    methylation_motif = 'CCWGG')
  linear <- find_restriction_sites(c(A = 'AGGACGTCC'), database = patterns,
    circular = FALSE, methylation = 'dam_dcm')
  circular <- find_restriction_sites(c(A = 'AGGACGTCC'), database = patterns,
    circular = TRUE, methylation = 'dam_dcm')
  stopifnot(all(linear$methylation_status == 'unblocked'),
    all(circular$methylation_status == 'blocked'))
  feature <- data.frame(accver = 'A', start = factor('100', levels = c('100','300')),
    end = factor('300', levels = c('100','300')), strand = '+', name = 'primer')
  source <- feature
  labels <- export_ggchord_layout(ggchord(seqs[1L, ]) + geom_seq() +
    geom_primer_label_repel(data = feature) + coord_circular(), include = 'labels')$labels
  labels <- labels[!is.na(labels$primer_start), ]
  stopifnot(all(labels$primer_start == 100), all(labels$primer_end == 300),
    identical(feature, source))
  # Explicit limitation: metadata does not change the standard-code matcher.
  db <- data.frame(name = 'nonstandard', sequence = 'ATGTGATAA', reference_protein = 'MW', geneticCode = 2)
  stopifnot(nrow(find_common_features(c(A = 'ATGTGATAA'), db, mode = 'protein', circular = FALSE)) == 0L,
    nrow(find_common_features(c(A = 'ATGTGATAA'), db, mode = 'dna', circular = FALSE)) > 0L)
})
cat('v0.16.2 overlap, merge, bundle, import, methylation and primer contracts passed.\n')
