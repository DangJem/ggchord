# Directed local reader check; optional mature decoders must already exist.
# BCF_PYTHON explicitly selects a prepared environment; no installation here.
pkgload::load_all(quiet = TRUE)
for (package in c('VariantAnnotation', 'StructuralVariantAnnotation'))
  if (!requireNamespace(package, quietly = TRUE)) stop('This check needs optional ', package)
python <- Sys.getenv('BCF_PYTHON', 'python3')
fail <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, 'error'), grepl(pattern, conditionMessage(error)))
}
path <- system.file('extdata', 'gridss.vcf', package = 'StructuralVariantAnnotation')
imported <- read_vcf(path, structural = TRUE, source_file = TRUE)
original <- VariantAnnotation::readVcf(path)
variant <- imported$features$variant
stopifnot(nrow(variant) == nrow(original), nrow(imported$links) == 2,
  identical(variant$feature_id, names(original)),
  all(variant$.source_file == path),
  identical(VariantAnnotation::info(imported$report$records[[1L]]$parsed), VariantAnnotation::info(original)),
  identical(VariantAnnotation::geno(imported$report$records[[1L]]$parsed), VariantAnnotation::geno(original)),
  !any(c('pident', 'Orientation', 'read_support') %in% names(imported$links)),
  any(grepl('unpaired', imported$report$records[[1L]]$warnings)),
  nrow(variant) > length(unique(c(imported$links$q_source_record, imported$links$s_source_record))))
# Same genuine source records transcode to modern BCF without editing evidence.
# Python script reads one JSON config so paths never become code.
temporary <- tempfile('ggchord-bcf-check-'); dir.create(temporary)
config <- file.path(temporary, 'input.json'); script <- file.path(temporary, 'encode.py')
binary <- file.path(temporary, 'source.bcf')
jsonlite::write_json(list(source = path, destination = binary), config, auto_unbox = TRUE)
writeLines(c('import json,sys,pysam',
  'with open(sys.argv[1]) as h: config=json.load(h)',
  'with pysam.VariantFile(config["source"]) as source:',
  '    with pysam.VariantFile(config["destination"], "wb", header=source.header) as target:',
  '        for record in source: target.write(record)'), script)
stopifnot(identical(system2(python, c(shQuote(script), shQuote(config))), 0L))
bcf <- read_vcf(binary, structural = TRUE, python = python, source_file = TRUE)
stopifnot(identical(bcf$links[names(imported$links)[names(imported$links) != '.source_file']],
                    imported$links[names(imported$links)[names(imported$links) != '.source_file']]),
  identical(bcf$features$variant$ALT, variant$ALT),
  !'.source_line' %in% names(bcf$features$variant),
  bcf$report$records[[1L]]$binary_bridge$package == 'pysam',
  identical(VariantAnnotation::geno(bcf$report$records[[1L]]$parsed), VariantAnnotation::geno(original)))
unlink(temporary, recursive = TRUE)
# Explicit synthetic format controls: not scientific or read-support evidence.
fixture <- tempfile(fileext = '.vcf')
header <- c('##fileformat=VCFv4.3', '##contig=<ID=track,length=1000>',
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
  '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Position interval">',
  '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="End interval">',
  '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample')
writeLines(c(header,
  'track\t10\t.\tA\tC,G\t.\t.\t.\tGT\t1/2',
  'track\t30\tsv\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=50;IMPRECISE;CIPOS=-2,3;CIEND=-4,5\tGT\t0/1'), fixture)
controlled <- read_vcf(fixture, structural = TRUE)
stopifnot(nrow(controlled$features$variant) == 2,
  is.na(controlled$features$variant$feature_id[1L]),
  identical(controlled$features$variant$ALT[[1L]], c('C', 'G')),
  controlled$features$variant$start[2L] == controlled$features$variant$end[2L],
  nrow(controlled$links) == 1,
  controlled$links$qend > controlled$links$qstart,
  controlled$links$send > controlled$links$sstart,
  all(is.na(controlled$features$variant$QUAL)))
outside <- gsub('CIPOS=-2,3', 'CIPOS=-40,3', readLines(fixture), fixed = TRUE)
writeLines(outside, fixture)
edge <- read_vcf(fixture, structural = TRUE)
stopifnot(any(!edge$breakends$.plottable), nrow(edge$links) == 0,
  nrow(edge$report$records[[1L]]$unplaced_links) == 1)
# Lack of full source lengths never becomes maximum observed position.
no_length <- gsub(',length=1000', '', readLines(fixture), fixed = TRUE)
writeLines(no_length, fixture)
unknown <- read_vcf(fixture)
stopifnot(all(is.na(unknown$report$records[[1L]]$sequence_info$length)))
writeLines(header, fixture)
stopifnot(nrow(read_vcf(fixture, structural = TRUE)$features$variant) == 0)
writeLines(c(header, 'track\t0\tterminal\tA\tC\t.\t.\t.\tGT\t0/1'), fixture)
telomeric <- read_vcf(fixture)
stopifnot(nrow(telomeric$features$variant) == 0,
  nrow(telomeric$report$records[[1L]]$unplaced_features) == 1)
# IDs label records, not unique biological calls. Default retains duplicates;
# explicit structural decode rejects ambiguous mate namespaces.
writeLines(c(header,
  'track\t10\trepeated\tA\tC\t.\t.\t.\tGT\t0/1',
  'track\t20\trepeated\tA\tG\t.\t.\t.\tGT\t0/1'), fixture)
stopifnot(nrow(read_vcf(fixture)$features$variant) == 2)
fail(read_vcf(fixture, structural = TRUE), 'unique record IDs')
unlink(fixture)
# Real caller's header supplies full lengths; no last-breakpoint inference.
known <- imported$report$records[[1L]]$sequence_info
known <- known[known$accver %in% c('chr1', 'chr12'), c('accver', 'length')]
stopifnot(nrow(known) == 2, all(is.finite(known$length)))
# Explicitly select the interchromosomal relationship for this overview.
# The intrachromosomal pair remains intact in the returned links/source.
plot_links <- imported$links[imported$links$qaccver != imported$links$saccver, , drop = FALSE]
points <- transform(plot_links, qpos = qstart, spos = sstart)
stopifnot(all(points$qstart == points$qend), all(points$sstart == points$send))
plot <- ggchord(known) + geom_seq(show.legend = FALSE) + geom_seq_label() +
  geom_link_line(data = points, colour = '#447C8D', linewidth = .5, alpha = .7) +
  coord_collinear() + labs(title = 'Observed structural breakpoints',
    subtitle = 'One interchromosomal pair; all eight source records retained in report')
ggplot2::ggsave('tools/long-term-audit/v0.27.0-variant-links.png', plot,
  width = 8, height = 4, dpi = 120, bg = 'white')
jsonlite::write_json(list(status = 'passed', source = path,
  package = 'StructuralVariantAnnotation', package_version = as.character(utils::packageVersion('StructuralVariantAnnotation')),
  source_sha256 = digest::digest(file = path, algo = 'sha256'),
  original_records = nrow(variant), reciprocal_pairs = nrow(imported$links),
  explicitly_selected_overview_pairs = nrow(plot_links),
  source_info_and_genotypes_preserved = TRUE,
  outside_reference_uncertainty_reported = TRUE, empty_and_telomeric_preserved = TRUE, unpaired_retained_in_source = TRUE,
  bcf_roundtrip = TRUE, multiallelic_records_not_split = TRUE,
  uncertainty_not_collapsed = TRUE, unknown_lengths_remain_unknown = TRUE,
  no_fabricated_support_identity_or_mate = TRUE),
  'tools/long-term-audit/v0.27.0-variant-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('Variant adapters: passed\n')
