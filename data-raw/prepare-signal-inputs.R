# Rebuild two representations of actual pUC19 GC, not read coverage.
# Optional development dependencies: pkgload, jsonlite, digest, rtracklayer.
# No original input is edited. Use the packaged fixture ID explicitly after
# verifying its complete sequence against the original FASTA.
for (package in c('pkgload', 'jsonlite', 'digest', 'rtracklayer'))
  if (!requireNamespace(package, quietly = TRUE)) stop('Requires optional package: ', package)
pkgload::load_all(quiet = TRUE)
input <- 'examples/plasmid/pUC19.fna'
genome <- read_fasta_lengths(input, keep_sequence = TRUE)
stopifnot(nrow(genome) == 1L)
data(plasmid_example_pUC19)
stopifnot(identical(genome$sequence[[1L]], plasmid_example_pUC19$sequence[[1L]]))
original_id <- genome$accver
genome$accver <- plasmid_example_pUC19$accver
sequence <- strsplit(toupper(genome$sequence[[1L]]), '', fixed = TRUE)[[1L]]
stopifnot(all(sequence %in% c('A', 'C', 'G', 'T')))
starts <- seq.int(1L, length(sequence), by = 100L)
ends <- pmin(starts+99L, length(sequence))
values <- vapply(seq_along(starts), function(i)
  100*mean(sequence[seq.int(starts[i], ends[i])] %in% c('G', 'C')), numeric(1))
signal <- data.frame(accver = genome$accver, start = starts, end = ends, value = values)
output <- 'inst/extdata'
dir.create(output, showWarnings = FALSE)
bedgraph <- file.path(output, 'pUC19-gc-100.bedGraph')
header <- c('# Derived from unchanged examples/plasmid/pUC19.fna; actual GC, not read coverage.',
  '# Non-overlapping 100-bp windows; last window retains its actual length.',
  'track type=bedGraph name="pUC19 GC (%)"')
writeLines(header, bedgraph)
old_digits <- getOption('digits'); options(digits = 17)
utils::write.table(data.frame(signal$accver, signal$start-1, signal$end, signal$value),
  bedgraph, append = TRUE, sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
options(digits = old_digits)
ranges <- GenomicRanges::GRanges(seqnames = signal$accver,
  ranges = IRanges::IRanges(signal$start, signal$end), score = signal$value,
  seqlengths = stats::setNames(genome$length, genome$accver))
bigwig <- file.path(output, 'pUC19-gc-100.bw')
rtracklayer::export(ranges, bigwig, format = 'BigWig')
metadata <- list(source_file = input, source_sha256 = digest::digest(file = input, algo = 'sha256'),
  original_header = readLines(input, n = 1L, warn = FALSE), original_identifier = original_id,
  identifier = genome$accver, identifier_rule = 'Explicit packaged fixture ID pUC19 after checking the complete original FASTA sequence. Readers never rename identifiers.',
  sequence_length = genome$length, window_bp = 100L, intervals = nrow(signal),
  rule = 'Count actual G/C bases in each non-overlapping 100-bp window; divide by its actual width and multiply by 100. No smoothing, zero imputation or biological inference.',
  value_unit = 'GC percent', biological_interpretation = 'Sequence composition; not read coverage.',
  transformation = 'Explicit fixture identifier mapping recorded separately; bases unchanged. bedGraph stores zero-based half-open intervals. BigWig is exported by rtracklayer and stores floating-point signal values, not original observation counts.',
  generator = 'data-raw/prepare-signal-inputs.R', rtracklayer = as.character(utils::packageVersion('rtracklayer')),
  files = lapply(c(bedgraph, bigwig), function(path) list(file = path,
    sha256 = digest::digest(file = path, algo = 'sha256'))))
jsonlite::write_json(metadata, file.path(output, 'signal-sources.json'), pretty = TRUE, auto_unbox = TRUE)
cat('Rebuilt actual GC signal:', nrow(signal), 'windows,', genome$length, 'bp\n')
