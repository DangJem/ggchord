# Optional mature-tool bridge, not a new BAM/CRAM reader or alignment engine.
.libPaths(c(Sys.getenv('ADAPTER_R_LIB', '/tmp/ggchord-adapter-r-lib'), .libPaths()))
pkgload::load_all(quiet = TRUE)
library(dplyr)
for (package in c('GenomicAlignments', 'Rsamtools', 'GenomicRanges', 'jsonlite', 'digest'))
  if (!requireNamespace(package, quietly = TRUE)) stop('Optional dependency required: ', package)
bam <- system.file('extdata', 'sm_treated1.bam', package = 'GenomicAlignments')
# Original fixed Bioconductor source bytes, not a generated format fixture.
expected <- 'ac2681fae928f670d597e0f95854602e1daea9aa34947e07335bd88a8be9ee53'
stopifnot(digest::digest(file = bam, algo = 'sha256') == expected)
param <- Rsamtools::ScanBamParam(
  flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE, isDuplicate = FALSE), mapqFilter = 0L)
reads <- GenomicAlignments::readGAlignments(bam, param = param)
coverage <- GenomicAlignments::coverage(reads, drop.D.ranges = TRUE)
blocks <- unlist(GenomicAlignments::grglist(reads, drop.D.ranges = TRUE), use.names = FALSE)
stopifnot(sum(sum(coverage)) == sum(IRanges::width(blocks)), any(GenomicAlignments::njunc(reads) > 0L))
# Convert the decoder's run lengths only; retain true zero runs as zero.
signal <- bind_rows(lapply(names(coverage), function(id) {
  runs <- coverage[[id]]
  ends <- cumsum(S4Vectors::runLength(runs))
  data.frame(accver = id, start = c(1, head(ends, -1L) + 1), end = ends,
    value = S4Vectors::runValue(runs))
}))
path <- tempfile(fileext = '.bedGraph')
encoded <- signal |> transmute(accver, start = start - 1, end, value)
write.table(encoded, path, sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
decoded <- read_bedgraph(path)
stopifnot(all(vapply(names(signal), function(key) all(decoded[[key]] == signal[[key]]), logical(1))),
  sum((decoded$end-decoded$start+1)*decoded$value) == sum(sum(coverage)))
unlink(path)
# Explicit CRAM no-reference encoding control from these real BAM records.
# No external reference/URL is requested; this does not certify reference-based
# CRAM files without their actual matching reference. Preserve original header.
temporary <- tempfile('ggchord-cram-'); dir.create(temporary)
roundtrip <- file.path(temporary, 'roundtrip.bam')
script <- c(
  'import sys,json,pysam,os',
  'from pathlib import Path',
  'cfg=json.loads(Path(sys.argv[1]).read_text())',
  'source=cfg["source"]; target=cfg["target"]; cram=target+".cram"',
  'def fingerprint(r):',
  ' return [r.query_name,r.flag,r.reference_id,r.reference_start,r.mapping_quality,r.cigarstring,r.next_reference_id,r.next_reference_start,r.template_length,r.query_sequence,list(r.query_qualities) if r.query_qualities is not None else None,sorted(r.get_tags())]',
  'with pysam.AlignmentFile(source,"rb") as src:',
  ' original=[fingerprint(r) for r in src.fetch(until_eof=True)]',
  ' header=src.header.to_dict()',
  'with pysam.AlignmentFile(source,"rb") as src, pysam.AlignmentFile(cram,"wc",template=src,format_options=[b"no_ref=1"]) as out:',
  ' for r in src.fetch(until_eof=True): out.write(r)',
  'with pysam.AlignmentFile(cram,"rc") as src, pysam.AlignmentFile(target,"wb",template=src) as out:',
  ' observed=[]',
  ' for r in src.fetch(until_eof=True): observed.append(fingerprint(r));out.write(r)',
  ' assert observed==original',
  ' assert src.header.to_dict()["SQ"]==header["SQ"]',
  'Path(sys.argv[2]).write_text(json.dumps({"records":len(original),"all_record_fields_and_tags_equal":True,"reference_dictionary_equal":True,"pysam":pysam.__version__,"samtools":pysam.__samtools_version__,"cram_mode":"explicit no_ref=1 re-encoding, not an original CRAM observation"}))')
bridge <- ggchord_parse_external_source(script, list(source = bam, target = roundtrip),
  Sys.getenv('BCF_PYTHON', '/tmp/ggchord-hts-python-venv/bin/python'), 'coverage prototype')
converted <- GenomicAlignments::readGAlignments(roundtrip, param = param)
stopifnot(identical(coverage, GenomicAlignments::coverage(converted, drop.D.ranges = TRUE)))
unlink(temporary, recursive = TRUE)
# Source lengths are the original BAM @SQ dictionary; never maximum read ends.
header <- Rsamtools::scanBamHeader(bam)[[1L]]
seqs <- data.frame(accver = names(header$targets), length = as.numeric(header$targets))
selected <- signal |> filter(accver == 'chr2L', start <= 10000) |>
  mutate(end = pmin(end, 10000))
# Explicitly show a source window; do not label it a complete chromosome.
focused <- focus_ggchord_data(seqs, loci = data.frame(accver = 'chr2L', start = 1, end = 10000))
plot <- ggchord(focused$seq_data) + geom_seq(show.legend = FALSE) +
  geom_track(aes(xmin = start - 1, xmax = end, ymin = 0, ymax = value), data = selected,
    geom = 'rect', fill = '#387D88', track = position_track('aligned_records', width = unit(12, 'mm'))) +
  scale_track_value_continuous('aligned_records', unit = 'Aligned records / base') +
  coord_collinear() + labs(title = 'Existing Drosophila RNA-seq alignments',
    subtitle = 'chr2L 1–10,000; reference skips and deletions excluded; zero runs retained') +
  theme(plot.title.position = 'plot', plot.subtitle = element_text(hjust = .5))
invisible(export_ggchord_layout(plot))
ggplot2::ggsave('tools/long-term-audit/v0.27.0-read-coverage.png', plot,
  width = 9, height = 4, dpi = 110, bg = 'white')
jsonlite::write_json(list(status = 'passed',
  source = list(package = 'GenomicAlignments', version = as.character(packageVersion('GenomicAlignments')),
    file = 'extdata/sm_treated1.bam', sha256 = expected,
    url = 'https://raw.githubusercontent.com/Bioconductor/GenomicAlignments/34a743a22b7ce4f754c3962bb06ae872581905b2/inst/extdata/sm_treated1.bam',
    metadata = header$text, original_geo_context = 'GSE18508; package summarizeOverlaps source-files section'),
  aligned_records = length(reads), signal_runs = nrow(signal),
  covered_record_bases = as.numeric(sum(sum(coverage))), original_reference_lengths = header$targets,
  semantics = 'Each retained alignment contributes one over M/= /X blocks; excludes D/N, unmapped/secondary/QC-fail/duplicate flags. Supplementary retained; no base-quality filter or fragment/molecule deduplication.',
  bedgraph_coordinate_and_value_roundtrip = TRUE, cram = bridge,
  no_new_reader_or_analysis_engine = TRUE, no_source_characters_or_identity_inferred = TRUE),
  'tools/long-term-audit/v0.27.0-read-coverage-results.json', auto_unbox = TRUE, pretty = TRUE)
cat('Optional read coverage bridge: passed\n')
