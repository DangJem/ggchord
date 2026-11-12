# Mature parsers own VCF/BCF and structural notation; this adapter owns tables.
# https://samtools.github.io/hts-specs/VCFv4.5.pdf
# https://bioconductor.org/packages/StructuralVariantAnnotation/
ggchord_variant_bcf_script <- function() r"---(
import json, sys
try:
    import pysam
except ImportError:
    raise SystemExit('The selected Python needs optional pysam for modern BCF. Select an existing environment or supply VCF text instead.')
with open(sys.argv[1], encoding='utf-8') as handle:
    config = json.load(handle)
with pysam.VariantFile(config['source']) as source:
    with pysam.VariantFile(config['destination'], 'w', header=source.header) as target:
        for record in source:
            target.write(record)
with open(sys.argv[2], 'w', encoding='utf-8') as handle:
    json.dump({'package':'pysam', 'version':pysam.__version__,
               'htslib_version':pysam.version.__htslib_version__}, handle)
)---"

ggchord_variant_bcf_text <- function(path, destination, python, caller) {
  result <- ggchord_parse_external_source(ggchord_variant_bcf_script(),
    list(source = normalizePath(path), destination = destination), python, caller)
  decoder <- result$data
  decoder$messages <- result$messages
  decoder
}

ggchord_variant_breakends <- function(vcf, caller) {
  ggchord_conversion_require('StructuralVariantAnnotation', caller)
  if (anyDuplicated(names(vcf))) ggchord_stop(caller,
    ': structural decoding needs unambiguous unique record IDs; duplicated IDs remain available with structural = FALSE')
  warnings <- character()
  capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
  })
  # Never infer absent mates or collapse positional uncertainty to a point.
  paired <- capture(StructuralVariantAnnotation::breakpointRanges(vcf,
    inferMissingBreakends = FALSE, nominalPosition = FALSE))
  single <- capture(StructuralVariantAnnotation::breakendRanges(vcf, nominalPosition = FALSE))
  ranges <- capture(c(paired, single))
  metadata <- as.data.frame(S4Vectors::mcols(ranges))
  ends <- data.frame(accver = as.character(GenomicRanges::seqnames(ranges)),
    start = as.numeric(GenomicRanges::start(ranges)), end = as.numeric(GenomicRanges::end(ranges)),
    breakend_orientation = as.character(GenomicRanges::strand(ranges)),
    decoder_id = names(ranges), stringsAsFactors = FALSE)
  ends$.source_record <- match(metadata$sourceId, names(vcf))
  if (anyNA(ends$.source_record)) ggchord_stop(caller,
    ': structural decoder source IDs cannot be matched to original records')
  ends$decoder_metadata <- I(lapply(seq_len(nrow(ends)), function(i) metadata[i, , drop = FALSE]))
  partner <- match(metadata$partner, ends$decoder_id)
  keep <- which(!is.na(partner) & seq_len(nrow(ends)) < partner)
  # Refuse inconsistent partner graphs instead of inventing/deduplicating events.
  if (length(keep) && any(is.na(partner[partner[keep]]) | partner[partner[keep]] != keep))
    ggchord_stop(caller, ': structural decoder returned nonreciprocal partners')
  links <- data.frame(qaccver = ends$accver[keep], saccver = ends$accver[partner[keep]],
    qstart = ends$start[keep], qend = ends$end[keep],
    sstart = ends$start[partner[keep]], send = ends$end[partner[keep]],
    q_breakend_orientation = ends$breakend_orientation[keep],
    s_breakend_orientation = ends$breakend_orientation[partner[keep]],
    q_decoder_id = ends$decoder_id[keep], s_decoder_id = ends$decoder_id[partner[keep]],
    q_source_record = ends$.source_record[keep], s_source_record = ends$.source_record[partner[keep]],
    stringsAsFactors = FALSE)
  information <- ggchord_conversion_seqinfo(GenomicRanges::seqinfo(vcf), caller)
  full_length <- information$length[match(ends$accver, information$accver)]
  ends$.plottable <- ends$start >= 1 & ends$end >= ends$start &
    (is.na(full_length) | ends$end <= full_length)
  placed <- ends$.plottable[keep] & ends$.plottable[partner[keep]]
  unplaced <- links[!placed, , drop = FALSE]
  links <- links[placed, , drop = FALSE]
  # Breakend orientation describes derivative sequence, not homology strand.
  # Keep interval endpoints ascending and do not manufacture pident/Orientation.
  list(breakends = ends, links = links, unplaced_links = unplaced,
    parsed_ranges = ranges, warnings = unique(warnings))
}

#' Read VCF or BCF with optional mature variant decoders
#'
#' VariantAnnotation reads local VCF text (including gzip). Files ending in
#' .bcf or .bcf.gz first use pysam/HTSlib in the explicitly selected Python
#' environment to transcode to temporary VCF. No parser is installed or run at
#' plot time. This is an in-memory adapter for selected local files, not an
#' indexed chromosome-scale query or variant caller.
#'
#' Each placed source record remains one row, including multiallelic records. Variant
#' features show the reference allele span at POS, including its padding base;
#' they are not claimed to be the inserted sequence or a symbolic SV footprint.
#' Original IDs, REF, ALT lists, QUAL, FILTER, INFO and sample/genotype arrays
#' are retained. Missing IDs stay missing. No alleles or samples are combined,
#' filtered, phased or normalised. No support, Identity or full lengths are
#' inferred from positions. Header length entries that are missing remain NA.
#' Telomeric or out-of-bounds reference spans remain in the complete parsed
#' source and `report$records[[i]]$unplaced_features` rather than being clipped
#' into invented reference features. Empty VCFs return empty feature tables.
#'
#' Explicit structural decoding uses optional StructuralVariantAnnotation,
#' with missing-mate inference disabled and positional uncertainty retained.
#' Each reciprocal decoded pair produces one link, with both source record
#' indices; single and unpaired records remain in the source/report. Breakend
#' orientation describes derivative sequence, not homology or feature strand.
#' Links retain ascending confidence intervals and separate orientation fields;
#' no pident, relative Orientation or connection to an absent mate is invented.
#' Decoder-derived metadata is labelled separately from original evidence.
#' Breakends outside known bounds have .plottable = FALSE; affected pairs stay
#' in the report's unplaced_links instead of being clipped to exact boundaries.
#'
#' @param file Optional local path vector or glob patterns, before files.
#' @param files Optional additional VCF/BCF paths or patterns.
#' @param structural Logical. Explicitly request optional structural decoding,
#'   default FALSE. Duplicate record IDs fail in this mode because mates would
#'   be ambiguous. All source records remain accessible in the default mode.
#' @param python Existing interpreter with optional pysam for BCF only.
#'   Ignored for VCF text; default "python3". No environment is modified.
#' @param source_file Add .source_file to returned table rows, default FALSE.
#' @return A list with features$variant (reference allele spans), breakends,
#'   links and report. The report retains each complete parsed VCF object,
#'   source record text, original headers/IDs, decoder versions, warnings and
#'   sequence dictionaries with unknown lengths preserved. Original INFO and
#'   genotype evidence stays in report$records; no caller-specific support
#'   field is silently promoted to a universal count. Table .source_input and
#'   .source_record refer to the report file and its record order, never an
#'   invented biological observation ID. BCF text is decoder output rather
#'   than original binary bytes; no physical BCF line numbers are assigned.
#' @seealso [geom_feature()], [geom_link_ribbon()], [read_bedpe()]
#' @export
#' @examples
#' if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
#'   path <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
#'   variants <- read_vcf(path)
#'   head(variants$features$variant[c("accver", "start", "end", "feature_id")])
#' }
read_vcf <- function(file = NULL, files = NULL, structural = FALSE,
                     python = 'python3', source_file = FALSE) {
  caller <- 'read_vcf()'
  ggchord_conversion_require('VariantAnnotation', caller)
  for (name in c('structural', 'source_file')) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) ggchord_stop(caller, ': ', name, ' must be TRUE or FALSE')
  }
  paths <- ggchord_resolve_files(file, files, 'read_vcf')
  temporary <- tempfile('ggchord-variants-'); dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  records <- list(); features <- list(); breakends <- list(); links <- list()
  for (i in seq_along(paths)) {
    path <- paths[i]; binary <- grepl('[.]bcf([.]gz)?$', path, ignore.case = TRUE)
    decoded <- if (binary) file.path(temporary, paste0(i, '.vcf')) else path
    bridge <- if (binary) ggchord_variant_bcf_text(path, decoded, python, caller) else NULL
    text <- ggchord_read_text(decoded, caller)
    rows <- which(nzchar(text) & !startsWith(text, '#'))
    fields <- lapply(text[rows], strsplit, split = '\t', fixed = TRUE)
    ids <- vapply(fields, function(x) x[[1L]][3L], character(1))
    ids[ids == '.'] <- NA_character_
    warnings <- character()
    vcf <- withCallingHandlers(VariantAnnotation::readVcf(decoded), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')
    })
    if (nrow(vcf) != length(rows)) ggchord_stop(caller, ': decoder record count differs from source text: ', path)
    ranges <- GenomicRanges::granges(vcf)
    out <- data.frame(accver = as.character(GenomicRanges::seqnames(ranges)),
      start = as.numeric(GenomicRanges::start(ranges)),
      end = as.numeric(GenomicRanges::start(ranges)) + nchar(as.character(VariantAnnotation::ref(vcf))) - 1,
      feature_id = ids, strand = rep('.', nrow(vcf)), directionality = rep('none', nrow(vcf)), type = rep('variant', nrow(vcf)),
      REF = as.character(VariantAnnotation::ref(vcf)), QUAL = VariantAnnotation::qual(vcf),
      FILTER = VariantAnnotation::filt(vcf), .source_input = rep(i, nrow(vcf)),
      .source_record = seq_len(nrow(vcf)), stringsAsFactors = FALSE)
    out$ALT <- I(lapply(as.list(VariantAnnotation::alt(vcf)), as.character))
    if (!binary) out$.source_line <- rows
    information <- ggchord_conversion_seqinfo(GenomicRanges::seqinfo(vcf), caller)
    full_length <- information$length[match(out$accver, information$accver)]
    placed <- out$start >= 1 & out$end >= out$start & (is.na(full_length) | out$end <= full_length)
    unplaced <- out[!placed, , drop = FALSE]
    if (nrow(unplaced)) unplaced$reason <- 'reference allele span outside known reference bounds'
    features[[i]] <- ggchord_add_source_file(out[placed, , drop = FALSE], path, source_file)
    sv <- if (structural) ggchord_variant_breakends(vcf, caller) else NULL
    for (name in c('breakends', 'links')) {
      table <- sv[[name]]
      if (is.null(table)) next
      table$.source_input <- rep(i, nrow(table))
      table <- ggchord_add_source_file(table, path, source_file)
      if (name == 'breakends') breakends[[i]] <- table else links[[i]] <- table
    }
    records[[i]] <- list(file = path, parsed = vcf, record_text = text[rows],
      header_text = text[startsWith(text, '#')], original_ids = ids,
      text_origin = if (binary) 'pysam transcoded BCF, not original physical lines' else 'original decompressed VCF',
      binary_bridge = bridge, warnings = unique(c(warnings, sv$warnings)),
      structural_ranges = sv$parsed_ranges,
      unplaced_links = sv$unplaced_links,
      unplaced_features = unplaced, sequence_info = information)
  }
  empty <- data.frame()
  list(features = list(variant = ggchord_rbind_fill(features)),
    breakends = ggchord_rbind_fill(breakends) %||% empty,
    links = ggchord_rbind_fill(links) %||% empty,
    report = list(records = records, structural_requested = structural,
      decoders = list(VariantAnnotation = as.character(utils::packageVersion('VariantAnnotation')),
        StructuralVariantAnnotation = if (structural) as.character(utils::packageVersion('StructuralVariantAnnotation')) else NULL),
      coordinates = '1-based-inclusive reference allele spans and decoded breakpoint uncertainty intervals; no additional shift'))
}
