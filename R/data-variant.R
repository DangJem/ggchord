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

#' @noRd
ggchord_variant_dictionary <- function(tables) {
  x <- ggchord_rbind_fill(tables)
  ids <- unique(x$accver)
  lengths <- vapply(ids, function(id) {
    values <- unique(x$length[x$accver == id & !is.na(x$length)])
    if (length(values) > 1L) ggchord_stop("read_vcf(): conflicting full length for ", id)
    if (length(values)) values else NA_real_
  }, numeric(1))
  data.frame(accver = ids %||% character(), length = unname(lengths),
    sequence_unit = rep("bp", length(ids)))
}
#' Read local variant and structural relationship records
#'
#' Use optional VariantAnnotation to decode variants, and
#' StructuralVariantAnnotation for requested structural relationships.
#' Preserve multiallelic records. A missing length in the source dictionary
#' remains unknown until a complete length is supplied explicitly.
#'
#' @param file Local VCF or BCF paths.
#' @param structural Decode supported structural breakend relationships.
#' @param python Python executable with pysam for BCF conversion.
#' @param source_file Whether to include row-level source paths.
#' @return A ggchord_import collection with a source dictionary,
#'   features$variant, optional links$structural, and a report with
#'   breakends and unplaced records.
#' @export
#' @examples
#' if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
#'   f <- tempfile(fileext = ".vcf")
#'   writeLines(c("##fileformat=VCFv4.3", "##contig=<ID=A>",
#'     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
#'     "A\t9\t.\tA\tG\t.\t.\t."), f)
#'   x <- read_vcf(f)
#'   as_ggchord_seq(x$sequences, lengths = c(A = 12))
#'   unlink(f)
#' }
read_vcf <- function(file, structural = FALSE,
                     python = 'python3', source_file = FALSE) {
  ggchord_local_console_errors()
  caller <- 'read_vcf()'
  ggchord_conversion_require('VariantAnnotation', caller)
  for (name in c('structural', 'source_file')) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) ggchord_stop(caller, ': ', name, ' must be TRUE or FALSE')
  }
  paths <- ggchord_resolve_files(file, NULL, 'read_vcf')
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
  sequences <- ggchord_variant_dictionary(lapply(records, `[[`, "sequence_info"))
  ggchord_import_collection(sequences, features = list(variant = ggchord_rbind_fill(features)),
    links = if (structural) list(structural = ggchord_rbind_fill(links) %||% empty) else list(),
    report = list(records = records, structural_requested = structural,
      breakends = ggchord_rbind_fill(breakends) %||% empty,
      decoders = list(VariantAnnotation = as.character(utils::packageVersion('VariantAnnotation')),
        StructuralVariantAnnotation = if (structural) as.character(utils::packageVersion('StructuralVariantAnnotation')) else NULL),
      coordinates = '1-based-inclusive reference allele spans and decoded breakpoint uncertainty intervals; no additional shift'))
}
