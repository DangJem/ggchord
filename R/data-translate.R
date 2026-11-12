#' Translate source-backed CDS with exact codon provenance
#'
#' Uses the complete observed codons in [build_ggchord_cds()], including codons
#' split across segments. Extraction follows biological segment order; negative
#' DNA/RNA strands are complemented. Biostrings supplies genetic codes and
#' translation, and is required only when this function is called.
#'
#' A known genetic_code is required for every model. No sequence, phase,
#' missing end or translation exception is inferred. Stops remain `*`, including
#' internal stops; they are not replaced or removed. The result is an observed
#' translation, not a validation of gene function or a predicted full protein.
#' Partial-end trimming follows the supplied CDS map and is reported in models.
#'
#' @param cds A map from [build_ggchord_cds()].
#' @param sequences Named character vector of complete forward source strings,
#'   or a data frame with unique accver and sequence (e.g. [read_fasta_lengths()]
#'   with keep_sequence = TRUE). DNAStringSet/RNAStringSet are also accepted.
#'   Every CDS source must exist and contain its declared intervals.
#' @param alphabet Explicit `"DNA"` (default) or `"RNA"`; ungapped IUPAC only.
#' @param initiation `"literal"` (default) translates every codon literally.
#'   `"start"` uses the genetic code's alternative initiation codons for the
#'   first codon only; it requires a complete biological 5-prime end.
#' @param ambiguous `"error"` (default) rejects ambiguous codons; `"solve"`
#'   uses Biostrings to resolve them, emitting X when not uniquely resolvable.
#' @return A list with models (coding_sequence, protein_sequence and policies),
#'   nucleotides (one row per observed coding base) and residues (one row per
#'   codon, with source_positions, source_cds_rows and source_cds_segments list
#'   columns in transcription order). Residues use the CDS feature_id as their
#'   observed-protein accver, not an inferred external accession. Nucleotides
#'   retain genomic source_accver/source_position and codon_index/codon_base;
#'   both tables can be displayed with [geom_residue()]. Source input and the
#'   original CDS map remain in provenance. No alignment is computed.
#' @seealso <https://bioconductor.org/packages/Biostrings>
#' @export
#' @examples
#' # Explicit synthetic codon contract, not a biological annotation claim.
#' if (requireNamespace("Biostrings", quietly = TRUE)) {
#'   coding <- data.frame(feature_id = "cds", accver = "source",
#'     start = 1, end = 6, strand = "+", phase = 0)
#'   cds <- build_ggchord_cds(coding, genetic_code = 11)
#'   translate_ggchord_cds(cds, c(source = "ATGTAA"))$residues
#' }
translate_ggchord_cds <- function(cds, sequences, alphabet = c("DNA", "RNA"),
    initiation = c("literal", "start"), ambiguous = c("error", "solve")) {
  caller <- "translate_ggchord_cds()"
  cds <- ggchord_project_model(cds, caller)
  alphabet <- match.arg(alphabet); initiation <- match.arg(initiation); ambiguous <- match.arg(ambiguous)
  if (!requireNamespace("Biostrings", quietly = TRUE))
    ggchord_stop(caller, ": install the optional Biostrings package for translation")
  original <- sequences
  if (inherits(sequences, c("DNAStringSet", "RNAStringSet"))) sequences <- as.character(sequences)
  if (is.data.frame(sequences)) {
    if (identical(ggchord_sequence_unit(sequences, caller), "aa"))
      ggchord_stop(caller, ": source strings require nucleotide bp units")
    ggchord_require_columns(sequences, c("accver", "sequence"), caller)
    ids <- as.character(sequences$accver); strings <- sequences$sequence
  } else { ids <- names(sequences); strings <- sequences }
  if (!is.character(strings) || !length(strings) || anyNA(strings) || is.null(ids) ||
      length(ids) != length(strings) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    ggchord_stop(caller, ": supply forward source strings with unique nonempty accession IDs")
  allowed <- if (alphabet == "DNA") "^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+$" else "^[ACGURYSWKMBDHVNacguryswkmbdhvn]+$"
  if (any(!grepl(allowed, strings))) ggchord_stop(caller, ": source strings must be ungapped IUPAC ", alphabet)
  models <- cds$models
  if (anyNA(models$genetic_code)) ggchord_stop(caller, ": every CDS needs an explicit genetic_code")
  if (initiation == "start" && any(models$partial5))
    ggchord_stop(caller, ": initiation = 'start' requires a complete 5-prime end; use 'literal' for partial CDS")
  nt <- aa <- vector("list", nrow(models)); coding <- protein <- character(nrow(models))
  for (i in seq_len(nrow(models))) {
    m <- models[i, , drop = FALSE]; s <- cds$segments[[m$feature_id]]
    at <- match(m$accver, ids)
    if (is.na(at)) ggchord_stop(caller, ": missing source sequence for '", m$accver, "'")
    # Validate even bases excluded by explicit partial trimming.
    source <- cds$expanded_data[cds$expanded_data$feature_id == m$feature_id, , drop = FALSE]
    if (max(source$end) > nchar(strings[at])) ggchord_stop(caller, ": CDS exceeds source string for '", m$feature_id, "'")
    positions <- unlist(lapply(seq_len(nrow(s)), function(j)
      if (m$strand == "+") seq.int(s$start[j], s$end[j]) else seq.int(s$end[j], s$start[j])), use.names = FALSE)
    symbols <- substring(strings[at], positions, positions)
    if (m$strand == "-") symbols <- ggchord_alignment_complement(symbols, alphabet)
    coding[i] <- paste0(symbols, collapse = "")
    code <- tryCatch(Biostrings::getGeneticCode(as.character(m$genetic_code)),
      error = function(e) ggchord_stop(caller, ": unsupported genetic_code ", m$genetic_code))
    input <- if (alphabet == "DNA") Biostrings::DNAString(coding[i]) else Biostrings::RNAString(coding[i])
    protein[i] <- tryCatch(as.character(Biostrings::translate(input, genetic.code = code,
      no.init.codon = initiation == "literal", if.fuzzy.codon = ambiguous)),
      error = function(e) ggchord_stop(caller, ": translation failed for CDS '", m$feature_id,
        "' (ambiguous = '", ambiguous, "'): ", conditionMessage(e)))
    count <- length(positions); codon <- rep(seq_len(count / 3), each = 3)
    sr <- rep(s$source_row, s$cds_end - s$cds_start + 1)
    ss <- rep(s$source_segment_index, s$cds_end - s$cds_start + 1)
    nt[[i]] <- data.frame(accver = m$feature_id, alignment_column = seq_len(count),
      source_accver = m$accver, source_position = positions, symbol = symbols,
      is_gap = FALSE, strand = m$strand, alphabet = alphabet, sequence_unit = "bp",
      source_verified = TRUE, codon_index = codon, codon_base = rep(1:3, count / 3),
      source_cds_row = sr, source_cds_segment = ss, stringsAsFactors = FALSE)
    residues <- strsplit(protein[i], "", fixed = TRUE)[[1L]]
    aa[[i]] <- data.frame(accver = m$feature_id, alignment_column = seq_along(residues),
      source_accver = m$feature_id, source_position = seq_along(residues),
      symbol = residues, is_gap = FALSE, strand = "+", alphabet = "AA", sequence_unit = "aa",
      source_verified = TRUE, genetic_code = m$genetic_code,
      codon = substring(coding[i], seq(1, count, 3), seq(3, count, 3)),
      is_stop = residues == "*", stringsAsFactors = FALSE)
    groups <- split(seq_len(count), codon)
    aa[[i]]$source_positions <- I(lapply(groups, function(g) positions[g]))
    aa[[i]]$source_cds_rows <- I(lapply(groups, function(g) sr[g]))
    aa[[i]]$source_cds_segments <- I(lapply(groups, function(g) ss[g]))
    aa[[i]]$genomic_accver <- m$accver; aa[[i]]$genomic_strand <- m$strand
  }
  models$coding_sequence <- coding; models$protein_sequence <- protein
  models$initiation <- initiation; models$ambiguous <- ambiguous
  models$translation_status <- ifelse(models$partial5 | models$partial3, "observed_partial", "observed_complete")
  result <- list(models = models, nucleotides = do.call(rbind, nt), residues = do.call(rbind, aa),
    provenance = list(cds = cds, source_sequences = original, alphabet = alphabet,
      engine = "Biostrings", engine_version = as.character(utils::packageVersion("Biostrings")),
      initiation = initiation, ambiguous = ambiguous, coordinate_system = "1-based-inclusive"))
  rownames(result$nucleotides) <- rownames(result$residues) <- NULL
  result
}
