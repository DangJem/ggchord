# Source characters and alignment columns are distinct coordinate spaces.
ggchord_alignment_complement <- function(x, alphabet) {
  if (alphabet == "AA") return(x)
  from <- if (alphabet == "DNA") "ACGTRYSWKMBDHVN" else "ACGURYSWKMBDHVN"
  to <- if (alphabet == "DNA") "TGCAYRSWMKVHDBN" else "UGCAYRSWMKVHDBN"
  chartr(paste0(from, tolower(from)), paste0(to, tolower(to)), x)
}

#' Adapt aligned characters with exact source positions
#'
#' Accepts an existing alignment, never computes one. Equal character counts
#' are necessary but do not establish homology. Each input character has an
#' `alignment_column`; `-` and `.` are retained gaps with missing
#' `source_position`. Case and ambiguous symbols are retained. Lowercase
#' characters consume positions: this is not an A3M insertion parser.
#'
#' Without `source_data`, positions count from one in each ungapped input
#' fragment on the forward strand. They are not complete-genome coordinates.
#' Optional source metadata declares the first nongap position and direction.
#' Negative-strand positions decrease. Already aligned strings must be in
#' alignment orientation: they are never automatically reversed. If a full
#' source string is supplied, characters are verified against that string;
#' DNA/RNA use the complement on the negative strand, AA only reverses order.
#' Missing full source strings mean unverified characters, not invented data.
#'
#' @param sequences Named character vector of aligned strings, or a data frame
#'   with unique `accver` and character `sequence`. All strings must have the
#'   same positive width and each record must contain a nongap symbol.
#' @param alphabet Explicit `"DNA"`, `"RNA"` or `"AA"`. IUPAC nucleotide
#'   ambiguity codes and protein ambiguity/rare-residue codes are accepted;
#'   `*` is accepted only for AA. The alphabet is never guessed.
#' @param source_data Optional data frame with one row per alignment `accver`,
#'   `start` (first nongap source coordinate) and `strand` (`"+"` or `"-"`).
#'   Optional `source_accver` declares a different source ID; optional character
#'   `sequence` supplies the full forward source for verification. Coordinates
#'   must remain positive exact integers. Input metadata is retained in the
#'   `ggchord_conversion` attribute.
#' @return An ordinary long data frame containing `accver`, `alignment_column`,
#'   `source_accver`, `source_position`, `symbol`, `is_gap`, `strand`, `alphabet`,
#'   `sequence_unit`, `source_verified` and `.source_row`. Conversion provenance
#'   retains the original input and normalized source metadata. No identity,
#'   consensus, translation, phase or full sequence length is inferred.
#'   source_verified records whether the record's nongap characters were
#'   checked, including on its gap rows; it is not independent biological
#'   confirmation of the supplied source.
#' @export
#' @examples
#' # Explicit synthetic format example, not a biological alignment claim.
#' aligned <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), "DNA")
#' aligned
as_ggchord_alignment <- function(sequences, alphabet, source_data = NULL) {
  caller <- "as_ggchord_alignment()"
  alphabet <- match.arg(alphabet, c("DNA", "RNA", "AA"))
  if (is.data.frame(sequences)) {
    ggchord_require_columns(sequences, c("accver", "sequence"), caller)
    ids <- as.character(sequences$accver); text <- sequences$sequence
  } else {
    ids <- names(sequences); text <- sequences
  }
  if (!is.character(text) || !length(text) || anyNA(text) || is.null(ids) ||
      length(ids) != length(text) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    ggchord_stop(caller, ": supply character strings with unique nonempty accession IDs")
  widths <- nchar(text, type = "chars")
  if (any(widths < 1) || length(unique(widths)) != 1L)
    ggchord_stop(caller, ": aligned strings must have the same positive width")
  symbols <- lapply(text, strsplit, split = "", fixed = TRUE)
  symbols <- lapply(symbols, `[[`, 1L)
  allowed <- switch(alphabet, DNA = "ACGTRYSWKMBDHVN-.", RNA = "ACGURYSWKMBDHVN-.",
    AA = "ACDEFGHIKLMNPQRSTVWYBXZJUO*-.")
  allowed <- strsplit(allowed, "", fixed = TRUE)[[1L]]
  if (any(!toupper(unlist(symbols, use.names = FALSE)) %in% allowed))
    ggchord_stop(caller, ": character outside the declared alphabet (whitespace is not a residue)")
  metadata <- source_data
  if (is.null(metadata)) metadata <- data.frame(accver = ids, start = 1, strand = "+")
  ggchord_require_columns(metadata, c("accver", "start", "strand"), caller)
  if (nrow(metadata) != length(ids) || anyNA(metadata$accver) ||
      anyDuplicated(as.character(metadata$accver)) || !setequal(as.character(metadata$accver), ids))
    ggchord_stop(caller, ": source_data must have exactly one matching row per accession")
  metadata <- as.data.frame(metadata[match(ids, as.character(metadata$accver)), , drop = FALSE])
  if (!is.numeric(metadata$start) || any(!is.finite(metadata$start)) ||
      any(metadata$start < 1 | metadata$start != floor(metadata$start) | metadata$start > 2^53 - 1) ||
      anyNA(metadata$strand) || any(!metadata$strand %in% c("+", "-")))
    ggchord_stop(caller, ": source starts must be positive integral coordinates and strands must be + or -")
  if (!"source_accver" %in% names(metadata)) metadata$source_accver <- ids
  if (anyNA(metadata$source_accver) || any(!nzchar(as.character(metadata$source_accver))))
    ggchord_stop(caller, ": source accession IDs must be nonempty")
  if ("sequence" %in% names(metadata) &&
      (!is.character(metadata$sequence) || anyNA(metadata$sequence)))
    ggchord_stop(caller, ": full source strings must be nonmissing character values")
  pieces <- lapply(seq_along(ids), function(i) {
    value <- symbols[[i]]; gap <- value %in% c("-", ".")
    if (all(gap)) ggchord_stop(caller, ": record contains only gaps: ", ids[i])
    position <- metadata$start[i] + (cumsum(!gap) - 1) * if (metadata$strand[i] == "+") 1 else -1
    position[gap] <- NA_real_
    if (any(position[!gap] < 1 | position[!gap] > 2^53 - 1))
      ggchord_stop(caller, ": source interval lies outside positive exact coordinates: ", ids[i])
    verified <- "sequence" %in% names(metadata)
    if (verified) {
      source <- metadata$sequence[i]
      if (max(position, na.rm = TRUE) > nchar(source))
        ggchord_stop(caller, ": source interval exceeds full source string: ", ids[i])
      observed <- substring(source, position[!gap], position[!gap])
      if (metadata$strand[i] == "-") observed <- ggchord_alignment_complement(observed, alphabet)
      if (any(toupper(observed) != toupper(value[!gap])))
        ggchord_stop(caller, ": aligned characters disagree with the declared full source: ", ids[i])
    }
    data.frame(accver = ids[i], alignment_column = seq_along(value),
      source_accver = as.character(metadata$source_accver[i]), source_position = position,
      symbol = value, is_gap = gap, strand = as.character(metadata$strand[i]),
      alphabet = alphabet, sequence_unit = if (alphabet == "AA") "aa" else "bp",
      source_verified = verified, .source_row = i, stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, pieces); rownames(result) <- NULL
  attr(result, "ggchord_conversion") <- list(format = "aligned characters", alphabet = alphabet,
    source_data = metadata, source_input = sequences,
    coordinates = "one-based alignment columns and nullable one-based source positions")
  result
}

#' Read aligned FASTA with source-to-column mapping
#'
#' Reads one alignment file. Whitespace wrapping is removed by the FASTA reader;
#' gaps and case remain unchanged. Records must have equal width. Does not run
#' an aligner or infer strand, genomic origins, mismatch classes or CDS phase.
#' @param file A single local aligned FASTA path; no multi-alignment concatenation.
#' @inheritParams as_ggchord_alignment
#' @param source_file Logical. Also retain `.source_file` on each character row.
#' @return The long frame from [as_ggchord_alignment()], with the physical FASTA
#'   header `.source_line`, optional `.source_file`, and import provenance.
#' @seealso [geom_residue()], [map_alignment_position()]
#' @export
#' @examples
#' \dontrun{
#' # Requires the source checkout or tutorial input bundle.
#' dna <- read_aligned_fasta("examples/alignment/opuntia-aligned.fasta", "DNA")
#' }
read_aligned_fasta <- function(file, alphabet, source_data = NULL, source_file = FALSE) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !file.exists(file) || dir.exists(file))
    ggchord_stop("read_aligned_fasta(): supply a single existing file")
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop("read_aligned_fasta(): source_file must be TRUE or FALSE")
  sequences <- read_fasta_lengths(file, keep_sequence = TRUE, source_file = source_file)
  result <- as_ggchord_alignment(sequences, alphabet, source_data)
  result$.source_line <- sequences$.source_line[result$.source_row]
  if (source_file) result$.source_file <- sequences$.source_file[result$.source_row]
  result <- ggchord_format_import(result, file, "aligned FASTA", nrow(result),
    "aligned characters; no genomic coordinate conversion")
  provenance <- attr(result, "ggchord_import")
  provenance$output_coordinates <- "1-based alignment columns and nullable 1-based source positions"
  attr(result, "ggchord_import") <- provenance
  result
}

#' Map exact alignment columns and source positions
#'
#' Looks up existing cells within one alignment row. Never interpolates across
#' gaps or discarded columns. Filter an alignment with ordinary data operations
#' before looking up positions to restrict the available window.
#' @param data Long alignment frame from [as_ggchord_alignment()].
#' @param accver One alignment-row ID, even if its `source_accver` differs.
#' @param position Numeric positions, positive integral or NA.
#' @param inverse FALSE maps alignment columns to source positions; TRUE maps
#'   source positions to alignment columns.
#' @return A frame with input_position, accver, alignment_column, source_accver,
#'   source_position, is_gap and matched. A known gap column has matched TRUE
#'   and source_position NA; an unavailable position has matched FALSE.
#' @export
#' @examples
#' # Explicit synthetic format example.
#' aligned <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), "DNA")
#' map_alignment_position(aligned, "A", c(2, 3))
#' map_alignment_position(aligned, "A", 3, inverse = TRUE)
map_alignment_position <- function(data, accver, position, inverse = FALSE) {
  caller <- "map_alignment_position()"
  ggchord_require_columns(data, c("accver", "alignment_column", "source_accver", "source_position", "is_gap"), caller)
  if (!is.character(accver) || length(accver) != 1L || is.na(accver) || !nzchar(accver) ||
      !is.logical(inverse) || length(inverse) != 1L || is.na(inverse))
    ggchord_stop(caller, ": select one nonempty accession and logical inverse")
  if (!is.numeric(position) || any(!is.na(position) &
      (!is.finite(position) | position < 1 | position != floor(position) | position > 2^53 - 1)))
    ggchord_stop(caller, ": positions must be positive exact integers or NA")
  part <- data[!is.na(data$accver) & data$accver == accver, , drop = FALSE]
  if (anyDuplicated(part$alignment_column) || anyDuplicated(part$source_position[!is.na(part$source_position)]))
    ggchord_stop(caller, ": ambiguous repeated columns or source positions within this accession")
  key <- if (inverse) part$source_position else part$alignment_column
  at <- match(position, key); at[is.na(position)] <- NA_integer_
  # Inverse NA must not match one of the gap cells.
  data.frame(input_position = position, accver = rep(accver, length(position)),
    alignment_column = part$alignment_column[at], source_accver = part$source_accver[at],
    source_position = part$source_position[at], is_gap = part$is_gap[at], matched = !is.na(at))
}
