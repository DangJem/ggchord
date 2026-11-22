# Source characters and alignment columns are distinct coordinate spaces.
ggchord_alignment_complement <- function(x, alphabet) {
  if (alphabet == "AA") return(x)
  from <- if (alphabet == "DNA") "ACGTRYSWKMBDHVN" else "ACGURYSWKMBDHVN"
  to <- if (alphabet == "DNA") "TGCAYRSWMKVHDBN" else "UGCAYRSWMKVHDBN"
  chartr(paste0(from, tolower(from)), paste0(to, tolower(to)), x)
}

ggchord_alignment_invalid_length <- function(value) {
  if (!is.numeric(value)) return(!is.na(value))
  is.nan(value) | (!is.na(value) &
    (!is.finite(value) | value < 1 | value > 2^53 - 1 | value != floor(value)))
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
#'   `sequence` supplies the full forward source for verification. Optional
#'   `source_length` gives the full source length as a positive exact integer.
#'   Linear nongap positions must be within a declared length. If both length
#'   and full sequence are supplied, they must agree. Input metadata is retained
#'   in the `ggchord_conversion` attribute.
#' @param circular Explicit logical declaration that every source is circular,
#'   default FALSE. Requires source_data with full source_length or full sequence.
#'   Coordinates wrap at that length in the declared strand direction; nongap
#'   characters may cover at most one complete tour. Adds source_circular and
#'   source_length columns. A circular plot never activates this conversion.
#' @return An ordinary long data frame containing `accver`, `alignment_column`,
#'   `source_accver`, `source_position`, `symbol`, `is_gap`, `strand`, `alphabet`,
#'   `sequence_unit`, `source_verified` and `.source_row`. Conversion provenance
#'   retains the original input and normalized source metadata. No identity,
#'   consensus, translation, phase or full sequence length is inferred.
#'   source_verified records whether the record's nongap characters were
#'   checked, including on its gap rows; it is not independent biological
#'   confirmation of the supplied source.
#'   `source_scope` is `"source"` with explicit source metadata and
#'   `"fragment"` for positions counted within each input fragment.
#'   A declared `source_length` is retained on the corresponding cells.
#' @export
#' @examples
#' library(ggplot2)
#'
#' cells <- as_ggchord_alignment(c(A = "AC-G", B = "ACTG"), alphabet = "DNA")
#' ggplot(cells, aes(alignment_column, accver, label = symbol)) +
#'   geom_residue(aes(fill = symbol))
as_ggchord_alignment <- function(sequences, alphabet, source_data = NULL, circular = FALSE) {
  caller <- "as_ggchord_alignment()"
  alphabet <- ggchord_match_arg(alphabet, c("DNA", "RNA", "AA"))
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
  if (!is.logical(circular) || length(circular) != 1L || is.na(circular))
    ggchord_stop(caller, ": circular must be TRUE or FALSE")
  if (circular && is.null(source_data))
    ggchord_stop(caller, ": circular coordinates require explicit source_data")
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
  declared_length <- metadata[["source_length"]]
  source_length <- rep(NA_real_, nrow(metadata))
  if (!is.null(declared_length)) {
    invalid <- ggchord_alignment_invalid_length(declared_length)
    if (any(invalid)) ggchord_stop(caller, ": source_length must be a positive exact integer for ",
      paste(ids[invalid], collapse = ", "))
    source_length <- as.numeric(declared_length)
  }
  if ("sequence" %in% names(metadata)) {
    actual_length <- nchar(metadata$sequence, type = "chars")
    conflict <- !is.na(source_length) & source_length != actual_length
    if (any(conflict)) ggchord_stop(caller, ": source_length disagrees with the full source string for ",
      paste(ids[conflict], collapse = ", "))
    if (circular) source_length <- actual_length
  }
  if (circular && (anyNA(source_length) || any(source_length < 1 | metadata$start > source_length)))
    ggchord_stop(caller, ": circular source_data requires exact full source_length (or full sequence) and start within it")
  if (circular) metadata$source_length <- source_length
  pieces <- lapply(seq_along(ids), function(i) {
    value <- symbols[[i]]; gap <- value %in% c("-", ".")
    if (all(gap)) ggchord_stop(caller, ": record contains only gaps: ", ids[i])
    position <- metadata$start[i] + (cumsum(!gap) - 1) * if (metadata$strand[i] == "+") 1 else -1
    if (circular) {
      n <- source_length[i]; start <- metadata$start[i]; k <- cumsum(!gap) - 1
      if (sum(!gap) > n)
        ggchord_stop(caller, ": circular record may not revisit a source position: ", ids[i])
      # Avoid an intermediate start + k above the exact-integer range.
      position <- if (metadata$strand[i] == "+") {
        ifelse(k <= n - start, start + pmin(k, n - start), k - (n - start))
      } else ifelse(k < start, start - pmin(k, start - 1), n - (k - start))
    }
    position[gap] <- NA_real_
    if (any(position[!gap] < 1 | position[!gap] > 2^53 - 1))
      ggchord_stop(caller, ": source interval lies outside positive exact coordinates: ", ids[i])
    if (!is.na(source_length[i]) && any(position[!gap] > source_length[i]))
      ggchord_stop(caller, ": source interval exceeds source_length (", source_length[i], ") for ", ids[i])
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
  result$source_scope <- if (is.null(source_data)) "fragment" else "source"
  if (circular) result$source_circular <- TRUE
  if (circular || !is.null(declared_length)) result$source_length <- source_length[result$.source_row]
  attr(result, "ggchord_conversion") <- list(format = "aligned characters", alphabet = alphabet,
    source_data = metadata, source_input = sequences,
    coordinates = "one-based alignment columns and nullable one-based source positions")
  result
}

#' Read an existing character alignment from FASTA
#'
#' All records must have an equal positive width. Gaps occupy alignment
#' columns without receiving source residue positions. Without explicit
#' source mapping, reported source positions refer to local input fragments.
#'
#' @param file One local aligned FASTA path.
#' @param alphabet Explicit DNA, RNA, or AA.
#' @param source_data Optional source-coordinate metadata as in
#'   [as_ggchord_alignment()].
#' @param circular Whether the explicit source positions can cross a known
#'   circular origin.
#' @param source_file Whether to include a row-level file path.
#' @return A long alignment data frame preserving character and source columns.
#' @export
#' @examples
#' f <- tempfile(fileext = ".fa")
#' writeLines(c(">A", "AC-G", ">B", "ACTG"), f)
#' cells <- read_aligned_fasta(f, alphabet = "DNA")
#' unlink(f)
read_aligned_fasta <- function(file, alphabet, source_data = NULL, circular = FALSE, source_file = FALSE) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !file.exists(file) || dir.exists(file))
    ggchord_stop("read_aligned_fasta(): supply a single existing file")
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop("read_aligned_fasta(): source_file must be TRUE or FALSE")
  sequences <- read_fasta_lengths_single(file, keep_sequence = TRUE, allow_gaps = TRUE)
  result <- as_ggchord_alignment(sequences, alphabet, source_data, circular = circular)
  result$.source_line <- sequences$.source_line[result$.source_row]
  if (source_file) result$.source_file <- rep(file, nrow(result))
  result <- ggchord_format_import(result, file, "aligned FASTA", nrow(result),
    "aligned characters; no genomic coordinate conversion")
  provenance <- attr(result, "ggchord_import")
  provenance$output_coordinates <- "1-based alignment columns and nullable 1-based source positions"
  attr(result, "ggchord_import") <- provenance
  result
}

#' Look up exact source positions and alignment columns
#' @param alignment A character table or a registered character plot.
#' @param positions Query table with source_accver/source_position or
#'   accver/alignment_column, according to from.
#' @param from,to The input and output spaces, source or alignment.
#' @param instance Optional accession of a specific alignment row.
#' @return A correspondence table retaining all matching rows, query IDs and
#'   explicit reasons for gaps and missing positions. For focused data it
#'   includes the original `alignment_row`, `locus_id` and `local_position`.
#' @export
#' @examples
#' cells <- as_ggchord_alignment(c(A = "AC-G"), alphabet = "DNA")
#' map_alignment_position(cells,
#'   data.frame(accver = "A", alignment_column = 3),
#'   from = "alignment", to = "source")
map_alignment_position <- function(alignment, positions, from = "source",
    to = "alignment", instance = NULL) {
  caller <- "map_alignment_position()"
  spaces <- c("source", "alignment")
  from <- ggchord_match_arg(from, spaces)
  to <- ggchord_match_arg(to, spaces)
  data <- alignment
  if (inherits(alignment, "ggplot")) {
    data <- alignment$ggchord_alignment
    if (is.null(data)) ggchord_stop(caller, ": register the character plot first")
    ggchord_alignment_layout(ggplot_build(alignment), data)
  }
  required <- c("accver", "alignment_column", "source_accver", "source_position", "is_gap")
  ggchord_require_columns(data, required, caller)
  keys <- if (from == "source") c("source_accver", "source_position") else c("accver", "alignment_column")
  ggchord_require_columns(positions, keys, caller)
  coords <- positions[[keys[2L]]]
  if (!is.numeric(coords) || any(!is.na(coords) &
      (!is.finite(coords) | coords < 1 | coords != floor(coords) | coords > 2^53 - 1)))
    ggchord_stop(caller, ": positions must be positive exact integers or NA")
  if (!is.null(instance)) {
    if (!is.character(instance) || length(instance) != 1L || is.na(instance) ||
        !instance %in% data$accver) ggchord_stop(caller, ": unknown alignment row")
    data <- data[data$accver == instance, , drop = FALSE]
  }
  if (anyDuplicated(data[c("accver", "alignment_column")]))
    ggchord_stop(caller, ": select a single alignment block")
  fields <- intersect(c(required, "strand", "sequence_unit", "symbol", "alphabet",
    "alignment_row", "locus_id", "local_position"), names(data))
  out <- lapply(seq_len(nrow(positions)), function(i) {
    at <- which(data[[keys[1L]]] == positions[[keys[1L]]][i] & data[[keys[2L]]] == coords[i])
    if (!length(at)) {
      row <- data[NA_integer_, fields, drop = FALSE]
      for (field in keys) row[[field]] <- positions[[field]][i]
      row$status <- "unavailable"
      row$reason <- "not_present"
    } else {
      row <- data[at, fields, drop = FALSE]
      gap <- to == "source" & row$is_gap
      row$status <- ifelse(gap, "unavailable", "available")
      row$reason <- ifelse(gap, "gap", NA_character_)
    }
    row$input_row <- i
    row
  })
  if (length(out)) return(do.call(rbind, out))
  empty <- data[FALSE, fields, drop = FALSE]
  empty$status <- empty$reason <- character()
  empty$input_row <- integer()
  empty
}
