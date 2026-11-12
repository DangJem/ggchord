# Existing protein annotations only: no domain caller or online lookup.
# https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html
# https://pfam-docs.readthedocs.io/en/latest/api.html

ggchord_signature_integer <- function(value, name, caller) {
  if (length(value) != 1L || !(is.numeric(value) || is.character(value)))
    ggchord_stop(caller, ': ', name, ' must be one positive source coordinate/length')
  number <- suppressWarnings(as.numeric(value))
  if (!is.finite(number) || number < 1 || number != floor(number) || number > 2^53-1)
    ggchord_stop(caller, ': ', name, ' must be an exact positive integer')
  number
}

ggchord_signature_text <- function(value, name, caller, optional = FALSE) {
  if (optional && is.null(value)) return(NA_character_)
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
    ggchord_stop(caller, ': missing/noncharacter ', name)
  value
}

ggchord_signature_row <- function(id, length, start, end, signature, description,
    analysis, type, model, score, fields, match, fragment, caller) {
  id <- ggchord_signature_text(id, 'protein accession', caller)
  signature <- ggchord_signature_text(signature, 'signature accession', caller)
  length <- ggchord_signature_integer(length, 'full protein length', caller)
  start <- ggchord_signature_integer(start, 'start', caller)
  end <- ggchord_signature_integer(end, 'end', caller)
  if (end < start || end > length) ggchord_stop(caller, ': signature interval exceeds its declared protein length')
  out <- data.frame(accver = id, start = start, end = end, strand = '.', directionality = 'none',
    feature_id = NA_character_, anno = signature, type = type,
    signature_id = signature, signature_description = description,
    analysis = analysis, model_id = model, source_length = length,
    .source_match = match, .source_fragment = fragment, stringsAsFactors = FALSE)
  # Scores are opaque native evidence. Neither strings nor JSON numbers imply
  # a shared metric, and repeated model IDs are not unique feature instances.
  out$source_score <- I(list(score))
  out$source_fields <- I(list(fields))
  out
}

ggchord_signature_tsv <- function(text, caller) {
  index <- which(nzchar(trimws(text)) & !grepl('^[[:space:]]*#', text))
  rows <- list(); original <- list(); sequences <- list()
  for (i in seq_along(index)) {
    # Retain even empty trailing fields; metadata is not whitespace-delimited.
    fields <- strsplit(paste0(text[index[i]], '\t.'), '\t', fixed = TRUE)[[1L]]
    fields <- head(fields, -1L)
    if (length(fields) < 11L || length(fields) > 15L) ggchord_stop(caller,
      ': signature TSV needs 11-15 tab-separated fields at line ', index[i])
    description <- if (nzchar(fields[6L])) fields[6L] else NA_character_
    out <- ggchord_signature_row(fields[1L], fields[3L], fields[7L], fields[8L],
      fields[5L], description, fields[4L], 'signature_match', NA_character_,
      fields[9L], fields, i, 1L, caller)
    out$source_md5 <- fields[2L]
    out$.source_line <- index[i]; out$.source_record <- i
    rows[[i]] <- out; original[[i]] <- fields
    sequences[[i]] <- out[c('accver', 'source_length', 'source_md5')]
  }
  list(rows = rows, sequences = sequences, parsed = original, pagination = NULL)
}

ggchord_signature_json <- function(parsed, caller) {
  entries <- parsed$results
  if (is.null(entries) && !is.null(parsed$metadata) && !is.null(parsed$proteins)) entries <- list(parsed)
  if (!is.list(entries)) ggchord_stop(caller, ': expected cached entry JSON with results containing metadata and proteins')
  rows <- list(); sequences <- list(); match <- 0L
  for (i in seq_along(entries)) {
    entry <- entries[[i]]; metadata <- entry$metadata
    if (!is.list(entry$proteins)) ggchord_stop(caller, ': entry JSON needs a proteins list')
    signature <- ggchord_signature_text(metadata$accession, 'signature accession', caller)
    description <- metadata$name
    if (is.list(description)) description <- description$name %||% description$short
    description <- ggchord_signature_text(description, 'signature name', caller, optional = TRUE)
    analysis <- ggchord_signature_text(metadata$source_database, 'source database', caller)
    type <- ggchord_signature_text(metadata$type, 'entry type', caller)
    for (j in seq_along(entry$proteins)) {
      protein <- entry$proteins[[j]]
      id <- ggchord_signature_text(protein$accession, 'protein accession', caller)
      length <- ggchord_signature_integer(protein$protein_length, 'full protein length', caller)
      sequences[[length(sequences)+1L]] <- data.frame(accver = id, source_length = length,
        source_md5 = NA_character_, stringsAsFactors = FALSE)
      for (k in seq_along(protein$entry_protein_locations)) {
        location <- protein$entry_protein_locations[[k]]
        match <- match+1L
        if (!length(location$fragments)) ggchord_stop(caller, ': source location lacks fragments; cannot infer a footprint')
        pieces <- lapply(seq_along(location$fragments), function(n) {
          fragment <- location$fragments[[n]]
          model <- ggchord_signature_text(location$model, 'model ID', caller, optional = TRUE)
          ggchord_signature_row(id, length, fragment$start, fragment$end,
            signature, description, analysis, type, model, location$score, fragment, match, n, caller)
        })
        out <- pieces[[1L]]
        segments <- do.call(rbind, lapply(seq_along(pieces), function(n) data.frame(
          start = pieces[[n]]$start, end = pieces[[n]]$end, strand = '.',
          segment_index = n, segment_type = 'standard',
          fragment_status = ggchord_signature_text(location$fragments[[n]][['dc-status']],
            'fragment status', caller, optional = TRUE), stringsAsFactors = FALSE)))
        out$start <- min(segments$start); out$end <- max(segments$end)
        out$segments <- I(list(segments))
        out$source_fields <- I(list(location$fragments))
        out$.source_fragment <- NULL
        out$source_md5 <- NA_character_
        out$.source_record <- i; out$.source_protein <- j; out$.source_location <- k
        representative <- if (is.null(location$representative)) NA else location$representative
        if (!is.logical(representative) || length(representative) != 1L)
          ggchord_stop(caller, ': representative must be a source logical flag')
        out$representative <- representative
        rows[[length(rows)+1L]] <- out
      }
    }
  }
  list(rows = rows, sequences = sequences, parsed = parsed,
    pagination = list(count = parsed$count, next_page = parsed[['next']], previous_page = parsed$previous))
}

#' Read existing protein signature matches as AA features
#'
#' Reads local signature TSV or cached entry JSON from InterPro/Pfam without
#' running an annotation engine or retrieving online data. Accession case and
#' versions remain verbatim, including lowercase accessions returned by the
#' API. Declared full protein lengths supply seq_data; feature endpoints never
#' stand in for full lengths. No source characters or verification are invented.
#'
#' Source score is a list-column preserving each native scalar/string/missing
#' value. Member-library scores and e-values cannot be assumed comparable;
#' no universal confidence, evalue, Identity or read support is manufactured.
#' Original TSV fields and JSON objects remain in the report. TSV field 10,
#' dates, optional integrated/GO/pathway metadata are retained without guessing
#' version-specific semantics or promoting them to plotting instructions.
#'
#' Each JSON location remains one feature row, reusing the existing segments
#' contract for actual fragments and their original indices/status. start/end
#' bound that location; gaps between its segments are not filled. Per-file
#' .source_match and entry/protein/location indices retain the source context.
#' Model and signature IDs are labels, not unique feature-instance IDs; feature_id
#' stays missing. All representative and nonrepresentative matches remain.
#' TSV rows retain their original interval and physical line. No grouping
#' across TSV rows is guessed. Match indices need .source_input when several
#' files are combined; none is an invented biological observation ID.
#'
#' @param file Optional local path vector or glob patterns, before files.
#' @param files Optional additional input paths or patterns.
#' @param format Explicit schema: "scan_tsv" (default, 11-15 columns) or
#'   "entry_json" (cached entry results with metadata/proteins/locations).
#'   JSON requires optional jsonlite. No pagination links are followed;
#'   declared counts and remaining-page links stay in the report.
#' @param source_file Add .source_file to feature rows, default FALSE.
#' @return A list with seq_data (aa), features$signature and report. Source
#'   lengths, native metadata, raw files and pagination remain traceable.
#'   Known length or declared MD5 conflicts for the same accession fail, rather
#'   than silently merging different sources. A declared MD5 is retained, not
#'   independently verified without actual source characters.
#' @seealso [geom_feature()], [geom_feature_label_repel()]
#' @export
#' @examples
#' # Explicit re-encoding of cached AA annotations, not a new analysis run.
#' path <- system.file("extdata", "P00517-signatures.tsv", package = "ggchord")
#' matches <- read_signature_matches(path)
#' matches$features$signature[c("accver", "start", "end", "signature_id")]
read_signature_matches <- function(file = NULL, files = NULL,
    format = c('scan_tsv', 'entry_json'), source_file = FALSE) {
  caller <- 'read_signature_matches()'; format <- match.arg(format)
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, files, 'read_signature_matches')
  if (format == 'entry_json') ggchord_conversion_require('jsonlite', caller)
  tables <- list(); reports <- list(); dictionary <- list(); raw <- list()
  for (i in seq_along(paths)) {
    text <- ggchord_read_text(paths[i], caller)
    input <- if (format == 'scan_tsv') ggchord_signature_tsv(text, caller) else
      ggchord_signature_json(jsonlite::fromJSON(paste(text, collapse = '\n'), simplifyVector = FALSE), caller)
    out <- ggchord_rbind_fill(input$rows)
    if (!is.null(out)) {
      out$.source_input <- rep(i, nrow(out))
      tables[[i]] <- ggchord_add_source_file(out, paths[i], source_file)
    }
    dictionary[[i]] <- ggchord_rbind_fill(input$sequences)
    reports[[i]] <- list(file = paths[i], parsed = input$parsed, pagination = input$pagination)
    raw[[i]] <- text
  }
  features <- ggchord_rbind_fill(tables) %||% data.frame()
  source <- ggchord_rbind_fill(dictionary)
  if (!is.null(source)) for (id in unique(source$accver)) {
    known <- source$source_md5[source$accver == id]
    known <- unique(known[!is.na(known) & nzchar(known) & known != '-'])
    if (length(unique(tolower(known))) > 1L) ggchord_stop(caller, ': conflicting source MD5 declarations for ', id)
  }
  sequence_rows <- lapply(Filter(Negate(is.null), dictionary), function(d)
    data.frame(accver = d$accver, length = d$source_length))
  list(seq_data = ggchord_alignment_dictionary(sequence_rows, 'aa', caller),
    features = list(signature = features),
    report = list(format = format, records = reports, raw_files = stats::setNames(raw, paths),
      coordinates = 'declared AA full lengths and 1-based-inclusive source intervals; no shift',
      score_semantics = 'opaque native member-library evidence; not a common confidence measure'))
}
