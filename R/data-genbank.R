# Mature optional parsing, with an explicit process boundary and provenance.
# Parser documentation: https://biopython.org/docs/latest/api/Bio.SeqIO.html
# Locations: https://biopython.org/docs/latest/api/Bio.SeqFeature.html

ggchord_record_parser_script <- function() r"---(
import json, sys
try:
    import Bio
    from Bio import SeqIO
except ImportError:
    raise SystemExit('The selected Python environment needs the optional biopython package; install it there, or import FASTA/GFF3/BED tables instead.')

with open(sys.argv[1], encoding='utf-8') as handle:
    config = json.load(handle)

def preserve(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {str(key): preserve(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [preserve(item) for item in value]
    result = {'source_class': type(value).__name__, 'representation': str(value)}
    if hasattr(value, '__dict__'):
        result['fields'] = preserve(vars(value))
    return result

def boundary(value, start=False):
    try:
        coordinate = int(value) + (1 if start else 0)
    except (TypeError, ValueError):
        coordinate = None
    return {'coordinate': coordinate, 'source': str(value),
            'source_class': type(value).__name__}

def parse_local(path):
    with open(path, encoding='utf-8-sig') as handle:
        yield from SeqIO.parse(handle, 'genbank')

records = []
for path in config['files']:
    for index, record in enumerate(parse_local(path), 1):
        features = []
        for number, feature in enumerate(record.features, 1):
            location = feature.location
            parts = []
            if location is not None:
                for part in location.parts:
                    parts.append({'start': boundary(part.start, True),
                        'end': boundary(part.end), 'strand': part.strand,
                        'ref': part.ref, 'ref_db': part.ref_db,
                        'source_location': str(part)})
            features.append({'source_feature': number, 'source_id': feature.id,
                'type': feature.type, 'qualifiers': preserve(feature.qualifiers),
                'source_location': str(location),
                'operator': getattr(location, 'operator', None), 'parts': parts})
        try:
            length = len(record.seq)
        except (TypeError, ValueError):
            length = None
        records.append({'accver': record.id, 'name': record.name,
            'description': record.description, 'length': length,
            'sequence_defined': getattr(record.seq, 'defined', None),
            'annotations': preserve(record.annotations),
            'dbxrefs': preserve(record.dbxrefs), 'features': features,
            'source_file': path, 'source_record': index})
with open(sys.argv[2], 'w', encoding='utf-8') as handle:
    json.dump({'parser_version': Bio.__version__, 'records': records}, handle,
              ensure_ascii=False, allow_nan=False)
)---"

#' Import local GenBank or GBFF records using an optional mature parser
#'
#' Uses the public `Bio.SeqIO.parse(..., "genbank")` API in an explicitly
#' selected Python environment with Biopython installed. The R package
#' `jsonlite` exchanges structured temporary JSON with that process. This
#' function never installs packages, downloads records, or runs during plot
#' building. It reads all selected records into memory and is intended for
#' small to moderate local annotation files, not chromosome-scale streaming.
#'
#' Record IDs, including accession versions, are retained exactly as returned
#' by the source parser. Full parsed record annotations, feature qualifiers,
#' source locations, and original file text are retained in `report`. Feature
#' IDs use explicit `ID`, then `protein_id` for CDS, then `locus_tag`, then a
#' non-placeholder parser feature ID. Missing IDs remain missing. Shared
#' `locus_tag` values do not invent gene/CDS parent relationships; explicit
#' `Parent` qualifiers become `parent_ids`.
#'
#' Simple and joined local locations retain one-based inclusive `segments`
#' in parser order. Before/after fuzzy bounds retain their source descriptors
#' and partial flags; they are displayed at the reported boundary rather than
#' treated as exact biology. Remote locations, missing or unsupported fuzzy
#' positions, zero-width sites, and mixed-strand joined features are never
#' flattened into guessed intervals. `unsupported = "error"` stops by default;
#' `"report"` retains these features in `report$unplaced_features` instead.
#' An absent ORIGIN sequence does not imply an unknown declared length.
#' Truly missing/zero lengths require explicit `lengths`; feature endpoints
#' are never used to infer full lengths. Protein records require an explicit
#' protein-coordinate workflow instead of this nucleotide reader.
#'
#' @param file Optional path vector or glob patterns; read before `files`.
#' @param files Optional paths or wildcard patterns for several files.
#' @param feature_types Optional feature types to select. `NULL` retains all.
#' @param lengths Optional named numeric vector of known full sequence lengths,
#'   following [as_ggchord_seq()].
#' @param python Interpreter name or path, default `"python3"`. It must already
#'   have the optional Biopython package installed. No environment is modified.
#' @param unsupported `"error"` or `"report"` for locations without a safe
#'   local interval representation.
#' @param source_file Add `.source_file` columns to output tables, default
#'   `FALSE`. File provenance is retained in the report regardless.
#' @return A list with `seq_data`, named `features$annotation`, and `report`.
#'   Tables can be supplied to [prepare_ggchord_data()]. Feature rows include
#'   `feature_id`, `parent_ids`, `segments`, source record/feature indices,
#'   `qualifiers`, `source_location`, `partial_start`/`partial_end`, and
#'   `location_operator`. `report$records` retains the complete parsed source
#'   CDS rows also expose source `codon_start`, `genetic_code`, biological
#'   `partial5`/`partial3`, and `translation_exception` qualifiers for explicit
#'   downstream coding-coordinate models. Missing CDS `codon_start` uses the
#'   INSDC default 1, marked by `codon_start_source = "format-default"`.
#'   Bound-only flags remain in `bound_partial5`/`bound_partial3`; CDS partial
#'   flags also record initial codon offset or terminal coding remainder, with
#'   their evidence in `partial5_source`/`partial3_source`. They describe the
#'   observed coding extent, not proof of a complete biological protein.
#'   No per-segment phase values are invented. `segment_order` preserves the
#'   source's declared coding-piece order independently of `segment_index`.
#'   metadata, `raw_files` retains original text, and `unplaced_features`
#'   identifies every selected feature not placed under `"report"`.
#' @seealso <https://biopython.org/docs/latest/api/Bio.SeqIO.html>
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' if (requireNamespace("jsonlite", quietly = TRUE) &&
#'     nzchar(Sys.which("python3"))) {
#'   available <- suppressWarnings(system2("python3",
#'     c("-c", shQuote("from Bio import SeqIO")), stdout = FALSE, stderr = FALSE))
#'   if (identical(available, 0L)) {
#'     path <- tempfile(fileext = ".gb")
#'     writeLines(c(
#'       "LOCUS       demo                      12 bp    DNA     linear   BCT 01-JAN-2000",
#'       "DEFINITION  Demonstration record.", "ACCESSION   demo", "VERSION     demo.1",
#'       "FEATURES             Location/Qualifiers",
#'       "     misc_feature    2..8", "                     /note=\"site\"",
#'       "ORIGIN", "        1 acgtacgtacgt", "//"), path)
#'     imported <- read_genbank(path)
#'     prepared <- prepare_ggchord_data(imported$seq_data, features = imported$features)
#'     prepared$seq_data
#'     unlink(path)
#'   }
#' }
read_genbank <- function(file = NULL, files = NULL, feature_types = NULL,
    lengths = NULL, python = "python3", unsupported = c("error", "report"),
    source_file = FALSE) {
  caller <- "read_genbank()"
  unsupported <- match.arg(unsupported)
  paths <- ggchord_resolve_files(file, files, "read_genbank")
  if (!is.null(feature_types) && (!is.character(feature_types) || !length(feature_types) ||
      anyNA(feature_types) || any(!nzchar(feature_types)))) ggchord_stop(
    caller, ": feature_types must be NULL or nonempty feature type names")
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ": source_file must be TRUE or FALSE")
  if (!is.character(python) || length(python) != 1L || is.na(python) || !nzchar(python))
    ggchord_stop(caller, ": python must name one existing interpreter")
  executable <- if (file.exists(python)) normalizePath(python) else unname(Sys.which(python))
  if (!nzchar(executable)) ggchord_stop(caller,
    ": selected Python interpreter was not found; select an existing environment with the optional parser, or use FASTA/GFF3/BED input")
  ggchord_conversion_require("jsonlite", caller)
  temporary <- tempfile("ggchord-records-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  script <- file.path(temporary, "parse_records.py")
  config <- file.path(temporary, "input.json")
  output <- file.path(temporary, "output.json")
  log <- file.path(temporary, "parser.log")
  writeLines(ggchord_record_parser_script(), script, useBytes = TRUE)
  jsonlite::write_json(list(files = unname(normalizePath(paths))), config, auto_unbox = FALSE)
  status <- suppressWarnings(system2(executable,
    args = c(shQuote(script), shQuote(config), shQuote(output)), stdout = log, stderr = log))
  messages <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
  if (!identical(status, 0L) || !file.exists(output)) ggchord_stop(caller,
    ": optional source parser failed in the selected environment: ", paste(messages, collapse = "\n"))
  parsed <- jsonlite::read_json(output, simplifyVector = FALSE)
  records <- parsed$records
  if (!length(records)) ggchord_stop(caller, ": no GenBank records were parsed")
  seqs <- list(); features <- list(); unplaced <- list(); feature_files <- character()
  strand_label <- function(x) if (is.null(x)) "." else if (x == 1) "+" else if (x == -1) "-" else "."
  qualifier <- function(x, key) {
    value <- x[[key]]
    if (is.null(value) || !length(value)) return(NA_character_)
    as.character(value[[1L]])
  }
  integer_qualifier <- function(x, key, maximum = Inf) {
    value <- qualifier(x, key)
    if (is.na(value)) return(NA_real_)
    number <- suppressWarnings(as.numeric(value))
    if (!grepl("^[0-9]+$", value) || !is.finite(number) || number < 1 || number > maximum)
      ggchord_stop(caller, ": invalid source ", key, " qualifier: ", value)
    number
  }
  for (record in records) {
    accver <- record$accver
    if (is.null(accver) || !nzchar(accver) || accver == "<unknown id>")
      ggchord_stop(caller, ": each source record needs an explicit identifier")
    if (identical(tolower(record$annotations$molecule_type %||% ""), "protein"))
      ggchord_stop(caller, ": protein record '", accver, "' requires an explicit protein-coordinate workflow")
    length <- record$length %||% NA_real_
    if (length == 0 || is.na(length)) length <- NA_real_
    topology <- record$annotations$topology %||% NA_character_
    seq_row <- data.frame(accver = accver, length = length,
      circular = if (is.na(topology)) NA else topology == "circular",
      name = record$name, description = record$description,
      sequence_defined = record$sequence_defined %||% NA, .source_record = record$source_record,
      stringsAsFactors = FALSE)
    seqs[[length(seqs) + 1L]] <- ggchord_add_source_file(seq_row, record$source_file, source_file)
    for (feature in record$features) {
      if (!is.null(feature_types) && !feature$type %in% feature_types) next
      reason <- character()
      parts <- feature$parts
      if (!length(parts)) reason <- "missing location"
      if (length(parts)) {
        if (any(vapply(parts, function(x) !is.null(x$ref) || !is.null(x$ref_db), logical(1))))
          reason <- c(reason, "remote location")
        classes <- unlist(lapply(parts, function(x) c(x$start$source_class, x$end$source_class)))
        if (any(!classes %in% c("ExactPosition", "BeforePosition", "AfterPosition")))
          reason <- c(reason, "unsupported fuzzy boundary")
        if (any(vapply(parts, function(x) is.null(x$start$coordinate) ||
            is.null(x$end$coordinate) || x$end$coordinate < x$start$coordinate, logical(1))))
          reason <- c(reason, "missing or zero-width interval")
        strands <- vapply(parts, function(x) strand_label(x$strand), character(1))
        if (length(unique(strands)) > 1L) reason <- c(reason, "mixed strands")
      }
      if (length(reason)) {
        entry <- list(accver = accver, source_file = record$source_file,
          source_record = record$source_record, source_feature = feature$source_feature,
          reason = paste(unique(reason), collapse = "; "), feature = feature)
        if (unsupported == "error") ggchord_stop(caller, ": ", accver,
          " feature ", feature$source_feature, ": ", entry$reason,
          "; use unsupported = 'report' to retain it without guessed geometry")
        unplaced[[length(unplaced) + 1L]] <- entry
        next
      }
      segments <- do.call(rbind, lapply(seq_along(parts), function(i) {
        part <- parts[[i]]
        data.frame(start = part$start$coordinate, end = part$end$coordinate,
          strand = strand_label(part$strand), segment_index = i, segment_order = i,
          segment_type = "standard", source_start = part$start$source,
          source_end = part$end$source, start_class = part$start$source_class,
          end_class = part$end$source_class, stringsAsFactors = FALSE)
      }))
      ids <- c(ID = qualifier(feature$qualifiers, "ID"),
        protein_id = if (feature$type == "CDS") qualifier(feature$qualifiers, "protein_id") else NA_character_,
        locus_tag = qualifier(feature$qualifiers, "locus_tag"),
        source_id = if (feature$source_id == "<unknown id>") NA_character_ else feature$source_id)
      ids <- ids[!is.na(ids) & nzchar(ids)]
      labels <- vapply(c("product", "gene", "locus_tag"), function(key)
        qualifier(feature$qualifiers, key), character(1))
      labels <- labels[!is.na(labels) & nzchar(labels)]
      row <- data.frame(accver = accver, start = min(segments$start), end = max(segments$end),
        strand = segments$strand[1L], type = feature$type,
        feature_id = if (length(ids)) unname(ids[1L]) else NA_character_,
        feature_id_source = if (length(ids)) names(ids)[1L] else NA_character_,
        anno = if (length(labels)) labels[1L] else NA_character_,
        source_location = feature$source_location,
        location_operator = feature$operator %||% NA_character_,
        partial_start = any(segments$start_class != "ExactPosition"),
        partial_end = any(segments$end_class != "ExactPosition"),
        .source_record = record$source_record, .source_feature = feature$source_feature,
        stringsAsFactors = FALSE)
      row$segments <- I(list(segments))
      row$codon_start <- integer_qualifier(feature$qualifiers, "codon_start", 3)
      row$codon_start_source <- if (is.na(row$codon_start)) "missing" else "qualifier"
      if (feature$type == "CDS" && is.na(row$codon_start) &&
          is.na(qualifier(feature$qualifiers, "codon_start"))) {
        row$codon_start <- 1L
        row$codon_start_source <- "format-default"
      }
      row$genetic_code <- integer_qualifier(feature$qualifiers, "transl_table")
      first_part <- segments[1L, ]
      last_part <- segments[nrow(segments), ]
      row$partial5 <- if (row$strand == "-") first_part$end_class != "ExactPosition" else
        if (row$strand == "+") first_part$start_class != "ExactPosition" else NA
      row$partial3 <- if (row$strand == "-") last_part$start_class != "ExactPosition" else
        if (row$strand == "+") last_part$end_class != "ExactPosition" else NA
      row$bound_partial5 <- row$partial5
      row$bound_partial3 <- row$partial3
      row$partial5_source <- if (isTRUE(row$partial5)) "boundary" else "no-recorded-partial-evidence"
      row$partial3_source <- if (isTRUE(row$partial3)) "boundary" else "no-recorded-partial-evidence"
      if (feature$type == "CDS" && !is.na(row$codon_start) && row$codon_start %in% 1:3) {
        if (row$codon_start > 1L) {
          row$partial5 <- TRUE
          row$partial5_source <- if (isTRUE(row$bound_partial5)) "boundary+codon-start" else "codon-start"
        }
        coding_width <- sum(segments$end - segments$start + 1) - row$codon_start + 1
        if (coding_width > 0 && coding_width %% 3 != 0) {
          row$partial3 <- TRUE
          row$partial3_source <- if (isTRUE(row$bound_partial3)) "boundary+coding-remainder" else "coding-remainder"
        }
      }
      exception <- unlist(feature$qualifiers[intersect(c("transl_except", "exception"),
        names(feature$qualifiers))], use.names = FALSE)
      row$translation_exception <- if (length(exception)) paste(exception, collapse = "; ") else NA_character_
      row$qualifiers <- I(list(feature$qualifiers))
      parents <- feature$qualifiers$Parent %||% list()
      row$parent_ids <- I(list(unlist(lapply(parents, function(x)
        strsplit(as.character(x), ",", fixed = TRUE)[[1L]]), use.names = FALSE)))
      features[[length(features) + 1L]] <- ggchord_add_source_file(row, record$source_file, source_file)
      feature_files <- c(feature_files, record$source_file)
    }
  }
  seq_data <- as_ggchord_seq(do.call(rbind, seqs), lengths = lengths)
  feature_data <- if (length(features)) ggchord_rbind_fill(features) else data.frame(
    accver = character(), start = numeric(), end = numeric(), strand = character(),
    type = character(), feature_id = character(), stringsAsFactors = FALSE)
  if (nrow(feature_data)) {
    reference <- seq_data$length[match(feature_data$accver, seq_data$accver)]
    if (any(feature_data$start < 1 | feature_data$end > reference))
      ggchord_stop(caller, ": parsed feature intervals lie outside their declared sequence lengths")
  }
  provenance <- list(format = "GenBank", files = paths,
    parser = "Bio.SeqIO.parse", parser_version = parsed$parser_version,
    input_coordinates = "INSDC", output_coordinates = "1-based-inclusive")
  record_files <- vapply(records, function(x) x$source_file, character(1))
  sequence_provenance <- feature_provenance <- provenance
  sequence_provenance$rows_per_file <- vapply(normalizePath(paths), function(path)
    sum(record_files == path), integer(1))
  feature_provenance$rows_per_file <- vapply(normalizePath(paths), function(path)
    sum(feature_files == path), integer(1))
  attr(seq_data, "ggchord_import") <- sequence_provenance
  attr(feature_data, "ggchord_import") <- feature_provenance
  list(seq_data = seq_data, features = list(annotation = feature_data), report = list(
    parser = provenance, records = records,
    raw_files = stats::setNames(lapply(paths, readLines, warn = FALSE), paths),
    parser_messages = messages, unplaced_features = unplaced))
}
