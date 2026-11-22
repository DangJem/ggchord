ggchord_validation_issues <- function() {
  data.frame(severity = character(), code = character(), source = character(),
    table = character(), row = integer(), column = character(), message = character())
}

ggchord_data_validation_result <- function(issues, action, caller, checks_run,
                                            checks_skipped, ...) {
  warnings <- issues$severity == "warning"
  errors <- issues$severity == "error"
  result <- structure(list(valid = !any(errors), issues = issues,
    summary = c(errors = sum(errors), warnings = sum(warnings)),
    checks_run = checks_run, checks_skipped = checks_skipped, ...),
    class = "ggchord_validation")
  if (action == "error" && any(errors)) {
    ggchord_abort(paste0(caller, ": ", issues$message[which(errors)[1L]]),
      class = "ggchord_validation_error", report = result)
  }
  if (action == "warn" && nrow(issues)) {
    rlang::warn(paste0(caller, ": ", nrow(issues), " data issue(s). ", issues$message[1L]),
      class = "ggchord_validation_warning", report = result)
  }
  result
}

ggchord_data_endpoints <- function(data, role) {
  switch(role,
    sequence = list(accver = "length"),
    feature = list(accver = c("start", "end")),
    link = if (all(c("qpos", "spos") %in% names(data)) &&
      !any(c("qstart", "qend", "sstart", "send") %in% names(data))) {
      list(qaccver = "qpos", saccver = "spos")
    } else list(qaccver = c("qstart", "qend"), saccver = c("sstart", "send")),
    track = if (any(c("start", "end") %in% names(data))) list(accver = c("start", "end")) else
      list(accver = "pos"),
    alignment = list(accver = "alignment_column", source_accver = "source_position"))
}

ggchord_data_validate_table <- function(data, role, references, enabled, add) {
  if (!is.data.frame(data)) {
    if (enabled("schema")) add("not_dataframe", "Expected a data frame")
    return(invisible(NULL))
  }
  endpoints <- ggchord_data_endpoints(data, role)
  required <- c(names(endpoints), unlist(endpoints, use.names = FALSE),
    switch(role, track = "value", alignment = c("symbol", "is_gap", "alphabet")))
  if (enabled("schema")) {
    for (column in setdiff(required, names(data))) {
      add("missing_column", paste("Missing required column", column), column = column)
    }
    if (anyDuplicated(names(data))) add("duplicate_column", "Column names must be unique")
    for (key in intersect(names(endpoints), names(data))) {
      values <- data[[key]]
      if (is.atomic(values)) {
        add("invalid_id", paste(key, "must contain nonempty IDs"),
          which(is.na(values) | !nzchar(as.character(values))), key)
      } else add("invalid_id", paste(key, "must be an identifier column"), column = key)
    }
    if (role == "sequence" && "accver" %in% names(data)) {
      add("duplicate_id", "Reference IDs must be unique", which(duplicated(data$accver)), "accver")
      if (!nrow(data)) add("empty_reference", "At least one complete reference is required")
    }
    if (role == "sequence") {
      topology <- data[["topology"]]
      circular <- data[["circular"]]
      add("invalid_topology", "topology must be linear or circular",
        which(!is.na(topology) & !topology %in% ggchord_sequence_topologies), "topology")
      if (!is.null(circular) && !is.logical(circular))
        add("invalid_circular", "circular must be logical", column = "circular")
      if (!is.null(topology) && is.logical(circular))
        add("conflicting_topology", "circular and topology disagree",
          which(topology %in% ggchord_sequence_topologies & !is.na(circular) &
            (topology == "circular") != circular), "topology")
    }
    if (role == "feature") {
      strand <- as.character(data[["strand"]])
      invalid <- !is.na(strand) & !strand %in% ggchord_feature_strands
      add("invalid_strand", "strand must be '+', '-', '*', '.', or '+/-'",
        which(invalid), "strand")
      if ("directionality" %in% names(data)) {
        expected <- ggchord_direction_strand(data[["directionality"]])
        add("invalid_directionality", "directionality must be forward, reverse, bidirectional or none",
          which(is.na(expected)), "directionality")
        if (length(strand)) {
          strand[is.na(strand) | strand == "*"] <- "."
          add("conflicting_strand", "strand and directionality describe conflicting directions",
            which(!invalid & !is.na(expected) & strand != expected), "strand")
        }
      }
    }
  }
  reliable_reference <- is.data.frame(references) &&
    all(c("accver", "length") %in% names(references)) && is.numeric(references$length)
  for (key in names(endpoints)) {
    columns <- intersect(endpoints[[key]], names(data))
    if (!length(columns)) next
    ref_match <- if (reliable_reference && key %in% names(data) && is.atomic(data[[key]]) &&
      role %in% c("link", "feature", "track")) match(data[[key]], references$accver)
    if (!is.null(ref_match) && enabled("references"))
      add("unknown_reference", "Observation reference is absent from the dictionary",
        which(is.na(ref_match)), key)
    if (!enabled("coordinates")) next
    for (column in columns) {
      values <- data[[column]]
      if (!is.numeric(values)) {
        add("invalid_coordinate_type", paste(column, "must be numeric"), column = column)
        next
      }
      bad <- !is.finite(values) | values < 1 | values > 2^53 | values != floor(values)
      if (column == "source_position") {
        gaps <- if (is.logical(data$is_gap)) data$is_gap %in% TRUE else rep(FALSE, length(values))
        bad <- (bad & !gaps) | (gaps & !is.na(values))
      }
      add("invalid_coordinate", paste(column, "must contain valid one-based coordinates"),
        which(bad), column)
      if (!is.null(ref_match)) {
        add("out_of_range", paste(column, "is outside its reference"),
          which(values > references$length[ref_match]), column)
      }
    }
    if (length(columns) == 2L && all(vapply(data[columns], is.numeric, logical(1)))) {
      add("reversed_interval", "Interval endpoints must be ordered", which(data[[columns[1L]]] >
        data[[columns[2L]]]), columns[1L])
    }
  }
  if (enabled("coordinates") && "pident" %in% names(data)) {
    values <- data$pident
    if (!is.numeric(values)) add("invalid_pident", "pident must be numeric", column = "pident")
    else add("invalid_pident", "pident must be a percentage in [0, 100]",
      which(!is.na(values) & (!is.finite(values) | values < 0 | values > 100)), "pident")
  }
  if (enabled("coordinates") && role == "feature" && "segments" %in% names(data) &&
      all(c("start", "end") %in% names(data)) && is.numeric(data$start) && is.numeric(data$end)) {
    error <- tryCatch({ ggchord_validate_feature_segments(data, "validate_ggchord_data()"); NULL },
      error = identity)
    if (!is.null(error)) add("invalid_segments", conditionMessage(error), column = "segments")
  }
  if (enabled("coordinates") && role == "alignment" && "source_length" %in% names(data)) {
    source_length <- data[["source_length"]]
    invalid <- ggchord_alignment_invalid_length(source_length)
    add("invalid_source_length", "source_length must be a positive exact integer",
      which(invalid), "source_length")
    position <- data[["source_position"]]
    if (is.numeric(source_length) && is.numeric(position)) {
      add("out_of_range", "source_position exceeds its declared source_length",
        which(!invalid & position > source_length), "source_position")
    }
  }
  if (enabled("schema") && role == "alignment") {
    if ("is_gap" %in% names(data)) {
      if (!is.logical(data$is_gap)) add("invalid_gap", "is_gap must be logical", column = "is_gap")
      else add("invalid_gap", "is_gap must be known for each cell", which(is.na(data$is_gap)), "is_gap")
    }
    if (all(c("symbol", "is_gap") %in% names(data)) && is.logical(data$is_gap)) {
      symbol <- as.character(data$symbol)
      bad <- is.na(symbol) | nchar(symbol) != 1L | ((symbol %in% c(".", "-")) != data$is_gap)
      add("invalid_symbol", "Characters and gap status must agree", which(bad), "symbol")
    }
    if ("alphabet" %in% names(data)) add("invalid_alphabet", "Alphabet must be DNA, RNA, or AA",
      which(!data$alphabet %in% c("DNA", "RNA", "AA")), "alphabet")
    if ("source_scope" %in% names(data)) add("invalid_source_scope", "source_scope must be source or fragment",
      which(!is.na(data$source_scope) & !data$source_scope %in% c("source", "fragment")), "source_scope")
    if (all(c("accver", "alignment_column") %in% names(data))) {
      context <- intersect(c("source_input", ".source_input", "source_block", "accver",
        "alignment_column"), names(data))
      add("duplicate_cell", "Each alignment row and column must identify one cell",
        which(duplicated(data[context])), "alignment_column")
    }
  }
  if (enabled("units")) {
    units <- data[["sequence_unit"]]
    if (!is.null(units)) {
      add("invalid_unit", "Sequence units must be bp or aa", which(!is.na(units) &
        !units %in% c("bp", "aa")), "sequence_unit")
      expected <- references$sequence_unit
      if (!is.null(expected)) {
        for (key in intersect(names(endpoints), names(data))) {
          if (!key %in% c("accver", "qaccver", "saccver") || role %in% c("sequence", "alignment")) next
          unit <- expected[match(data[[key]], references$accver)]
          add("unit_mismatch", "Observation and reference sequence units disagree",
            which(!is.na(units) & !is.na(unit) & units != unit), "sequence_unit")
        }
      }
    }
  }
  if (enabled("provenance") && ".source_row" %in% names(data)) {
    rows <- data$.source_row
    if (!is.numeric(rows)) add("invalid_source_row", "Source rows must be numeric", column = ".source_row")
    else add("invalid_source_row", "Source rows must be positive integers",
      which(!is.na(rows) & (!is.finite(rows) | rows < 1 | rows != floor(rows))), ".source_row")
  }
  if (enabled("duplicates") && nrow(data)) {
    keys <- data[setdiff(names(data), c(".source_row", ".source_line", ".source_file"))]
    if (length(keys)) add("duplicate_row", "Duplicate data record", which(duplicated(keys)),
      severity = "warning")
  }
  invisible(NULL)
}

#' Validate explicit biological data contracts
#'
#' @param x A collection or one data frame.
#' @param role The role of a single table.
#' @param sequences Optional reference data for a single observation table.
#' @param action Report, warn, or signal an error with the complete report.
#' @param checks Which checks to run. The schema check includes declared
#'   reference topology and feature strand and directionality.
#'   The coordinates check includes declared alignment source lengths.
#' @return A validation result with source and row information for each issue.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' bundle <- ggchord_data(seqs, features=list(genes=genes), links=list(blast=links))
#' report <- validate_ggchord_data(bundle)
#' report$issues
#' validate_ggchord_data(genes, role="feature", sequences=seqs, action="error")
validate_ggchord_data <- function(x, role = NULL, sequences = NULL, action = "report",
                                  checks = c("schema", "coordinates", "units", "references", "provenance")) {
  ggchord_local_console_errors()
  caller <- "validate_ggchord_data()"
  action <- ggchord_match_arg(action, c("report", "warn", "error"))
  roles <- c(features = "feature", links = "link", tracks = "track", alignments = "alignment")
  allowed_checks <- c("schema", "coordinates", "units", "references", "provenance", "duplicates")
  if (!is.character(checks) || anyNA(checks) || any(!checks %in% allowed_checks))
    ggchord_stop(caller, ": unknown check")
  checks <- unique(checks)
  grouped <- ggchord_is_data_collection(x)
  if (grouped && (!is.null(role) || !is.null(sequences))) {
    ggchord_stop(caller, ": a collection defines its own roles and references")
  }
  if (grouped) sequences <- x$sequences
  if (!grouped && is.null(role)) ggchord_stop(caller, ": role is required for a single table")
  if (!grouped) role <- ggchord_match_arg(role, c("sequence", unname(roles)))
  if (!grouped && identical(role, "sequence")) sequences <- x
  skipped <- if (is.null(sequences)) intersect(checks, "references") else character()
  checks_run <- setdiff(checks, skipped)
  enabled <- function(check) check %in% checks_run
  issues <- list()
  check_table <- function(data, role, table, source = table) {
    add <- function(code, message, row = NA_integer_, column = NA_character_, severity = "error") {
      if (!length(row)) return(invisible(NULL))
      issues[[length(issues) + 1L]] <<- data.frame(severity = severity, code = code,
        source = source, table = table, row = as.integer(row), column = column, message = message)
      invisible(NULL)
    }
    ggchord_data_validate_table(data, role, sequences, enabled, add)
  }
  if (grouped) {
    check_table(sequences, "sequence", "sequences")
    for (collection in names(roles)) {
      tables <- x[[collection]]
      for (table in names(tables)) check_table(tables[[table]], roles[[collection]],
        table, paste(collection, table, sep = "/"))
    }
  } else {
    check_table(x, role, role)
  }
  issues <- if (length(issues)) do.call(rbind, issues) else ggchord_validation_issues()
  ggchord_data_validation_result(issues, action, caller, checks_run, skipped)
}

#' Validate an explicit genome and contig hierarchy
#' @param sequences Reference dictionary with genome_id foreign keys.
#' @param genomes Genome table keyed by genome_id.
#' @param action Report, warn, or signal an error with the complete report.
#' @return A validation result and a table of contigs belonging to each genome.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' ref <- seqs
#' ref$genome_id <- c("genomeA", "genomeB")
#' genomes <- data.frame(genome_id=c("genomeA", "genomeB"))
#' validate_ggchord_hierarchy(ref, genomes)
validate_ggchord_hierarchy <- function(sequences, genomes, action = "report") {
  ggchord_local_console_errors()
  action <- ggchord_match_arg(action, c("report", "warn", "error"))
  reference <- validate_ggchord_data(sequences, role = "sequence")
  seq_key <- sequences
  genome_key <- genomes
  if (is.data.frame(seq_key)) {
    seq_key$genome <- seq_key[["genome_id"]]
  }
  if (is.data.frame(genome_key)) {
    genome_key$genome <- genome_key[["genome_id"]]
  }
  relations <- ggchord_validate_hierarchy_keys(seq_key, genome_key)
  issues <- relations$issues
  issues$column[issues$column %in% "genome"] <- "genome_id"
  issues$code <- issues$category
  issues$source <- ifelse(issues$table == "seq", "sequences", "genomes")
  issues$table <- issues$source
  hierarchy <- relations$contigs
  names(hierarchy)[names(hierarchy) == "genome"] <- "genome_id"
  names(hierarchy)[names(hierarchy) == "n_contigs"] <- "n_sequences"
  names(hierarchy)[names(hierarchy) == "sequence_rows"] <- "source_rows"
  ggchord_data_validation_result(rbind(reference$issues, issues[names(reference$issues)]),
    action, "validate_ggchord_hierarchy()", c(reference$checks_run, "hierarchy"),
    reference$checks_skipped, hierarchy = hierarchy)
}
