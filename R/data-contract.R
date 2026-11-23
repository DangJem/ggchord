#' Biological data and source coordinates
#'
#' Sequence tables carry stable reference IDs, complete lengths, coordinate
#' units and biological topology. Reference lengths are supplied by the
#' source or declared explicitly. Interval annotations and links use ordered,
#' one-based closed coordinates. Site annotations use a one-based residue
#' position. Cut annotations carry a residue anchor and may have separate
#' interbase boundaries. Biological strands and relative alignment direction
#' remain separate from endpoints and display orientation.
#'
#' [ggchord_data()] optionally organizes a reference and named observation
#' tables. Each drawing layer uses its own observation table. [data_report()]
#' retrieves source information and operation results, including the members
#' of aggregated relationships and the original positions of focused windows.
#' Annotation reports distinguish original rows, entities and segments.
#' Without a source-backed or supplied identifier, entity keys are local to
#' the input table. Track location tables can describe intervals or points
#' independently of their optional values.
#'
#' [validate_ggchord_data()] reports input issues and [clean_ggchord_data()]
#' applies explicitly chosen repairs. [focus_ggchord_data()] builds new local
#' data coordinates. A layout display window leaves the reference dictionary
#' and source statistics in their original coordinate system.
#'
#' Character alignments use real alignment columns and nullable source
#' positions. A gap consumes a column without occupying a source residue.
#' Restriction cut boundaries are interbase coordinates, from zero to the
#' complete reference length. Source coordinate metadata records unit and
#' origin without relying on the selected plotting coordinate system.
#'
#' @name ggchord_data_contract
NULL

# Resolve the native coordinate domain without inferring it from alphabet or IDs.
ggchord_sequence_unit <- function(seq_data, caller = "ggchord") {
  if (!is.data.frame(seq_data) || !"sequence_unit" %in% names(seq_data)) return("bp")
  units <- as.character(seq_data$sequence_unit)
  if (!length(units) || anyNA(units) || any(!units %in% c("bp", "aa"))) {
    ggchord_stop(caller, ": sequence_unit must be 'bp' or 'aa'")
  }
  if (length(unique(units)) != 1L) ggchord_stop(caller,
    ": bp and aa cannot share one native coordinate domain; project protein intervals to CDS coordinates or use separate plots")
  units[1L]
}

# Apply explicit interval declarations to the sequence coordinate domain.
ggchord_check_interval_unit <- function(data, unit, caller) {
  if (is.data.frame(data) && nrow(data) && "sequence_unit" %in% names(data) &&
      !identical(ggchord_sequence_unit(data, caller), unit)) {
    ggchord_stop(caller, ": interval sequence_unit (", ggchord_sequence_unit(data, caller),
      ") disagrees with the sequence table (", unit, ")")
  }
  invisible(NULL)
}

ggchord_sequence_topologies <- c("linear", "circular")
ggchord_feature_strands <- c("+", "-", "*", ".", "+/-")

ggchord_annotation_locations <- c("interval", "site", "cut")

ggchord_annotation_location <- function(data, location = "auto", caller = "Annotation") {
  location <- ggchord_match_choice(location, c("auto", ggchord_annotation_locations),
    paste0(caller, ": location"))
  declared <- data[["location_kind"]]
  if (!is.null(declared) && (anyNA(declared) ||
      any(!declared %in% ggchord_annotation_locations) || length(unique(declared)) > 1L)) {
    ggchord_stop(caller, ": location_kind must declare one interval, site or cut type per table")
  }
  declared <- if (length(declared)) as.character(declared[1L]) else
    attr(data, "ggchord_annotation_location", exact = TRUE)
  if (location != "auto" && !is.null(declared) && location != declared)
    ggchord_stop(caller, ": location disagrees with the table's location_kind: ", declared)
  interval <- any(c("start", "end") %in% names(data))
  point <- "pos" %in% names(data)
  cut <- any(c("cut_top", "cut_bottom", "cut_status") %in% names(data))
  if (location == "auto" && !is.null(declared)) location <- declared
  if (location == "auto") {
    if (interval && (point || cut)) ggchord_stop(caller,
      ": ambiguous annotation location: interval and site or cut columns coexist; declare location explicitly")
    location <- if (cut) "cut" else if (point) "site" else "interval"
  }
  if (location != "cut" && cut) ggchord_stop(caller,
    ": interbase cut fields require location = 'cut'")
  location
}

ggchord_annotation_coordinates <- function(location) {
  if (location == "interval") c("start", "end") else "pos"
}

ggchord_annotation_cut_columns <- c("cut_top", "cut_bottom")

ggchord_interbase_numeric <- function(x) {
  is.numeric(x) || (is.logical(x) && all(is.na(x)))
}

ggchord_invalid_interbase <- function(x) {
  is.nan(x) | (!is.na(x) & (!is.finite(x) | x < 0 | x > 2^53 | x != floor(x)))
}

ggchord_annotation_cut_status <- function(data) {
  top <- data[["cut_top"]] %||% rep(NA_real_, nrow(data))
  bottom <- data[["cut_bottom"]] %||% rep(NA_real_, nrow(data))
  c("unknown", "partial", "complete")[1L + (!is.na(top)) + (!is.na(bottom))]
}

ggchord_annotation_import <- function(data) {
  import <- attr(data, "ggchord_import", exact = TRUE)
  conversion <- attr(data, "ggchord_conversion", exact = TRUE)
  while (is.null(import) && !is.null(conversion)) {
    import <- conversion$source_import
    conversion <- conversion$source_conversion
  }
  import
}

ggchord_annotation_source_files <- function(data) {
  import <- ggchord_annotation_import(data)
  files <- import$annotation_files %||% import$files
  path <- rep(NA_character_, nrow(data))
  input <- data[[".source_input"]]
  if (is.numeric(input) && length(files)) {
    legal <- !is.na(input) & input >= 1 & input <= length(files) & input == floor(input)
    path[legal] <- files[input[legal]]
  } else if (length(files) == 1L) {
    path[] <- files
  }
  supplied <- data[[".source_file"]]
  if (!is.null(supplied)) {
    supplied <- as.character(supplied)
    known <- !is.na(supplied) & nzchar(supplied)
    if (any(known)) {
      distinct <- unique(supplied[known])
      resolved <- normalizePath(distinct, winslash = "/", mustWork = FALSE)
      # Reader metadata already contains absolute paths even if the row's
      # optional path was originally supplied relative to a different folder.
      resolved_rows <- resolved[match(supplied[known], distinct)]
      absolute <- grepl("^(/|[[:alpha:]]:/)", resolved_rows)
      path[known] <- ifelse(absolute | is.na(path[known]), resolved_rows, path[known])
    }
  }
  path
}

ggchord_annotation_key <- function(...) {
  fields <- list(...)
  if (!length(fields[[1L]])) return(character())
  parts <- lapply(fields, function(x) {
    x <- enc2utf8(as.character(x))
    result <- paste0(nchar(x, type = "bytes"), ":", x)
    result[is.na(x)] <- "-:"
    result
  })
  do.call(paste0, parts)
}

ggchord_annotation_identity <- function(data, location, namespace = NULL,
                                        caller = "Annotation") {
  n <- nrow(data)
  source_ref <- data[[".source_accver"]] %||% data[["source_accver"]]
  original_ref <- as.character(source_ref %||% data$accver)
  unresolved_ref <- is.na(original_ref) | !nzchar(original_ref)
  original_ref[unresolved_ref] <- as.character(data$accver[unresolved_ref])
  displayed <- !is.null(source_ref) && (".source_origin" %in% names(data) ||
    "source_origin" %in% names(data))
  source_row <- seq_len(n)
  for (column in intersect(c(".source_row", ".source_record", ".source_feature", ".source_line"),
      names(data))) {
    value <- data[[column]]
    supplied <- !is.na(value)
    source_row[supplied] <- value[supplied]
  }
  file <- ggchord_annotation_source_files(data)
  known_namespace <- file
  record <- data[[".source_record"]]
  feature <- data[[".source_feature"]]
  if (!is.null(feature) && !is.null(record)) {
    known <- !is.na(file) & !is.na(record) & !is.na(feature)
    known_namespace[known] <- ggchord_annotation_key(file[known], record[known])
  }
  scope <- known_namespace
  scope[is.na(scope)] <- namespace %||% ""
  id <- id_column <- rep(NA_character_, n)
  preferred <- if (location == "interval") c("feature_id", "site_id") else c("site_id", "feature_id")
  for (column in intersect(c("entity_id", preferred), names(data))) {
    value <- data[[column]]
    if (!is.atomic(value) || !is.null(dim(value))) ggchord_stop(caller,
      ": ", column, " must be an identifier column")
    present <- is.na(id) & !is.na(value) & nzchar(as.character(value))
    id[present] <- as.character(value[present])
    id_column[present] <- column
  }
  key <- ggchord_annotation_key(scope, original_ref, location, id)
  ambiguous <- rep(FALSE, n)
  repeated <- !is.na(id) & (duplicated(key) | duplicated(key, fromLast = TRUE))
  groups <- split(which(repeated), key[repeated])
  gff3 <- identical(ggchord_annotation_import(data)$format, "GFF3")
  for (rows in groups) {
    if (displayed && !anyDuplicated(data.frame(reference = data$accver[rows],
        source_row = source_row[rows]))) {
      rows_to_check <- rows[!duplicated(source_row[rows])]
      if (length(rows_to_check) == 1L) next
    } else {
      rows_to_check <- rows
    }
    # Rows remain separate. Establish a common identity only for explicitly
    # declared entities or source-backed, consistent discontinuous GFF3 parts.
    same_parts <- location == "interval" && all(c("start", "end") %in% names(data)) &&
      (all(id_column[rows] == "entity_id") || (gff3 && all(!is.na(file[rows]))))
    if (same_parts) {
      type <- data[["type"]] %||% rep(NA_character_, n)
      strand <- data[["source_strand"]] %||% data[["strand"]] %||% rep(NA_character_, n)
      source <- data[["source"]] %||% rep(NA_character_, n)
      parents <- lapply((data[["parent_ids"]] %||% rep(list(NULL), n))[rows_to_check], function(x) {
        sort(unique(as.character(x)))
      })
      start <- (data[[".source_start"]] %||% data$start)[rows_to_check]
      end <- (data[[".source_end"]] %||% data$end)[rows_to_check]
      order <- order(start, end)
      same_parts <- length(unique(type[rows_to_check])) == 1L && length(unique(strand[rows_to_check])) == 1L &&
        length(unique(source[rows_to_check])) == 1L && length(unique(parents)) == 1L &&
        all(is.finite(start) & is.finite(end) & start <= end) &&
        all(head(end[order], -1L) < tail(start[order], -1L))
    }
    if (!same_parts) ambiguous[rows] <- TRUE
  }
  temporary <- is.na(id) | ambiguous
  if (any(temporary)) {
    row_key <- as.character(source_row[temporary])
    distinct_rows <- ggchord_annotation_key(scope[temporary], data$accver[temporary], row_key)
    duplicate_row <- is.na(row_key) | duplicated(distinct_rows) | duplicated(distinct_rows, fromLast = TRUE)
    row_key[duplicate_row] <- ggchord_annotation_key(row_key[duplicate_row], which(temporary)[duplicate_row])
    key[temporary] <- ggchord_annotation_key(scope[temporary], original_ref[temporary], location,
      NA_character_, row_key)
  }
  data.frame(key = key, id = id, id_column = id_column, namespace = known_namespace,
    source_file = file, source_row = source_row, ambiguous = ambiguous, stringsAsFactors = FALSE)
}

ggchord_annotation_lineage <- function(data, identity) {
  out <- data.frame(annotation_row = seq_len(nrow(data)), entity_key = identity$key,
    entity_id = identity$id, id_column = identity$id_column, ambiguous_id = identity$ambiguous,
    accver = data$accver, source_file = identity$source_file, source_row = identity$source_row,
    stringsAsFactors = FALSE)
  fields <- intersect(c(".source_input", ".source_record", ".source_feature", ".source_line",
    ".source_row", ".source_accver", "source_accver", "source_members", "feature_id", "site_id",
    "parent_ids", "phase", "source_start", "source_end", "source_position", "source_cut_top",
    "source_cut_bottom", "source_cut_status"), names(data))
  for (field in fields) out[[field]] <- data[[field]]
  out
}

ggchord_annotation_segment_index <- function(data, identity) {
  nested <- if (is.null(data$segments)) integer() else which(lengths(data$segments) > 0L)
  plain <- setdiff(seq_len(nrow(data)), nested)
  rows <- function(index) {
    data.frame(entity_key = identity$key[index], annotation_row = index,
      source_file = identity$source_file[index], source_row = identity$source_row[index],
      stringsAsFactors = FALSE)
  }
  basic <- rows(plain)
  basic$start <- data$start[plain]
  basic$end <- data$end[plain]
  basic$strand <- data$strand[plain]
  basic$segment_index <- rep(1L, length(plain))
  if ("phase" %in% names(data)) basic$phase <- data$phase[plain]
  if (length(nested)) {
    specific <- lapply(nested, function(i) {
      parts <- data$segments[[i]]
      out <- rows(rep(i, nrow(parts)))
      for (name in setdiff(names(parts), names(out))) out[[name]] <- parts[[name]]
      out$segment_index <- out$segment_index %||% seq_len(nrow(out))
      if (is.null(out$strand)) out$strand <- rep(data$strand[i], nrow(out))
      out
    })
    basic <- ggchord_rbind_fill(c(list(basic), specific))
    basic <- basic[order(basic$annotation_row), , drop = FALSE]
  }
  # Keep original segment numbers. A separate ordinal connects parts across
  # different source rows when one entity comprises more than one row.
  basic$entity_segment_index <- basic$segment_index
  multiple <- duplicated(identity$key) | duplicated(identity$key, fromLast = TRUE)
  for (entity in split(which(multiple[basic$annotation_row]), basic$entity_key[multiple[basic$annotation_row]])) {
    source <- data.frame(row = basic$source_row[entity], index = basic$segment_index[entity])
    key <- ggchord_annotation_key(source$row, source$index)
    basic$entity_segment_index[entity] <- match(key, unique(key[order(source$row, source$index)]))
  }
  basic
}

# Check declared units against a scalar domain or the reference of each row.
# A missing declaration leaves the domain to the calling operation.
ggchord_check_declared_unit <- function(units, expected, caller, column = "sequence_unit") {
  units <- as.character(units)
  invalid <- which(!is.na(units) & !units %in% c("bp", "aa"))
  if (length(invalid)) ggchord_stop(caller, ": ", column, " must be 'bp' or 'aa' at row ", invalid[1L])
  if (length(units) && !is.null(expected)) {
    expected <- rep_len(expected, length(units))
    mismatch <- which(!is.na(units) & !is.na(expected) & units != expected)
    if (length(mismatch)) {
      row <- mismatch[1L]
      ggchord_stop(caller, ": ", column, " (", units[row], ") disagrees with the expected unit (",
        expected[row], ") at row ", row)
    }
  }
  invisible(NULL)
}
