# Explicit adapters keep optional interval classes outside the drawing engine.

ggchord_conversion_require <- function(package, caller) {
  if (!requireNamespace(package, quietly = TRUE)) ggchord_stop(
    caller, ": converting this object requires the optional package '", package,
    "'; install it in your R environment, or supply a data.frame"
  )
}

ggchord_conversion_mapping <- function(mapping, columns, allowed, caller) {
  if (is.null(mapping)) return(setNames(character(), character()))
  if (!is.character(mapping) || is.null(names(mapping)) || anyNA(mapping) ||
      anyNA(names(mapping)) || any(!nzchar(mapping)) || any(!nzchar(names(mapping))) ||
      anyDuplicated(names(mapping)) || any(!names(mapping) %in% allowed)) {
    ggchord_stop(caller, ": columns must be a named character vector from output ",
      "roles to source columns; supported roles: ", paste(allowed, collapse = ", "))
  }
  missing <- setdiff(unname(mapping), columns)
  if (length(missing)) ggchord_stop(caller, ": mapped source column(s) absent: ",
    paste(missing, collapse = ", "))
  mapping
}

ggchord_conversion_append <- function(out, extra, prefix = "metadata_") {
  map <- data.frame(source = character(), output = character(), stringsAsFactors = FALSE)
  if (!ncol(extra)) return(list(data = out, map = map))
  for (i in seq_along(extra)) {
    source <- names(extra)[i]
    target <- source
    if (!nzchar(target) || target %in% c(names(out), "seq_id")) target <- paste0(prefix, source)
    while (target %in% names(out)) target <- paste0(prefix, target)
    out[[target]] <- extra[[i]]
    map <- rbind(map, data.frame(source = source, output = target, stringsAsFactors = FALSE))
  }
  list(data = out, map = map)
}

ggchord_conversion_frame <- function(x, mapping, roles, caller) {
  if (anyDuplicated(names(x)) || anyNA(names(x)) || any(!nzchar(names(x)))) {
    ggchord_stop(caller, ": source columns must have unique non-empty names")
  }
  mapping <- ggchord_conversion_mapping(mapping, names(x), roles, caller)
  selected <- setNames(intersect(roles, names(x)), intersect(roles, names(x)))
  selected[names(mapping)] <- mapping
  out <- data.frame(.source_row = x[[".source_row"]] %||% seq_len(nrow(x)))
  for (role in names(selected)) out[[role]] <- x[[selected[[role]]]]
  # Renamed source roles survive under their original metadata names when an
  # explicit mapping replaces them; one source may intentionally map twice.
  extra <- x[setdiff(names(x), c(unname(selected), ".source_row"))]
  appended <- ggchord_conversion_append(out, extra)
  list(data = appended$data,
    map = rbind(data.frame(source = unname(selected), output = names(selected),
      stringsAsFactors = FALSE), appended$map))
}

ggchord_conversion_finish <- function(out, x, map, metadata = NULL) {
  rownames(out) <- NULL
  sequence_info <- NULL
  if (inherits(x, "GRanges") || inherits(x, "GInteractions") || inherits(x, "Seqinfo")) {
    info_source <- if (inherits(x, "GInteractions")) InteractionSet::regions(x) else x
    sequence_info <- as.data.frame(if (inherits(info_source, "Seqinfo")) info_source else
      GenomicRanges::seqinfo(info_source))
  }
  attr(out, "ggchord_conversion") <- list(
    source_class = class(x), coordinate_system = "1-based-inclusive",
    input_rows = if (is.data.frame(x)) nrow(x) else if (inherits(x, "XString")) 1L else length(x),
    output_rows = nrow(out), column_map = map, source_metadata = metadata,
    source_names = if (is.data.frame(x)) rownames(x) else names(x),
    sequence_info = sequence_info,
    source_import = attr(x, "ggchord_import"),
    source_projection = attr(x, "ggchord_projection"),
    source_conversion = attr(x, "ggchord_conversion")
  )
  out
}

ggchord_conversion_ids <- function(ids, caller) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(ids))) ggchord_stop(
    caller, ": sequence identifiers must be non-missing and non-empty")
  ids
}

ggchord_conversion_coordinates <- function(out, columns, caller) {
  for (column in columns) {
    value <- out[[column]]
    if (!is.numeric(value) || any(!is.finite(value)) || any(value < 1 | value != floor(value))) {
      ggchord_stop(caller, ": ", column, " must contain positive integer ",
        "1-based inclusive coordinates; convert zero-based inputs explicitly")
    }
  }
  invisible(out)
}

ggchord_declare_sequence_unit <- function(data, unit, caller, default = NULL) {
  allowed <- c("bp", "aa")
  if (!is.null(unit) && (!is.character(unit) || length(unit) != 1L || !unit %in% allowed))
    ggchord_stop(caller, ": sequence_unit must be 'bp' or 'aa'")
  existing <- data$sequence_unit
  if (!is.null(existing)) {
    existing <- as.character(existing)
    ggchord_check_declared_unit(existing, unit, caller)
  }
  selected <- unit %||% default
  if (!is.null(selected)) {
    if (is.null(existing)) existing <- rep(NA_character_, nrow(data))
    existing[is.na(existing)] <- selected
  }
  if (!is.null(existing)) data$sequence_unit <- existing
  data
}

ggchord_validate_feature_segments <- function(data, caller) {
  if (!is.list(data$segments)) ggchord_stop(caller, ": segments must be a list of interval tables")
  for (i in seq_len(nrow(data))) {
    parts <- data$segments[[i]]
    if (is.null(parts)) next
    if (!is.data.frame(parts)) ggchord_stop(caller, ": segments must be interval tables, row ", i)
    ggchord_require_columns(parts, c("start", "end"), caller)
    ggchord_conversion_coordinates(parts, c("start", "end"), caller)
    if (any(parts$start > parts$end | parts$start < data$start[i] | parts$end > data$end[i]))
      ggchord_stop(caller, ": segments must be ordered within their parent interval, row ", i)
  }
}

ggchord_link_endpoint_strand <- function(start, end, strand, caller, field,
                                          directional_order = FALSE) {
  direction <- as.character(strand %||% rep(NA_character_, length(start)))
  if (any(!is.na(direction) & !direction %in% c("+", "-", ".", "*", "?")))
    ggchord_stop(caller, ": invalid ", field)
  descending <- is.finite(start) & is.finite(end) & start > end
  ascending <- is.finite(start) & is.finite(end) & start < end
  if (any(descending & direction %in% "+") ||
      (directional_order && any(ascending & direction %in% "-")))
    ggchord_stop(caller, ": endpoint order conflicts with ", field)
  implicit <- is.na(direction)
  direction[descending & implicit] <- "-"
  direction[ascending & implicit & directional_order] <- "+"
  direction[!direction %in% c("+", "-")] <- "."
  direction
}

ggchord_order_link_endpoints <- function(data, caller, directional_order = FALSE) {
  for (prefix in c("q", "s")) {
    start <- paste0(prefix, "start")
    end <- paste0(prefix, "end")
    strand <- paste0(prefix, "strand")
    before <- data[[start]]
    after <- data[[end]]
    if (!nrow(data)) {
      data[[strand]] <- character()
      next
    }
    direction <- ggchord_link_endpoint_strand(before, after, data[[strand]], caller, strand,
      directional_order)
    for (field in c(start, end)) {
      original <- paste0("source_", field)
      if (!original %in% names(data)) data[[original]] <- data[[field]]
    }
    data[[start]] <- pmin(before, after)
    data[[end]] <- pmax(before, after)
    data[[strand]] <- direction
  }
  data
}

ggchord_conversion_ranges <- function(x, caller) {
  ggchord_conversion_require("GenomicRanges", caller)
  ggchord_conversion_require("S4Vectors", caller)
  out <- data.frame(
    accver = as.character(GenomicRanges::seqnames(x)),
    start = as.numeric(GenomicRanges::start(x)),
    end = as.numeric(GenomicRanges::end(x)),
    strand = as.character(GenomicRanges::strand(x)),
    .source_row = seq_along(x), stringsAsFactors = FALSE
  )
  if (any(out$end < out$start)) ggchord_stop(
    caller, ": zero-width source intervals cannot be drawn as genomic features or links")
  out$directionality <- unname(c("+" = "forward", "-" = "reverse", "*" = "none")[out$strand])
  if (!is.null(names(x))) out$.source_name <- names(x)
  extras <- as.data.frame(S4Vectors::mcols(x), stringsAsFactors = FALSE)
  appended <- ggchord_conversion_append(out, extras)
  list(data = appended$data, map = appended$map,
    metadata = S4Vectors::metadata(x))
}

ggchord_conversion_seqinfo <- function(x, caller) {
  if (!inherits(x, "Seqinfo")) ggchord_conversion_require("GenomicRanges", caller)
  info <- if (inherits(x, "Seqinfo")) x else GenomicRanges::seqinfo(x)
  # as.data.frame(Seqinfo) is the public class method and works with both the
  # current standalone Seqinfo package and earlier GenomeInfoDb ownership.
  info <- as.data.frame(info)
  if (!"seqlengths" %in% names(info)) ggchord_stop(
    caller, ": sequence information does not expose seqlengths")
  sequence_names <- if ("seqnames" %in% names(info)) as.character(info$seqnames) else rownames(info)
  out <- data.frame(accver = sequence_names, length = as.numeric(info$seqlengths),
    stringsAsFactors = FALSE)
  if ("isCircular" %in% names(info)) out$circular <- info$isCircular
  if ("genome" %in% names(info)) out$genome <- info$genome
  out
}

#' Adapt a complete reference dictionary
#'
#' A reference dictionary has one row per stable reference ID, with a complete
#' positive integral length and an explicit native coordinate unit. Optional
#' sequence and feature range classes supply their declared whole-reference
#' lengths. Observation endpoints cannot establish a full sequence length.
#'
#' @param x A reference data frame, named length vector, or supported
#'   `Seqinfo`, `GRanges`, `GInteractions` or Biostrings sequence object.
#' @param columns Named mapping from canonical roles to data-frame fields,
#'   such as `c(accver = "chrom", length = "size")`.
#' @param lengths Named positive integral lengths. They may fill missing
#'   dictionary values and must agree with known lengths.
#' @param ids Explicit reference IDs for unnamed input. Established source
#'   IDs retain their original identity.
#' @param sequence_unit Native `"bp"` or `"aa"`. The input declaration takes
#'   precedence over the documented default of `"bp"`. Explicit declarations
#'   must agree.
#' @param topology `"linear"`, `"circular"`, or a named configuration by
#'   reference ID. Absent source topology defaults to linear.
#' @return An ordinary data frame with `accver`, `length`, `sequence_unit`,
#'   and `topology`. [data_report()] records source metadata and defaults.
#' @export
#' @examples
#' raw <- data.frame(chrom = c("A", "B"), size = c(1000, 800))
#' ref <- as_ggchord_seq(raw, columns = c(accver = "chrom", length = "size"),
#'                        sequence_unit = "bp")
#' ggchord(ref) + geom_seq()
as_ggchord_seq <- function(x, columns = NULL, lengths = NULL, ids = NULL,
                           sequence_unit = NULL, topology = NULL) {
  caller <- "as_ggchord_seq()"
  map <- data.frame(source = character(), output = character(), stringsAsFactors = FALSE)
  metadata <- NULL
  sequence_object <- inherits(x, "XStringSet") || inherits(x, "XString")
  if (is.numeric(x) && is.null(dim(x))) {
    if (!is.null(columns)) ggchord_stop(caller, ": columns are only for data-frame inputs")
    vector_ids <- names(x) %||% ids
    if (is.null(vector_ids) || length(vector_ids) != length(x))
      ggchord_stop(caller, ": provide one reference ID per length")
    x <- data.frame(accver = vector_ids, length = unname(x))
  }
  if (is.data.frame(x)) {
    converted <- ggchord_conversion_frame(x, columns, c("accver", "length", "sequence_unit",
      "topology", "circular", "source_accver", "source_origin", "source_length"), caller)
    out <- converted$data; map <- converted$map
    if (!"accver" %in% names(out) && !is.null(ids) && length(ids) == nrow(out)) out$accver <- ids
    ggchord_require_columns(out, "accver", caller)
    if (!"length" %in% names(out)) out$length <- rep(NA_real_, nrow(out))
  } else {
    if (!is.null(columns)) ggchord_stop(caller, ": columns are for data-frame inputs")
    if (sequence_object) {
      ggchord_conversion_require("Biostrings", caller)
      count <- if (inherits(x, "XString")) 1L else length(x)
      if (!count) ggchord_stop(caller, ": sequence input must contain at least one sequence")
      source_ids <- if (inherits(x, "XString")) NULL else names(x)
      if (is.null(ids)) ids <- source_ids
      if (is.null(ids) || !is.character(ids) || length(ids) != count ||
          anyNA(ids) || any(!nzchar(ids))) ggchord_stop(caller,
        ": unnamed sequence objects require one explicit nonempty ids value per sequence")
      if (!is.null(source_ids)) {
        known <- !is.na(source_ids) & nzchar(source_ids)
        if (any(source_ids[known] != ids[known])) ggchord_stop(caller,
          ": ids disagree with existing source sequence names; identifiers are not renamed")
      }
      widths <- if (inherits(x, "XString")) length(x) else Biostrings::width(x)
      out <- data.frame(accver = ids, length = as.numeric(widths),
        sequence_unit = if (inherits(x, "AAString") || inherits(x, "AAStringSet")) "aa" else
          if (inherits(x, "DNAString") || inherits(x, "DNAStringSet") ||
              inherits(x, "RNAString") || inherits(x, "RNAStringSet")) "bp" else "character",
        .source_row = seq_len(count), stringsAsFactors = FALSE)
      if (inherits(x, "Annotated")) metadata <- S4Vectors::metadata(x)
      if (inherits(x, "XStringSet")) {
        appended <- ggchord_conversion_append(out,
          as.data.frame(S4Vectors::mcols(x), stringsAsFactors = FALSE))
        out <- appended$data; map <- appended$map
      }
    } else if (inherits(x, "GInteractions")) {
      ggchord_conversion_require("InteractionSet", caller)
      first <- InteractionSet::anchors(x, type = "first", id = FALSE)
      second <- InteractionSet::anchors(x, type = "second", id = FALSE)
      used <- unique(c(as.character(GenomicRanges::seqnames(first)),
        as.character(GenomicRanges::seqnames(second))))
      out <- ggchord_conversion_seqinfo(InteractionSet::regions(x), caller)
      out <- out[match(used, out$accver), , drop = FALSE]
      metadata <- S4Vectors::metadata(x)
    } else if (inherits(x, "GRanges") || inherits(x, "Seqinfo")) {
      out <- ggchord_conversion_seqinfo(x, caller)
      if (inherits(x, "GRanges")) {
        used <- unique(as.character(GenomicRanges::seqnames(x)))
        out <- out[match(used, out$accver), , drop = FALSE]
      }
      metadata <- if (inherits(x, "GRanges")) S4Vectors::metadata(x) else NULL
    } else ggchord_stop(caller, ": x must be lengths, a data frame, or a supported sequence object")
  }
  if (!nrow(out)) ggchord_stop(caller, ": at least one reference is required")
  if (!is.null(ids)) {
    if (!is.character(ids) || length(ids) != nrow(out) || anyNA(ids) || any(!nzchar(ids)))
      ggchord_stop(caller, ": ids must provide one nonempty ID per reference")
    known <- !is.na(out$accver) & nzchar(out$accver)
    if (any(out$accver[known] != ids[known])) ggchord_stop(caller, ": ids cannot replace source IDs")
    out$accver <- ids
  }
  out$accver <- ggchord_conversion_ids(out$accver, caller)
  if (anyDuplicated(out$accver)) ggchord_stop(caller, ": every sequence ID must occur once")
  if (!is.numeric(out$length)) ggchord_stop(caller, ": sequence lengths must be numeric")
  if (!is.null(lengths)) {
    if (!is.numeric(lengths) || is.null(names(lengths)) || anyNA(names(lengths)) ||
        any(!nzchar(names(lengths))) || anyDuplicated(names(lengths)) ||
        any(!is.finite(lengths)) || any(lengths < 1 | lengths != floor(lengths)) ||
        any(!names(lengths) %in% out$accver)) {
      ggchord_stop(caller, ": lengths must be a named vector of positive integer lengths")
    }
    supplied <- lengths[match(out$accver, names(lengths))]
    conflict <- is.finite(out$length) & !is.na(supplied) & out$length != supplied
    if (any(conflict)) ggchord_stop(caller, ": explicit lengths disagree with source lengths for: ",
      paste(out$accver[conflict], collapse = ", "))
    fill <- is.na(out$length) & !is.na(supplied)
    out$length[fill] <- supplied[fill]
  }
  bad <- !is.finite(out$length) | out$length < 1 | out$length > 2^53 | out$length != floor(out$length)
  if (any(bad)) {
    ggchord_stop(caller, ": unknown or invalid full sequence lengths for: ",
      paste(out$accver[bad], collapse = ", "),
      "; set source seqlengths or supply lengths = c(accver = ...); ",
      "feature endpoints are not full sequence lengths")
  }
  default_unit <- if (is.null(sequence_unit)) {
    if (is.null(out$sequence_unit)) out$accver else out$accver[is.na(out$sequence_unit)]
  } else character()
  out <- ggchord_declare_sequence_unit(out, sequence_unit, caller, default = "bp")
  known_topology <- out$topology
  if (is.null(known_topology)) known_topology <- rep(NA_character_, nrow(out))
  if (!is.null(out$circular)) {
    if (!is.logical(out$circular)) ggchord_stop(caller, ": circular must be logical")
    circular_topology <- ifelse(out$circular, "circular", "linear")
    conflict <- !is.na(known_topology) & !is.na(circular_topology) & known_topology != circular_topology
    if (any(conflict)) ggchord_stop(caller, ": circular and topology disagree")
    known_topology[is.na(known_topology)] <- circular_topology[is.na(known_topology)]
  }
  if (!is.null(topology)) {
    if (!is.character(topology) || any(!topology %in% ggchord_sequence_topologies))
      ggchord_stop(caller, ": topology must be 'linear' or 'circular'")
    if (!is.null(names(topology))) {
      if (anyDuplicated(names(topology)) || any(!names(topology) %in% out$accver))
        ggchord_stop(caller, ": topology keys must be unique reference IDs")
      topology <- unname(topology[out$accver])
    } else if (length(topology) != 1L) {
      ggchord_stop(caller, ": use a scalar topology or values named by reference ID")
    }
    topology <- rep_len(topology, nrow(out))
    if (any(!is.na(known_topology) & !is.na(topology) & known_topology != topology))
      ggchord_stop(caller, ": topology disagrees with the source")
    known_topology[is.na(known_topology)] <- topology[is.na(known_topology)]
  }
  default_topology <- is.na(known_topology)
  known_topology[default_topology] <- "linear"
  if (any(!known_topology %in% ggchord_sequence_topologies))
    ggchord_stop(caller, ": topology must be 'linear' or 'circular'")
  out$topology <- as.character(known_topology)
  out$circular <- out$topology == "circular"
  out <- ggchord_conversion_finish(out, x, map, metadata)
  attr(out, "ggchord_conversion")$defaults <- list(sequence_unit = if (length(default_unit))
    setNames(rep("bp", length(default_unit)), default_unit),
    topology = out$accver[default_topology])
  out
}

#' Adapt one-based sequence feature intervals
#'
#' Coordinates are ordered and closed. A missing strand is reported as `"."`.
#' A `segments` list of explicit ordered interval tables can record
#' discontinuous features without filling gaps in their parent range.
#'
#' @param x A data frame or `GRanges` object.
#' @param columns Named mapping from canonical roles to input column names.
#' @param sequence_unit Optional native `"bp"` or `"aa"` declaration. An
#'   existing declaration must agree. Undeclared observation units can be
#'   resolved against a reference later.
#' @return A feature data frame with `accver`, `start`, `end`, strand
#'   and available feature metadata. Use [data_report()] for provenance.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' raw <- data.frame(chrom = "A", begin = 100, finish = 240, direction = "+")
#' features <- as_ggchord_feature(raw, columns = c(
#'   accver = "chrom", start = "begin", end = "finish", strand = "direction"))
#' ggchord(seqs) + geom_feature(data = features)
as_ggchord_feature <- function(x, columns = NULL, sequence_unit = NULL) {
  caller <- "as_ggchord_feature()"
  out <- ggchord_convert_feature(x, columns)
  if (any(out$start > out$end)) ggchord_stop(caller,
    ": start must not exceed end; use segments for discontinuous or circular-origin intervals")
  if ("segments" %in% names(out)) ggchord_validate_feature_segments(out, caller)
  ggchord_declare_sequence_unit(out, sequence_unit, caller)
}

ggchord_convert_feature <- function(x, mapping, unified = FALSE) {
  caller <- "as_ggchord_feature()"
  metadata <- NULL
  if (is.data.frame(x)) {
    roles <- c("accver", "start", "end", "strand", "directionality", "feature_id", "parent_ids",
      "type", "label", "anno", "feature_color", "display_priority", "segments", "sequence_unit")
    converted <- ggchord_conversion_frame(x, mapping, roles, caller)
  } else if (inherits(x, "GRanges")) {
    if (!is.null(mapping)) ggchord_stop(caller, ": mapping is for data-frame inputs")
    converted <- ggchord_conversion_ranges(x, caller)
    metadata <- converted$metadata
  } else ggchord_stop(caller, ": x must be a data.frame or GRanges")
  out <- converted$data
  ggchord_require_columns(out, c("accver", "start", "end"), caller)
  out$accver <- ggchord_conversion_ids(out$accver, caller)
  ggchord_conversion_coordinates(out, c("start", "end"), caller)
  if (unified) out <- ggchord_normalize_feature_input(out)
  strand <- if ("strand" %in% names(out)) as.character(out$strand) else NULL
  direction_roles <- if (unified) c("directionality", "direction") else "directionality"
  direction_column <- intersect(direction_roles, names(out))[1L]
  expected <- if (!is.na(direction_column)) ggchord_direction_strand(out[[direction_column]]) else NULL
  if (anyNA(expected)) {
    ggchord_stop(caller, ": directionality must be forward, reverse, bidirectional or none")
  }
  if (is.null(strand)) strand <- expected %||% rep(".", nrow(out))
  if (any(!is.na(strand) & !strand %in% ggchord_feature_strands)) {
    ggchord_stop(caller, ": strand must be '+', '-', '*', '.', or '+/-'")
  }
  source_strand <- out$source_strand %||% strand
  strand[is.na(strand) | strand == "*"] <- "."
  if (!is.null(expected) && any(strand != expected)) {
    ggchord_stop(caller, ": strand and directionality describe conflicting directions")
  }
  out$strand <- strand
  out$directionality <- unname(c("+" = "forward", "-" = "reverse", "." = "none", "+/-" = "bidirectional")[strand])
  map <- converted$map
  out$source_strand <- source_strand
  map <- rbind(map, data.frame(source = "strand", output = "source_strand",
    stringsAsFactors = FALSE))
  ggchord_conversion_finish(out, x, map, metadata)
}

#' Adapt point or interval relationships
#'
#' Query and subject roles and biological direction evidence are kept
#' separate. Interval endpoints are ordered within each reference, and
#' original endpoint order is retained as source metadata. Equal endpoints
#' carry no direction evidence of their own.
#'
#' @param x Data frame, or a `GInteractions` object for interval links.
#' @param columns Named mapping from canonical roles to input fields.
#' @param kind `"interval"` consumes qstart/qend and sstart/send.
#'   `"point"` consumes qpos and spos. Both require qaccver and saccver.
#' @param sequence_unit Optional declaration of native `"bp"` or `"aa"`.
#'   This operation preserves the scientific meaning of existing units.
#' @return A point or interval relationship table and conversion provenance.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' paired <- as_ggchord_link(links, kind = "interval", sequence_unit = "bp")
#' ggchord(seqs) + geom_link_ribbon(data = paired)
as_ggchord_link <- function(x, columns = NULL, kind = "interval", sequence_unit = NULL) {
  caller <- "as_ggchord_link()"
  kind <- ggchord_match_arg(kind, c("interval", "point"))
  endpoints <- if (kind == "interval") c("qstart", "qend", "sstart", "send") else c("qpos", "spos")
  required <- c("qaccver", "saccver", endpoints)
  metadata <- NULL
  if (is.data.frame(x)) {
    converted <- ggchord_conversion_frame(x, columns,
      c(required, "qstrand", "sstrand", "pident", "length", "direction", "orientation_known",
        "link_id", "qlen", "slen", "sequence_unit", "qsegments", "ssegments"), caller)
    out <- converted$data; map <- converted$map
  } else if (inherits(x, "GInteractions")) {
    if (kind == "point") ggchord_stop(caller, ": provide explicit point endpoints for kind = 'point'")
    if (!is.null(columns)) ggchord_stop(caller, ": columns are for data-frame inputs")
    ggchord_conversion_require("InteractionSet", caller)
    first <- InteractionSet::anchors(x, type = "first", id = FALSE)
    second <- InteractionSet::anchors(x, type = "second", id = FALSE)
    query <- ggchord_conversion_ranges(first, caller)$data
    subject <- ggchord_conversion_ranges(second, caller)$data
    out <- data.frame(qaccver = query$accver, saccver = subject$accver,
      qstart = ifelse(query$strand == "-", query$end, query$start),
      qend = ifelse(query$strand == "-", query$start, query$end),
      sstart = ifelse(subject$strand == "-", subject$end, subject$start),
      send = ifelse(subject$strand == "-", subject$start, subject$end),
      qstrand = query$strand, sstrand = subject$strand,
      .source_row = seq_along(x),
      .query_anchor = InteractionSet::anchors(x, type = "first", id = TRUE),
      .subject_anchor = InteractionSet::anchors(x, type = "second", id = TRUE),
      stringsAsFactors = FALSE)
    out$orientation_known <- out$qstrand %in% c("+", "-") & out$sstrand %in% c("+", "-")
    out$direction <- ifelse(out$orientation_known,
      ifelse(out$qstrand == out$sstrand, "same", "reverse"), NA_character_)
    if (!is.null(names(x))) out$.source_name <- names(x)
    map <- data.frame(source = character(), output = character(), stringsAsFactors = FALSE)
    anchor_roles <- c("accver", "start", "end", "strand", ".source_row")
    for (pair in list(list(data = query, prefix = "query_"),
        list(data = subject, prefix = "subject_"))) {
      extra <- pair$data[setdiff(names(pair$data), anchor_roles)]
      names(extra) <- paste0(pair$prefix, names(extra))
      appended <- ggchord_conversion_append(out, extra)
      out <- appended$data; map <- rbind(map, appended$map)
    }
    appended <- ggchord_conversion_append(out,
      as.data.frame(S4Vectors::mcols(x), stringsAsFactors = FALSE))
    out <- appended$data; map <- rbind(map, appended$map)
    metadata <- S4Vectors::metadata(x)
  } else ggchord_stop(caller, ": x must be a data.frame or GInteractions")
  ggchord_require_columns(out, required, caller)
  out$qaccver <- ggchord_conversion_ids(out$qaccver, caller)
  out$saccver <- ggchord_conversion_ids(out$saccver, caller)
  ggchord_conversion_coordinates(out, endpoints, caller)
  if (kind == "interval") out <- ggchord_order_link_endpoints(out, caller)
  ggchord_link_orientation(out, caller)
  out <- ggchord_declare_sequence_unit(out, sequence_unit, caller)
  ggchord_conversion_finish(out, x, map, metadata)
}
