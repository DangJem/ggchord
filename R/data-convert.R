# Explicit adapters keep optional interval classes outside the drawing engine.

ggchord_conversion_require <- function(package, caller) {
  if (!requireNamespace(package, quietly = TRUE)) ggchord_stop(
    caller, ": converting this object requires the optional package '", package,
    "'; install it in your R environment, or supply a data.frame"
  )
}

ggchord_conversion_mapping <- function(mapping, columns, allowed, caller) {
  if (is.null(mapping)) return(stats::setNames(character(), character()))
  if (!is.character(mapping) || is.null(names(mapping)) || anyNA(mapping) ||
      anyNA(names(mapping)) || any(!nzchar(mapping)) || any(!nzchar(names(mapping))) ||
      anyDuplicated(names(mapping)) || any(!names(mapping) %in% allowed)) {
    ggchord_stop(caller, ": mapping must be a named character vector from output ",
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
  selected <- stats::setNames(intersect(roles, names(x)), intersect(roles, names(x)))
  selected[names(mapping)] <- mapping
  out <- data.frame(.source_row = seq_len(nrow(x)))
  for (role in names(selected)) out[[role]] <- x[[selected[[role]]]]
  # Renamed source roles survive under their original metadata names when an
  # explicit mapping replaces them; one source may intentionally map twice.
  extra <- x[setdiff(names(x), unique(unname(selected)))]
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
    input_rows = if (is.data.frame(x)) nrow(x) else length(x),
    output_rows = nrow(out), column_map = map, source_metadata = metadata,
    source_names = if (is.data.frame(x)) rownames(x) else names(x),
    sequence_info = sequence_info
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

#' Convert sequence information to ggchord input
#'
#' Explicit adapters return ordinary data frames. Optional interval packages
#' are loaded only when their objects are supplied; drawing never installs or
#' runs external data preparation. Supported S4 inputs supplied directly to
#' ggchord or its layers use the same adapters. Identifiers are kept verbatim.
#'
#' For a `GRanges` object, lengths come from its `Seqinfo`, never from the
#' maximum observed feature endpoint. A `GInteractions` object uses the same
#' information from its common regions. Only sequence names used by the
#' supplied intervals are returned; a `Seqinfo` object returns all its entries.
#' Unknown lengths require explicit `lengths` or corrected source information.
#'
#' @param x A data frame, `GRanges`, `Seqinfo`, or `GInteractions` object.
#' @param mapping Optional named character vector `output_role = "input_column"`
#'   for data-frame inputs. For example `c(accver = "chrom", length = "size")`.
#' @param lengths Optional named numeric vector of known sequence lengths.
#'   Fills unknown lengths and must agree with any known source lengths.
#' @return An ordinary data frame with `accver`, `length` and retained metadata.
#'   Its `ggchord_conversion` attribute records input class, coordinate system,
#'   column renaming and object-level metadata. No genomic coordinates or
#'   accession versions are guessed or rewritten.
#' @md
#' @export
#' @examples
#' as_ggchord_seq(data.frame(chrom = c("host", "virus"), size = c(1000, 300)),
#'   mapping = c(accver = "chrom", length = "size"))
as_ggchord_seq <- function(x, mapping = NULL, lengths = NULL) {
  caller <- "as_ggchord_seq()"
  map <- data.frame(source = character(), output = character(), stringsAsFactors = FALSE)
  metadata <- NULL
  if (is.data.frame(x)) {
    converted <- ggchord_conversion_frame(x, mapping, c("accver", "length"), caller)
    out <- converted$data; map <- converted$map
    ggchord_require_columns(out, "accver", caller)
    if (!"length" %in% names(out)) out$length <- rep(NA_real_, nrow(out))
  } else {
    if (!is.null(mapping)) ggchord_stop(caller, ": mapping is for data-frame inputs")
    if (inherits(x, "GInteractions")) {
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
    } else ggchord_stop(caller, ": x must be a data.frame, GRanges, Seqinfo or GInteractions")
  }
  out$accver <- ggchord_conversion_ids(out$accver, caller)
  if (anyDuplicated(out$accver)) ggchord_stop(caller, ": every sequence ID must occur once")
  if (!is.numeric(out$length)) ggchord_stop(caller, ": sequence lengths must be numeric")
  if (!is.null(lengths)) {
    if (!is.numeric(lengths) || is.null(names(lengths)) || anyNA(names(lengths)) ||
        any(!nzchar(names(lengths))) || anyDuplicated(names(lengths)) ||
        any(!is.finite(lengths)) || any(lengths < 1 | lengths != floor(lengths))) {
      ggchord_stop(caller, ": lengths must be a named vector of positive integer lengths")
    }
    supplied <- lengths[match(out$accver, names(lengths))]
    conflict <- is.finite(out$length) & !is.na(supplied) & out$length != supplied
    if (any(conflict)) ggchord_stop(caller, ": explicit lengths disagree with source lengths for: ",
      paste(out$accver[conflict], collapse = ", "))
    fill <- is.na(out$length) & !is.na(supplied)
    out$length[fill] <- supplied[fill]
  }
  if (any(!is.finite(out$length)) || any(out$length < 1 | out$length != floor(out$length))) {
    bad <- !is.finite(out$length) | out$length < 1 | out$length != floor(out$length)
    ggchord_stop(caller, ": unknown or invalid full sequence lengths for: ",
      paste(out$accver[bad], collapse = ", "),
      "; set source seqlengths or supply lengths = c(accver = ...); ",
      "feature endpoints are not full sequence lengths")
  }
  ggchord_conversion_finish(out, x, map, metadata)
}

#' Convert intervals to generic ggchord features
#'
#' `GRanges` coordinates are already one-based inclusive and are not shifted.
#' Strand `*` is retained as `source_strand`, while the output uses
#' `strand = "."` and `directionality = "none"`. Use the returned table in
#' `geom_feature(data = ...)` to preserve nondirectional intervals; `geom_gene()`
#' requires directional gene input. No unknown strand is invented as `+`.
#'
#' @param x A data frame or `GRanges` object.
#' @param mapping Optional named character vector mapping output roles to
#'   source columns, as in [as_ggchord_seq()].
#' @return An ordinary feature data frame with `accver`, `start`, `end`,
#'   `strand`, `directionality`, retained metadata and conversion provenance.
#'   Conflicting metadata names receive a `metadata_` prefix; the attribute's
#'   `column_map` records the output name.
#' @md
#' @export
#' @examples
#' features <- as_ggchord_feature(data.frame(chrom = "host", from = 100,
#'   to = 200, orientation = "*", label = "integration site"),
#'   mapping = c(accver = "chrom", start = "from", end = "to",
#'     strand = "orientation"))
#' ggchord(data.frame(accver = "host", length = 1000)) + geom_seq() +
#'   geom_feature(data = features) + coord_collinear()
as_ggchord_feature <- function(x, mapping = NULL) {
  caller <- "as_ggchord_feature()"
  metadata <- NULL
  if (is.data.frame(x)) {
    roles <- c("accver", "start", "end", "strand", "directionality", "feature_id",
      "type", "label", "anno", "feature_color", "display_priority")
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
  strand <- if ("strand" %in% names(out)) as.character(out$strand) else NULL
  direction <- if ("directionality" %in% names(out)) as.character(out$directionality) else NULL
  direction_map <- c(forward = "+", reverse = "-", bidirectional = "+/-", none = ".", nondirectional = ".")
  if (!is.null(direction) && (anyNA(direction) || any(!direction %in% names(direction_map)))) {
    ggchord_stop(caller, ": directionality must be forward, reverse, bidirectional or none")
  }
  if (is.null(strand)) strand <- if (is.null(direction)) rep(".", nrow(out)) else unname(direction_map[direction])
  if (anyNA(strand) || any(!strand %in% c("+", "-", "*", ".", "+/-"))) {
    ggchord_stop(caller, ": strand must be '+', '-', '*', '.', or '+/-'")
  }
  source_strand <- strand
  strand[strand == "*"] <- "."
  if (!is.null(direction) && any(strand != unname(direction_map[direction]))) {
    ggchord_stop(caller, ": strand and directionality describe conflicting directions")
  }
  out$strand <- strand
  out$directionality <- unname(c("+" = "forward", "-" = "reverse", "." = "none", "+/-" = "bidirectional")[strand])
  map <- converted$map
  if ("source_strand" %in% names(out) && !identical(out$source_strand, source_strand)) {
    appended <- ggchord_conversion_append(out, out["source_strand"])
    out <- appended$data; map <- rbind(map, appended$map)
  }
  out$source_strand <- source_strand
  map <- rbind(map, data.frame(source = "strand", output = "source_strand",
    stringsAsFactors = FALSE))
  ggchord_conversion_finish(out, x, map, metadata)
}

#' Convert paired genomic intervals to ggchord links
#'
#' For `GInteractions`, the first stored anchor becomes query and the second
#' becomes subject. Their order is never sorted or swapped by this adapter.
#' A source object that already canonicalized its anchors cannot recover their
#' earlier biological host/virus roles; establish those roles before conversion.
#' Negative-strand anchor endpoints run from end to start. Unknown-strand
#' anchors retain increasing coordinates, their `*` strand and
#' `orientation_known = FALSE`; increasing order is not direction evidence.
#' Interval interactions do not imply a pairwise alignment: identity and
#' alignment length are not inferred from anchor widths.
#'
#' @param x A data frame or `GInteractions` object.
#' @param mapping Optional named character vector mapping canonical link roles
#'   to source columns, as in [as_ggchord_seq()].
#' @return An ordinary data frame with `qaccver`, `saccver`, `qstart`, `qend`,
#'   `sstart`, `send` and source metadata. Paired objects also retain anchor
#'   indices, source row/name, anchor strands, anchor metadata prefixed with
#'   `query_`/`subject_`, and object-level metadata in the conversion attribute.
#' @md
#' @export
#' @examples
#' as_ggchord_link(data.frame(host = "host", virus = "virus",
#'   host_start = 100, host_end = 120, virus_start = 10, virus_end = 30),
#'   mapping = c(qaccver = "host", saccver = "virus", qstart = "host_start",
#'     qend = "host_end", sstart = "virus_start", send = "virus_end"))
as_ggchord_link <- function(x, mapping = NULL) {
  caller <- "as_ggchord_link()"
  required <- c("qaccver", "saccver", "qstart", "qend", "sstart", "send")
  metadata <- NULL
  if (is.data.frame(x)) {
    converted <- ggchord_conversion_frame(x, mapping,
      c(required, "qstrand", "sstrand", "pident", "length", "direction", "orientation_known"), caller)
    out <- converted$data; map <- converted$map
  } else if (inherits(x, "GInteractions")) {
    if (!is.null(mapping)) ggchord_stop(caller, ": mapping is for data-frame inputs")
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
  ggchord_conversion_coordinates(out, c("qstart", "qend", "sstart", "send"), caller)
  ggchord_conversion_finish(out, x, map, metadata)
}
