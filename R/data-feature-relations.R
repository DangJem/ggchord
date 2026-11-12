# Offline precomputed feature relationships; never infer pairwise edges from groups.

ggchord_relation_features <- function(data, caller) {
  if (!is.data.frame(data) || !all(c('feature_id', 'accver', 'start', 'end') %in% names(data)))
    ggchord_stop(caller, ': feature_data needs feature_id, accver, start and end')
  for (key in c('feature_id', 'accver')) {
    if (!is.character(data[[key]]) || anyNA(data[[key]]) || any(!nzchar(data[[key]])))
      ggchord_stop(caller, ': feature_data ', key, ' must contain nonempty character IDs')
  }
  if (anyDuplicated(data$feature_id)) ggchord_stop(caller, ': ambiguous duplicate feature_id in feature_data')
  for (key in c('start', 'end')) {
    x <- data[[key]]
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x) | x < 1 | x > 2^53-1 | x != floor(x)))
      ggchord_stop(caller, ': feature_data ', key, ' must contain exact positive source coordinates')
  }
  if (any(data$end < data$start)) ggchord_stop(caller, ': feature_data needs ascending 1-based-inclusive intervals')
  if ('segments' %in% names(data) && any(lengths(data$segments) > 0L))
    ggchord_stop(caller, ': segmented features need an explicit contiguous anchor table; bounding spans would fill gaps')
  if ('strand' %in% names(data) &&
      (!is.character(data$strand) || anyNA(data$strand) || any(!data$strand %in% c('+', '-', '.'))))
    ggchord_stop(caller, ': feature_data strand must be +, - or .')
  data
}

#' Read precomputed feature pairs or group memberships
#'
#' Reads local whitespace-separated pairs (two IDs and an optional opaque
#' score), named groups (group ID followed by members), or unnamed groups
#' (members only). Schemas are explicit; headers are not guessed. Blank lines
#' and whole-line # comments are skipped. IDs, native score strings, physical
#' lines, member order, duplicates and original fields remain traceable.
#' Numeric-looking IDs remain strings. Group indices are local to each file;
#' unnamed groups have no invented biological group ID. Membership never
#' generates pairwise links or asserts orthology.
#'
#' Optional feature_data matches unique feature_id values exactly, without
#' editing case, prefixes or versions or guessing gene/CDS parentage. Supply
#' ascending 1-based-inclusive contiguous source intervals, with optional
#' strand (+, - or .). Negative feature strand stays annotation metadata:
#' paired anchor endpoints remain ascending and direction stays missing.
#' These are feature relationships, not basewise alignments. Segmented
#' features require an explicit contiguous anchor table rather than silently
#' using their bounding spans. No full sequence length, identity, support,
#' alignment direction or source characters are inferred.
#'
#' @param file Optional local paths or patterns, read before files.
#' @param files Optional additional local paths or patterns.
#' @param format Explicit schema: "pairs" (default), "named_groups", or "groups".
#' @param feature_data Optional data frame with unique character feature_id,
#'   character accver and numeric start/end in declared source coordinates.
#'   Original columns remain in report$feature_data; optional strand is retained.
#' @param unmatched "error" (default) rejects IDs absent from supplied
#'   feature_data; "report" retains incomplete rows with mapping flags and
#'   reports each missing ID. Without feature_data all rows stay unplaced.
#' @param source_file Add .source_file to pairs, membership and mapped links.
#' @return A list with pairs, membership, ribbon_data (only complete mapped
#'   original pairs), and report. Group input returns no pair links. Mapped
#'   members carry source intervals in membership. report retains raw files,
#'   original fields, missing IDs, unplaced pairs and supplied feature_data.
#'   No seq_data is generated from annotation endpoints. File-local record
#'   indices must be interpreted with .source_input, not as biological IDs.
#' @seealso [geom_feature()], [geom_link_line()], [read_gff3()]
#' @export
#' @examples
#' path <- system.file("extdata", "plant-feature-groups.txt", package = "ggchord")
#' relations <- read_feature_relations(path, format = "named_groups")
#' relations$membership[c("group_id", "feature_id")]
read_feature_relations <- function(file = NULL, files = NULL,
    format = c('pairs', 'named_groups', 'groups'), feature_data = NULL,
    unmatched = c('error', 'report'), source_file = FALSE) {
  caller <- 'read_feature_relations()'; format <- match.arg(format); unmatched <- match.arg(unmatched)
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, files, 'read_feature_relations')
  if (!is.null(feature_data)) feature_data <- ggchord_relation_features(feature_data, caller)
  pairs <- list(); members <- list(); raw <- list(); records <- list()
  for (input in seq_along(paths)) {
    lines <- ggchord_read_text(paths[input], caller); raw[[input]] <- lines
    kept <- which(nzchar(trimws(lines)) & !grepl('^[[:space:]]*#', lines))
    fields <- strsplit(trimws(lines[kept]), '[[:space:]]+')
    records[[input]] <- list(file = paths[input], lines = kept, fields = fields)
    for (row in seq_along(fields)) {
      x <- fields[[row]]
      provenance <- data.frame(.source_input = input, .source_record = row, .source_line = kept[row])
      if (format == 'pairs') {
        if (!length(x) %in% c(2L, 3L)) ggchord_stop(caller, ': pairs need 2 or 3 fields in ', paths[input], ' at line ', kept[row])
        out <- cbind(data.frame(qfeature_id = x[1L], sfeature_id = x[2L],
          source_score = if (length(x) == 3L) x[3L] else NA_character_), provenance)
        pairs[[length(pairs)+1L]] <- ggchord_add_source_file(out, paths[input], source_file)
      } else {
        if (format == 'named_groups' && length(x) < 2L)
          ggchord_stop(caller, ': named_groups need a group ID and at least one member at line ', kept[row])
        group_id <- if (format == 'named_groups') x[1L] else NA_character_
        ids <- if (format == 'named_groups') x[-1L] else x
        out <- cbind(data.frame(group_id = group_id, feature_id = ids,
          .source_group = row, .source_member = seq_along(ids)), provenance[rep(1L, length(ids)), , drop = FALSE])
        rownames(out) <- NULL
        members[[length(members)+1L]] <- ggchord_add_source_file(out, paths[input], source_file)
      }
    }
  }
  pair_table <- ggchord_rbind_fill(pairs) %||% data.frame(qfeature_id = character(),
    sfeature_id = character(), source_score = character(), .source_input = integer(),
    .source_record = integer(), .source_line = integer())
  membership <- ggchord_rbind_fill(members) %||% data.frame(group_id = character(),
    feature_id = character(), .source_group = integer(), .source_member = integer(),
    .source_input = integer(), .source_record = integer(), .source_line = integer())
  ids <- unique(c(pair_table$qfeature_id, pair_table$sfeature_id, membership$feature_id))
  missing <- if (is.null(feature_data)) ids else setdiff(ids, feature_data$feature_id)
  if (!is.null(feature_data) && length(missing) && unmatched == 'error')
    ggchord_stop(caller, ': unmatched feature_id: ', paste(head(missing, 8L), collapse = ', '),
      '; supply exact IDs or use unmatched="report"')
  anchor <- function(ids) {
    index <- match(ids, feature_data$feature_id)
    data.frame(accver = if (is.null(feature_data)) rep(NA_character_, length(ids)) else feature_data$accver[index],
      start = if (is.null(feature_data)) rep(NA_real_, length(ids)) else feature_data$start[index],
      end = if (is.null(feature_data)) rep(NA_real_, length(ids)) else feature_data$end[index],
      strand = if (is.null(feature_data) || !'strand' %in% names(feature_data)) rep('.', length(ids)) else feature_data$strand[index],
      mapped = !is.na(index))
  }
  q <- anchor(pair_table$qfeature_id); s <- anchor(pair_table$sfeature_id)
  pair_table$q_mapped <- q$mapped; pair_table$s_mapped <- s$mapped
  links <- cbind(data.frame(qaccver = q$accver, saccver = s$accver,
    qstart = q$start, qend = q$end, sstart = s$start, send = s$end,
    qfeature_strand = q$strand, sfeature_strand = s$strand,
    orientation_known = rep(FALSE, nrow(pair_table)), direction = rep(NA_character_, nrow(pair_table))), pair_table)
  placed <- q$mapped & s$mapped
  a <- anchor(membership$feature_id)
  membership <- cbind(membership, a)
  list(pairs = pair_table, membership = membership, ribbon_data = links[placed, , drop = FALSE],
    report = list(format = format, records = records, raw_files = stats::setNames(raw, paths),
      feature_data = feature_data, unmapped_ids = missing, unplaced_pairs = links[!placed, , drop = FALSE],
      coordinates = 'user-declared ascending 1-based-inclusive feature anchors; no source length inferred',
      relationship_semantics = 'original pairs or memberships only; no group expansion or alignment orientation inference'))
}
