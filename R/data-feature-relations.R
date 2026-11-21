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

#' Read existing feature pairs or group memberships
#'
#' Keep IDs exactly as supplied. A pair can be placed using feature
#' coordinates passed explicitly. Group inputs retain memberships and
#' can be used to construct separately declared comparisons.
#'
#' @param file Local relationship-file paths.
#' @param format `"pairs"`, `"named_groups"`, or `"groups"`.
#' @param features Optional feature table to place exact IDs.
#' @param unmatched `"error"` rejects unknown IDs when features are given.
#'   `"report"` records unplaced observations.
#' @param source_file Whether to include row-level source paths.
#' @return A list with pairs, membership, links$feature_relations and report.
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
#' f <- tempfile(fileext=".txt")
#' writeLines("g1 g3", f)
#' relations <- read_feature_relations(f, features=genes)
#' ggchord(seqs) + geom_link_ribbon(data=relations$links$feature_relations)
#' unlink(f)
read_feature_relations <- function(file, format = "pairs", features = NULL,
    unmatched = "error", source_file = FALSE) {
  caller <- "read_feature_relations()"
  format <- match.arg(format, c("pairs", "named_groups", "groups"))
  unmatched <- match.arg(unmatched, c("error", "report"))
  feature_data <- features
  if (!is.logical(source_file) || length(source_file) != 1L || is.na(source_file))
    ggchord_stop(caller, ': source_file must be TRUE or FALSE')
  paths <- ggchord_resolve_files(file, NULL, 'read_feature_relations')
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
  list(pairs = pair_table, membership = membership, links = list(feature_relations = links[placed, , drop = FALSE]),
    report = list(format = format, records = records, raw_files = setNames(raw, paths),
      feature_data = feature_data, unmapped_ids = missing, unplaced_pairs = links[!placed, , drop = FALSE],
      coordinates = 'user-declared ascending 1-based-inclusive feature anchors; no source length inferred',
      relationship_semantics = 'original pairs or memberships only; no group expansion or alignment orientation inference'))
}
