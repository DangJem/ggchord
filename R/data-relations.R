# Relational preparation never infers alignments or concatenates contigs.

ggchord_relation_ids <- function(x, name, caller) {
  if (!is.atomic(x) || is.list(x)) ggchord_stop(caller, ": ", name, " must be an atomic identifier column")
  value <- as.character(x)
  if (anyNA(value) || any(!nzchar(value))) ggchord_stop(caller, ": ", name,
    " must contain non-missing, non-empty identifiers")
  value
}

#' Validate explicit genome and contig relationships
#'
#' Checks relational keys without deriving genome membership, changing sequence
#' order or joining contigs into a continuous sequence. Original input tables
#' and their metadata are returned unchanged. Sequence lengths and annotation
#' coordinates can be checked separately with [validate_ggchord_data()].
#'
#' A sequence belongs to exactly one genome. Sequence `accver` values must be
#' globally unique, including when local contig names would be reused by two
#' genomes. Supply distinct explicit identifiers before validation; this
#' function never prefixes or rewrites them. Genome rows with no contigs are
#' reported as warnings, and remain in the returned genome table.
#'
#' @param seq_data A data frame containing `accver` and `genome` foreign keys.
#' @param genomes A data frame containing a unique `genome` primary key.
#' @param strict Logical. Stop if relational errors are present; warnings about
#'   genomes without contigs do not cause failure.
#' @return A list with `valid`, unchanged `seq_data` and `genomes`, an `issues`
#'   data frame (`table`, `row`, `column`, `category`, `severity`, `message`),
#'   and `contigs`, a genome-order table with `genome`, `n_contigs`, and
#'   `sequence_rows`/`accvers` list-columns for each valid unique genome key.
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' sequences <- data.frame(accver = c("sample1_c1", "sample1_c2", "sample2_c1"),
#'   length = c(1000, 800, 1200), genome = c("sample1", "sample1", "sample2"))
#' genomes <- data.frame(genome = c("sample1", "sample2"), treatment = c("A", "B"))
#' hierarchy <- validate_ggchord_hierarchy(sequences, genomes)
#' hierarchy$valid
#' hierarchy$contigs
validate_ggchord_hierarchy <- function(seq_data, genomes, strict = FALSE) {
  caller <- "validate_ggchord_hierarchy()"
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) ggchord_stop(
    caller, ": strict must be TRUE or FALSE")
  issues <- list()
  add <- function(table, rows, column, category, severity, message) {
    if (!length(rows)) return(invisible(NULL))
    issues[[length(issues) + 1L]] <<- data.frame(table = table, row = as.integer(rows),
      column = column, category = category, severity = severity, message = message,
      stringsAsFactors = FALSE)
    invisible(NULL)
  }
  check <- function(data, table, required) {
    if (!is.data.frame(data)) {
      add(table, NA_integer_, NA_character_, "not_dataframe", "error", "Input must be a data.frame")
      return(FALSE)
    }
    absent <- setdiff(required, names(data))
    for (name in absent) add(table, NA_integer_, name, "missing_column", "error", "Required key column is absent")
    if (anyDuplicated(names(data))) add(table, NA_integer_, NA_character_,
      "duplicate_column", "error", "Column names must be unique")
    !length(absent) && !anyDuplicated(names(data))
  }
  seq_ok <- check(seq_data, "seq", c("accver", "genome"))
  genome_ok <- check(genomes, "genomes", "genome")
  keys <- function(x, table, column, unique_required) {
    if (!is.atomic(x) || is.list(x)) {
      add(table, NA_integer_, column, "invalid_key", "error", "Keys must be an atomic identifier column")
      return(rep(NA_character_, length(x)))
    }
    x <- as.character(x)
    invalid <- is.na(x) | !nzchar(x)
    add(table, which(invalid), column, "missing_key", "error", "Keys must be non-missing and non-empty")
    if (unique_required) {
      duplicate <- !invalid & (duplicated(x) | duplicated(x, fromLast = TRUE))
      add(table, which(duplicate), column, "duplicate_key", "error", "Primary keys must be unique")
    }
    x
  }
  sid <- membership <- genome <- character()
  if (seq_ok) {
    sid <- keys(seq_data$accver, "seq", "accver", TRUE)
    membership <- keys(seq_data$genome, "seq", "genome", FALSE)
    if (!nrow(seq_data)) add("seq", NA_integer_, "accver", "empty", "error", "At least one contig is required")
  }
  if (genome_ok) {
    genome <- keys(genomes$genome, "genomes", "genome", TRUE)
    if (!nrow(genomes)) add("genomes", NA_integer_, "genome", "empty", "error", "At least one genome is required")
  }
  usable <- !is.na(genome) & nzchar(genome) & !duplicated(genome) & !duplicated(genome, fromLast = TRUE)
  contigs <- data.frame(genome = genome[usable], n_contigs = integer(sum(usable)), stringsAsFactors = FALSE)
  sequence_rows <- accvers <- vector("list", nrow(contigs))
  if (seq_ok && genome_ok) {
    missing <- !is.na(membership) & nzchar(membership) & !membership %in% genome[!is.na(genome)]
    add("seq", which(missing), "genome", "unknown_genome", "error", "Genome foreign key is absent from genomes")
    for (i in seq_len(nrow(contigs))) {
      rows <- which(!is.na(membership) & membership == contigs$genome[i])
      sequence_rows[[i]] <- rows
      accvers[[i]] <- sid[rows]
      contigs$n_contigs[i] <- length(rows)
      if (!length(rows)) add("genomes", which(genome == contigs$genome[i]), "genome",
        "orphan_genome", "warning", "Genome has no contigs in seq_data")
    }
  }
  contigs$sequence_rows <- I(sequence_rows)
  contigs$accvers <- I(accvers)
  issue_table <- if (length(issues)) do.call(rbind, issues) else data.frame(
    table = character(), row = integer(), column = character(), category = character(),
    severity = character(), message = character(), stringsAsFactors = FALSE)
  valid <- !any(issue_table$severity == "error")
  if (strict && !valid) ggchord_stop(caller, ": invalid genome/contig relationships; ",
    "run with strict = FALSE to inspect issues (first error: ",
    issue_table$message[which(issue_table$severity == "error")[1L]], ")")
  list(valid = valid, seq_data = seq_data, genomes = genomes, issues = issue_table, contigs = contigs)
}

#' Convert explicit homology-cluster membership to feature links
#'
#' Connects feature intervals belonging to the same declared cluster for
#' explicitly requested sequence pairs. Membership is an input assertion, not
#' an inferred alignment: no identity, alignment length or score is invented.
#' Both sequence and feature identifiers are kept exactly as supplied.
#'
#' The default `paralogs = "drop"` retains only comparisons containing exactly
#' one member on each side. Ambiguous one-to-many/many-to-many comparisons are
#' omitted and reported. `"all"` explicitly requests their Cartesian product.
#' Its total size is checked against `max_links` before any links are expanded.
#' Pair row order and its query/subject roles are retained. Self-pairs and
#' repeated pairs (including the reversed pair) must be resolved by the caller.
#'
#' Feature IDs normally identify one row. If repeated rows deliberately
#' describe one discontinuous feature, `segment_policy = "span"` requires
#' one sequence, one strand, and non-overlapping intervals. The link then
#' anchors the feature's outer genomic bounds; real segments remain in the
#' `qsegments` and `ssegments` list-columns. This span is a feature anchor,
#' never an exon-to-exon alignment, and does not remove intervening sequence.
#' Nested `segments` fields in a single feature row are retained unchanged.
#'
#' Unknown feature IDs in membership are reported and omitted. Malformed
#' identifiers, conflicting cluster assignments, ambiguous repeated features,
#' unsupported origin-crossing intervals and unknown pair sequence IDs stop
#' with an actionable error. Feature IDs not assigned to a cluster are also
#' listed in the report. Negative strands reverse their link endpoints;
#' unknown strands retain increasing coordinates and `orientation_known = FALSE`.
#'
#' @param feature_data Feature data frame with `feature_id`, `accver`, `start`,
#'   `end`, and optional `strand`, accepted by [as_ggchord_feature()].
#' @param membership Data frame with `feature_id` and `cluster_id`. A feature
#'   can belong to one cluster; extra columns are retained in the report.
#' @param pairs Data frame with explicit `qaccver` and `saccver` comparisons.
#' @param paralogs How to handle multiple members per cluster and sequence:
#'   `"drop"` (default) or explicit `"all"`.
#' @param segment_policy `"error"` (default) rejects repeated feature IDs;
#'   `"span"` explicitly uses a validated discontinuous feature's outer bounds.
#' @param max_links Positive integer upper limit for generated rows, default
#'   10000. Exceeding it stops before Cartesian products are allocated.
#' @return A list with `data`, an ordinary link data frame, and `report`.
#'   Links retain cluster/member IDs, strands, source feature/membership rows
#'   and segment tables. `report$comparisons` lists member counts, output counts
#'   and retention/omission reasons for each requested pair and clusters
#'   represented on either side (clusters absent from both sides are omitted);
#'   `report$members` retains membership metadata with resolution status.
#'   `report$valid` is false when membership contains unknown feature IDs.
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' genes <- data.frame(feature_id = c("gA", "gB"), accver = c("A", "B"),
#'   start = c(100, 400), end = c(200, 500), strand = c("+", "-"))
#' membership <- data.frame(feature_id = c("gA", "gB"), cluster_id = "family1")
#' pairs <- data.frame(qaccver = "A", saccver = "B")
#' homologs <- link_ggchord_clusters(genes, membership, pairs)
#' homologs$data[, c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
link_ggchord_clusters <- function(feature_data, membership, pairs,
    paralogs = c("drop", "all"), segment_policy = c("error", "span"), max_links = 10000) {
  caller <- "link_ggchord_clusters()"
  paralogs <- match.arg(paralogs)
  segment_policy <- match.arg(segment_policy)
  if (!is.numeric(max_links) || length(max_links) != 1L || !is.finite(max_links) ||
      max_links < 1 || max_links != floor(max_links)) ggchord_stop(caller, ": max_links must be a positive integer")
  for (input in c("feature_data", "membership", "pairs")) {
    value <- get(input)
    if (!is.data.frame(value) || anyDuplicated(names(value))) ggchord_stop(
      caller, ": ", input, " must be a data.frame with unique column names")
  }
  ggchord_require_columns(feature_data, "feature_id", caller)
  ggchord_require_columns(membership, c("feature_id", "cluster_id"), caller)
  ggchord_require_columns(pairs, c("qaccver", "saccver"), caller)
  feature_id <- ggchord_relation_ids(feature_data$feature_id, "feature_id", caller)
  member_id <- ggchord_relation_ids(membership$feature_id, "membership$feature_id", caller)
  cluster_id <- ggchord_relation_ids(membership$cluster_id, "cluster_id", caller)
  qid <- ggchord_relation_ids(pairs$qaccver, "pairs$qaccver", caller)
  sid <- ggchord_relation_ids(pairs$saccver, "pairs$saccver", caller)
  if (!nrow(pairs)) ggchord_stop(caller, ": request at least one explicit sequence pair")
  if (any(qid == sid)) ggchord_stop(caller, ": pairs must compare distinct sequences")
  pair_keys <- lapply(seq_along(qid), function(i) unname(sort(c(qid[i], sid[i]))))
  if (anyDuplicated(pair_keys)) ggchord_stop(caller, ": repeated or reversed sequence pairs are not allowed")
  if (anyDuplicated(member_id)) ggchord_stop(caller,
    ": each feature_id must occur once in membership; resolve duplicates or conflicting clusters explicitly")
  features <- as_ggchord_feature(feature_data)
  if (any(features$start > features$end)) ggchord_stop(caller,
    ": feature intervals must be ordered; split origin-crossing features explicitly")
  unknown_pair <- setdiff(unique(c(qid, sid)), features$accver)
  if (length(unknown_pair)) ggchord_stop(caller, ": pair sequence ID absent from feature_data: ", paste(unknown_pair, collapse = ", "))

  ids <- unique(feature_id)
  groups <- split(seq_along(feature_id), factor(feature_id, levels = ids))
  representative <- lapply(seq_along(ids), function(i) {
    rows <- groups[[i]]
    if (length(rows) > 1L && segment_policy == "error") ggchord_stop(caller,
      ": repeated feature_id '", ids[i], "'; use segment_policy = 'span' only for one discontinuous feature")
    data <- features[rows, , drop = FALSE]
    if (length(unique(data$accver)) != 1L || length(unique(data$strand)) != 1L) ggchord_stop(caller,
      ": repeated feature_id '", ids[i], "' must have one sequence and strand")
    ordered <- data[order(data$start, data$end), , drop = FALSE]
    if (nrow(ordered) > 1L && any(ordered$start[-1L] <= utils::head(ordered$end, -1L))) ggchord_stop(caller,
      ": repeated feature_id '", ids[i], "' has overlapping or duplicate segments")
    list(feature_id = ids[i], accver = data$accver[1L], strand = data$strand[1L],
      start = min(data$start), end = max(data$end), source_rows = rows,
      segments = feature_data[rows, , drop = FALSE])
  })
  names(representative) <- ids
  members <- membership
  # Report columns deliberately use reserved names; reject collisions instead
  # of overwriting user metadata during a relational join.
  reserved <- c(".membership_row", ".feature_rows", ".resolved_accver", ".status")
  if (any(reserved %in% names(members))) ggchord_stop(caller, ": membership contains reserved report columns: ",
    paste(intersect(reserved, names(members)), collapse = ", "))
  member_index <- match(member_id, ids)
  resolved <- !is.na(member_index)
  members$.membership_row <- seq_len(nrow(members))
  members$.status <- ifelse(resolved, "resolved", "unknown_feature")
  members$.resolved_accver <- vapply(member_index, function(index) {
    if (!is.na(index)) representative[[index]]$accver else NA_character_
  }, character(1))
  members$.feature_rows <- I(lapply(member_index, function(index) {
    if (!is.na(index)) representative[[index]]$source_rows else integer()
  }))
  clusters <- unique(cluster_id)
  comparisons <- list()
  total <- 0
  for (pair_row in seq_len(nrow(pairs))) {
    relevant <- clusters[clusters %in% cluster_id[resolved &
      members$.resolved_accver %in% c(qid[pair_row], sid[pair_row])]]
    for (cluster in relevant) {
      query <- which(resolved & cluster_id == cluster & members$.resolved_accver == qid[pair_row])
      subject <- which(resolved & cluster_id == cluster & members$.resolved_accver == sid[pair_row])
      potential <- as.double(length(query)) * as.double(length(subject))
      reason <- if (!length(query) || !length(subject)) "missing_member" else
        if (paralogs == "drop" && (length(query) != 1L || length(subject) != 1L)) "multiple_members" else "retained"
      n_output <- if (reason == "retained") potential else 0
      total <- total + n_output
      if (!is.finite(total) || total > max_links) ggchord_stop(caller,
        ": requested membership expansion exceeds max_links = ", max_links,
        "; restrict pairs/members or raise the limit explicitly")
      comparisons[[length(comparisons) + 1L]] <- list(pair_row = pair_row, cluster_id = cluster,
        query = query, subject = subject, n_query = length(query), n_subject = length(subject),
        n_possible = potential, n_output = n_output, reason = reason)
    }
  }
  output <- list()
  for (comparison in comparisons) {
    if (!comparison$n_output) next
    for (qrow in comparison$query) for (srow in comparison$subject) {
      query <- representative[[member_id[qrow]]]
      subject <- representative[[member_id[srow]]]
      known <- query$strand %in% c("+", "-") && subject$strand %in% c("+", "-")
      row <- data.frame(qaccver = query$accver, saccver = subject$accver,
        qstart = if (query$strand == "-") query$end else query$start,
        qend = if (query$strand == "-") query$start else query$end,
        sstart = if (subject$strand == "-") subject$end else subject$start,
        send = if (subject$strand == "-") subject$start else subject$end,
        qstrand = query$strand, sstrand = subject$strand,
        orientation_known = known,
        direction = if (known) if (query$strand == subject$strand) "same" else "reverse" else NA_character_,
        cluster_id = comparison$cluster_id, qfeature_id = query$feature_id, sfeature_id = subject$feature_id,
        .pair_row = comparison$pair_row, .query_membership_row = qrow, .subject_membership_row = srow,
        stringsAsFactors = FALSE)
      row$.query_feature_rows <- I(list(query$source_rows))
      row$.subject_feature_rows <- I(list(subject$source_rows))
      row$qsegments <- I(list(query$segments))
      row$ssegments <- I(list(subject$segments))
      output[[length(output) + 1L]] <- row
    }
  }
  empty <- data.frame(qaccver = character(), saccver = character(), qstart = numeric(), qend = numeric(),
    sstart = numeric(), send = numeric(), qstrand = character(), sstrand = character(),
    orientation_known = logical(), direction = character(), cluster_id = character(),
    qfeature_id = character(), sfeature_id = character(), .pair_row = integer(),
    .query_membership_row = integer(), .subject_membership_row = integer(), stringsAsFactors = FALSE)
  for (name in c(".query_feature_rows", ".subject_feature_rows", "qsegments", "ssegments")) empty[[name]] <- I(list())
  data <- if (length(output)) do.call(rbind, output) else empty
  rownames(data) <- NULL
  comparison_report <- if (length(comparisons)) do.call(rbind, lapply(comparisons, function(x) {
    data.frame(pair_row = x$pair_row, qaccver = qid[x$pair_row], saccver = sid[x$pair_row],
      cluster_id = x$cluster_id, n_query = x$n_query, n_subject = x$n_subject,
      n_possible = x$n_possible, n_output = x$n_output, reason = x$reason, stringsAsFactors = FALSE)
  })) else data.frame(pair_row = integer(), qaccver = character(), saccver = character(), cluster_id = character(),
    n_query = integer(), n_subject = integer(), n_possible = numeric(), n_output = numeric(), reason = character())
  list(data = data, report = list(valid = all(resolved), members = members,
    comparisons = comparison_report, unassigned_feature_ids = setdiff(ids, member_id),
    n_output = nrow(data), paralogs = paralogs, segment_policy = segment_policy,
    source_pairs = pairs, source_features = feature_data,
    source_conversion = attr(features, "ggchord_conversion")))
}
