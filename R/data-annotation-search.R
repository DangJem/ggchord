ggchord_search_topology <- function(data, circular, ids, caller) {
  if (!is.null(circular)) {
    if (!is.logical(circular) || length(circular) != 1L || is.na(circular))
      ggchord_stop(caller, ": circular must be TRUE, FALSE, or NULL")
    return(setNames(rep(circular, length(ids)), ids))
  }
  known <- if (is.data.frame(data)) data$topology
  if (!is.null(known) && any(!is.na(known) & !known %in% c("linear", "circular")))
    ggchord_stop(caller, ": topology must be 'linear' or 'circular'")
  setNames((known %||% rep(NA_character_, length(ids))) %in% "circular", ids)
}

ggchord_segment_search_matches <- function(start, end, length) {
  lapply(seq_along(start), function(i) {
    if (start[i] <= end[i]) return(NULL)
    data.frame(start = c(start[i], 1), end = c(length[i], end[i]),
      segment_index = 1:2, segment_type = "standard")
  })
}

ggchord_feature_search_codes <- function(sequence, genetic_code, ids, caller) {
  source <- if (is.data.frame(sequence)) sequence$genetic_code
  if (is.null(source)) source <- rep(NA_real_, length(ids))
  if (!is.null(genetic_code)) {
    ggchord_project_integer(genetic_code, "genetic_code", caller)
    if (is.null(names(genetic_code))) {
      if (length(genetic_code) != 1L) ggchord_stop(caller, ": supply a scalar or name genetic_code by reference")
      supplied <- rep(genetic_code, length(ids))
    } else {
      if (anyDuplicated(names(genetic_code)) || any(!names(genetic_code) %in% ids))
        ggchord_stop(caller, ": genetic_code keys must match unique reference IDs")
      supplied <- unname(genetic_code[ids])
    }
    if (any(!is.na(source) & !is.na(supplied) & source != supplied))
      ggchord_stop(caller, ": genetic_code conflicts with the source")
    source[is.na(source)] <- supplied[is.na(source)]
  }
  if (anyNA(source)) ggchord_stop(caller, ": protein mode requires explicit genetic_code")
  ggchord_project_integer(source, "genetic_code", caller)
}

#' Select restriction sites using explicit criteria
#' @param sites Restriction-site table.
#' @param enzymes Restrict accepted enzyme names.
#' @param min_site_length Minimum recognition sequence length.
#' @param cuts Allowed full-reference cut counts.
#' @param window Reference interval table to select by recognition position.
#' @param commercial_only Keep enzymes with explicit commercial availability.
#' @param representatives Keep all enzymes or the preferred member of equivalent sets.
#' @param profile Explicit enzyme names or the source-backed reference profile.
#' @return A row subset with a selection report.
#' @export
#' @examples
#' found <- find_restriction_sites(c(A="ACGGAATTCGTT"), enzymes="EcoRI", circular=FALSE)
#' shown <- filter_restriction_sites(found, cuts=1, min_site_length=6)
#' shown[c("label", "pos", "cut_top", "cut_bottom")]
filter_restriction_sites <- function(sites, enzymes = NULL, min_site_length = NULL,
    cuts = NULL, window = NULL, commercial_only = FALSE, representatives = "all", profile = NULL) {
  caller <- "filter_restriction_sites()"
  representatives <- match.arg(representatives, c("all", "preferred"))
  ggchord_require_columns(sites, c("accver", "pos", "enzyme", "site_id", "enzyme_site_count"), caller)
  if (!is.null(profile) && (!is.character(profile) || anyNA(profile) || any(!nzchar(profile))))
    ggchord_stop(caller, ": profile must contain nonempty enzyme names or 'reference'")
  if (!is.null(window)) {
    ggchord_require_columns(window, c("accver", "start", "end"), caller)
    ggchord_conversion_coordinates(window, c("start", "end"), caller)
    if (any(window$start > window$end)) ggchord_stop(caller, ": window endpoints must be ordered")
  }
  source <- sites
  sites$position <- sites$pos
  if (!identical(profile, "reference") && !is.null(profile)) {
    enzymes <- if (is.null(enzymes)) profile else intersect(enzymes, profile)
  }
  sets <- if (identical(profile, "reference")) "reference" else "all"
  selected <- ggchord_filter_restriction_profile(sites, set = sets, enzymes = enzymes,
    min_site_length = min_site_length, cuts = cuts, commercial_only = commercial_only)
  if (identical(profile, "reference") && !is.null(enzymes))
    selected <- selected[selected$enzyme %in% enzymes, , drop = FALSE]
  if (!is.null(window)) {
    inside <- vapply(seq_len(nrow(selected)), function(i) {
      relevant <- which(window$accver == selected$accver[i])
      any(selected$pos[i] >= window$start[relevant] & selected$pos[i] <= window$end[relevant])
    }, logical(1))
    selected <- selected[inside, , drop = FALSE]
  }
  if (representatives == "preferred" && nrow(selected)) {
    ggchord_require_columns(selected, c("preferred_enzyme", "is_preferred_enzyme",
      "cut_top", "cut_bottom"), caller)
    key <- selected[c("accver", "preferred_enzyme", "pos", "cut_top", "cut_bottom")]
    index <- order(!selected$is_preferred_enzyme, selected$enzyme, selected$site_id)
    selected <- selected[sort(index[!duplicated(key[index, , drop = FALSE])]), , drop = FALSE]
  }
  keep <- which(source$site_id %in% selected$site_id)
  out <- source[keep, , drop = FALSE]
  ggchord_attach_data_report(out, "filter_restriction_sites", inputs = data_report(source),
    lineage = keep, changes = data.frame(source_row = setdiff(seq_len(nrow(source)), keep),
      action = rep("drop", nrow(source) - length(keep))),
    metadata = list(enzymes = enzymes, min_site_length = min_site_length, cuts = cuts,
      window = window, commercial_only = commercial_only, representatives = representatives, profile = profile))
}
