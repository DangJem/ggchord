# Common-feature database access and deterministic sequence matching.

ggchord_normalize_sequence_input <- function(sequence, caller) {
  if (is.data.frame(sequence)) {
    ggchord_require_columns(sequence, c("accver", "sequence"), caller)
    ids <- as.character(sequence$accver)
    seqs <- as.character(sequence$sequence)
  } else if (is.character(sequence) && length(sequence)) {
    seqs <- as.character(sequence)
    ids <- names(sequence)
    if (is.null(ids)) {
      ids <- if (length(seqs) == 1L) "sequence" else
        paste0("sequence_", seq_along(seqs))
    }
    missing_ids <- is.na(ids) | !nzchar(ids)
    ids[missing_ids] <- paste0("sequence_", which(missing_ids))
  } else {
    ggchord_stop(caller, ": invalid sequence input")
  }
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    ggchord_stop(caller, ": sequence IDs must be non-empty and unique")
  }
  seqs <- toupper(gsub("[[:space:]]", "", seqs))
  if (anyNA(seqs) || any(!nzchar(seqs)) ||
      any(!grepl("^[ACGTRYSWKMBDHVN]+$", seqs))) {
    ggchord_stop(caller, ": invalid DNA symbols")
  }
  list(ids = ids, sequences = seqs)
}

ggchord_builtin_common_features <- function() {
  if (!exists("ggchord_common_feature_database", inherits = TRUE)) {
    ggchord_stop(
      "The internal common-feature database is missing; reinstall ggchord"
    )
  }
  get("ggchord_common_feature_database", inherits = TRUE)
}

ggchord_common_database <- function(database) {
  if (is.null(database)) return(ggchord_builtin_common_features())
  if (is.data.frame(database)) {
    x <- as.data.frame(database, stringsAsFactors = FALSE)
    name_col <- intersect(c("name", "anno", "feature_name"), names(x))[1L]
    dna_col <- intersect(c(
      "reference_dna", "sequence", "dna_sequence",
      "reference_dna_feature_5to3"
    ), names(x))[1L]
    if (is.na(name_col) || is.na(dna_col)) {
      ggchord_stop(
        "find_common_features(): a custom data frame needs a name/anno and ",
        "sequence/reference_dna column"
      )
    }
    if (!"common_feature_id" %in% names(x)) {
      id_col <- intersect(c("feature_id", "id"), names(x))[1L]
      x$common_feature_id <- if (is.na(id_col)) {
        sprintf("custom_feature_%04d", seq_len(nrow(x)))
      } else as.character(x[[id_col]])
    }
    x$name <- as.character(x[[name_col]])
    if (!"type" %in% names(x)) x$type <- "misc_feature"
    if (!"directionality_label" %in% names(x)) {
      strand_col <- intersect(c("strand", "directionality"), names(x))[1L]
      x$directionality_label <- if (is.na(strand_col)) "forward" else
        ifelse(x[[strand_col]] == "-", "reverse",
          ifelse(x[[strand_col]] %in% c(".", "0", "nondirectional"),
            "nondirectional", "forward"))
    }
    x$reference_dna_feature_5to3 <- as.character(x[[dna_col]])
    if (!"reference_protein" %in% names(x)) x$reference_protein <- NA_character_
    if (!"detectionMode" %in% names(x)) x$detectionMode <- NA_character_
    if (!"geneticCode" %in% names(x)) x$geneticCode <- NA_integer_
    segments <- data.frame(
      feature_id = x$common_feature_id,
      segment_index = 1L, segment_type = "standard", start = 1L,
      end = nchar(x$reference_dna_feature_5to3),
      length_bp = nchar(x$reference_dna_feature_5to3), translated = FALSE,
      segment_name = NA_character_,
      dna_sequence_top_strand = x$reference_dna_feature_5to3,
      stringsAsFactors = FALSE
    )
    return(list(
      metadata = list(
        database_id = "custom", database_version = "custom",
        source_workbook_sha256 = NA_character_
      ),
      features = x, segments = segments,
      qualifiers = data.frame(), qualifier_links = data.frame(),
      feature_type_summary = data.frame(), search_indexes = list()
    ))
  }
  if (!is.list(database) ||
      !all(c("metadata", "features", "segments") %in% names(database))) {
    ggchord_stop(
      "find_common_features(): database must be a normalized database list ",
      "or a custom data frame"
    )
  }
  database
}

ggchord_common_feature_columns <- function(features) {
  if (!"common_feature_id" %in% names(features)) {
    features$common_feature_id <- as.character(features$feature_id)
  }
  if (!"name" %in% names(features) && "anno" %in% names(features)) {
    features$name <- as.character(features$anno)
  }
  features
}

ggchord_common_pattern_parts <- function(segments) {
  segments <- segments[order(segments$segment_index), , drop = FALSE]
  dna <- toupper(gsub(
    "[[:space:]]", "", as.character(segments$dna_sequence_top_strand)
  ))
  dna[is.na(dna)] <- ""
  length_bp <- as.integer(segments$length_bp)
  length_bp[is.na(length_bp)] <- nchar(dna[is.na(length_bp)])
  data.frame(
    segment_index = as.integer(segments$segment_index),
    segment_type = as.character(segments$segment_type),
    length_bp = length_bp, dna = dna,
    segment_name = if ("segment_name" %in% names(segments))
      as.character(segments$segment_name) else NA_character_,
    color = if ("color" %in% names(segments))
      as.character(segments$color) else NA_character_,
    translated = if ("translated" %in% names(segments))
      as.logical(segments$translated) else FALSE,
    stringsAsFactors = FALSE
  )
}

ggchord_common_regex <- function(parts, reverse = FALSE) {
  if (reverse) {
    parts <- parts[nrow(parts):1L, , drop = FALSE]
    fixed <- parts$segment_type != "gap"
    parts$dna[fixed] <- vapply(
      parts$dna[fixed], ggchord_reverse_complement, character(1)
    )
  }
  paste0(vapply(seq_len(nrow(parts)), function(i) {
    if (parts$segment_type[i] == "gap") {
      paste0("[ACGT]{", parts$length_bp[i], "}")
    } else ggchord_iupac_regex(parts$dna[i])
  }, character(1)), collapse = "")
}

ggchord_overlapping_regex_starts <- function(pattern, target, max_start) {
  hit <- gregexpr(paste0("(?=", pattern, ")"), target, perl = TRUE)[[1L]]
  if (length(hit) == 1L && hit[1L] < 0L) return(integer())
  hit[hit <= max_start]
}

ggchord_wrap_position <- function(x, length, circular) {
  if (circular) ((x - 1) %% length) + 1 else x
}

ggchord_common_target_segments <- function(parts, hit, seq_length,
                                            reverse, circular) {
  offsets <- cumsum(c(0L, utils::head(parts$length_bp, -1L)))
  if (!reverse) {
    starts <- hit + offsets
    ends <- starts + parts$length_bp - 1L
  } else {
    span <- sum(parts$length_bp)
    starts <- hit + span - offsets - parts$length_bp
    ends <- starts + parts$length_bp - 1L
  }
  data.frame(
    segment_index = parts$segment_index,
    segment_type = parts$segment_type,
    start = ggchord_wrap_position(starts, seq_length, circular),
    end = ggchord_wrap_position(ends, seq_length, circular),
    length_bp = parts$length_bp,
    translated = parts$translated,
    segment_name = parts$segment_name,
    color = parts$color,
    stringsAsFactors = FALSE
  )
}

ggchord_common_strand <- function(directionality, observed) {
  directionality <- tolower(directionality %||% "forward")
  if (directionality == "nondirectional") return(".")
  if (directionality == "bidirectional") return("+/-")
  if (directionality == "reverse") {
    return(if (observed == "+") "-" else "+")
  }
  observed
}

ggchord_common_candidate <- function(accver, feature, method, start, end,
                                     strand, cross_origin, segments,
                                     identity = 1, coverage = 1,
                                     metadata = list()) {
  id <- as.character(feature$common_feature_id)
  db_id <- as.character(metadata$database_id %||% "ggchord-common-features")
  db_version <- as.character(metadata$database_version %||%
    metadata$export_version %||% NA_character_)
  source_sha <- as.character(metadata$source_workbook_sha256 %||% NA_character_)
  feature_fill <- if ("segment_colors" %in% names(feature)) {
    strsplit(as.character(feature$segment_colors)[1L], ",", fixed = TRUE)[[1L]][1L]
  } else NA_character_
  if (is.na(feature_fill) || !nzchar(feature_fill)) feature_fill <- "#B8BDC3"
  feature_label_colour <- ggchord_contrast_colour(feature_fill)
  feature_class <- ggchord_feature_classes(feature$type)
  feature_shape <- ggchord_feature_class_shape(feature_class)
  if (strand == "." && feature_shape %in% c(
      "arrow", "compact_arrow", "promoter_arrow", "primer_arrow")) {
    feature_shape <- "block"
  }
  preferred_lane <- ggchord_feature_preferred_lane(feature_class)
  data.frame(
    accver = accver,
    start = as.integer(start), end = as.integer(end), strand = strand,
    cross_origin = isTRUE(cross_origin),
    type = as.character(feature$type), anno = as.character(feature$name),
    feature_color = feature_fill, feature_label_colour = feature_label_colour,
    feature_class = feature_class, feature_shape = feature_shape,
    preferred_lane = preferred_lane, .semantic_lane_hint = TRUE,
    common_feature_id = id,
    match_id = paste(accver, id, strand, start, end, method, sep = ":"),
    match_method = method, identity = as.numeric(identity),
    coverage = as.numeric(coverage),
    confidence = if (method == "protein_approx") {
      if (identity >= .97 && coverage >= .95) "high" else "medium"
    } else "high",
    segment_count = as.integer(nrow(segments)),
    segments = I(list(segments)), database_id = db_id,
    database_version = db_version, source_sha256 = source_sha,
    stringsAsFactors = FALSE
  )
}

ggchord_common_reference_candidates <- function(target, accver, database,
                                                 types = NULL,
                                                 features = NULL) {
  if (!all(c("reference_sequences", "reference_features") %in%
      names(database))) return(list())
  sequence_rows <- database$reference_sequences[
    database$reference_sequences$sequence == target, , drop = FALSE]
  if (!nrow(sequence_rows)) return(list())
  reference_ids <- sequence_rows$reference_id
  annotations <- database$reference_features[
    database$reference_features$reference_id %in% reference_ids, , drop = FALSE]
  if (!is.null(types)) {
    annotations <- annotations[tolower(annotations$type) %in%
      tolower(types), , drop = FALSE]
  }
  if (!is.null(features)) {
    annotations <- annotations[
      tolower(annotations$name) %in% tolower(features) |
        annotations$common_feature_id %in% features, , drop = FALSE]
  }
  if (!nrow(annotations)) return(list())
  lapply(seq_len(nrow(annotations)), function(i) {
    annotation <- annotations[i, , drop = FALSE]
    feature <- data.frame(
      common_feature_id = annotation$common_feature_id,
      name = annotation$name, type = annotation$type,
      segment_colors = annotation$color, stringsAsFactors = FALSE
    )
    candidate <- ggchord_common_candidate(
      accver, feature, "reference_exact", annotation$start, annotation$end,
      annotation$strand, annotation$start > annotation$end,
      annotation$segments[[1L]], metadata = database$metadata
    )
    source_row <- sequence_rows[
      match(annotation$reference_id, sequence_rows$reference_id), , drop = FALSE]
    candidate$source_sha256 <- source_row$source_sha256
    candidate
  })
}

ggchord_match_common_dna <- function(target, accver, feature, segments,
                                     circular, metadata) {
  parts <- ggchord_common_pattern_parts(segments)
  if (!nrow(parts) || any(parts$length_bp < 0L) ||
      any(parts$segment_type != "gap" & !nzchar(parts$dna))) return(list())
  span <- sum(parts$length_bp)
  target_length <- nchar(target)
  if (!span || span > target_length) return(list())
  search <- if (circular && span > 1L) {
    paste0(target, substr(target, 1L, span - 1L))
  } else target
  plus_pattern <- ggchord_common_regex(parts, FALSE)
  minus_pattern <- ggchord_common_regex(parts, TRUE)
  patterns <- list("+" = plus_pattern)
  if (!identical(plus_pattern, minus_pattern)) patterns[["-"]] <- minus_pattern
  out <- list()
  for (observed in names(patterns)) {
    hits <- ggchord_overlapping_regex_starts(
      patterns[[observed]], search, target_length
    )
    for (hit in hits) {
      target_segments <- ggchord_common_target_segments(
        parts, hit, target_length, observed == "-", circular
      )
      finish <- hit + span - 1L
      out[[length(out) + 1L]] <- ggchord_common_candidate(
        accver, feature, "dna_exact", hit,
        ggchord_wrap_position(finish, target_length, circular),
        ggchord_common_strand(feature$directionality_label, observed),
        finish > target_length, target_segments, metadata = metadata
      )
    }
  }
  out
}

ggchord_match_common_dna_near_exact <- function(target, accver, feature,
                                                 metadata) {
  reference <- toupper(gsub("[[:space:]]", "",
    as.character(feature$reference_dna_feature_5to3)))
  if (length(reference) != 1L || is.na(reference) ||
      nchar(reference) < 100L || nchar(reference) > min(1500L, nchar(target)) ||
      !grepl("^[ACGTRYSWKMBDHVN]+$", reference)) return(list())
  patterns <- list("+" = reference,
    "-" = ggchord_reverse_complement(reference))
  out <- list()
  for (observed in names(patterns)) {
    pattern <- patterns[[observed]]
    seed_starts <- unique(round(seq(1L, nchar(pattern) - 19L,
      length.out = min(12L, max(1L, floor(nchar(pattern) / 20L))))))
    if (!any(vapply(seed_starts, function(i) {
      grepl(substr(pattern, i, i + 19L), target, fixed = TRUE)
    }, logical(1)))) next
    distance <- utils::adist(pattern, target,
      partial = TRUE, counts = TRUE, fixed = TRUE)
    counts <- attr(distance, "counts")
    offsets <- attr(distance, "offsets")
    substitutions <- counts[1L, 1L, "sub"]
    if (counts[1L, 1L, "ins"] != 0L || counts[1L, 1L, "del"] != 0L ||
        substitutions > max(1L, floor(nchar(reference) * .01))) next
    start <- offsets[1L, 1L, "first"]
    end <- offsets[1L, 1L, "last"]
    segments <- data.frame(
      segment_index = 1L, segment_type = "standard", start = start, end = end,
      length_bp = end - start + 1L, translated = FALSE,
      segment_name = NA_character_,
      color = if ("segment_colors" %in% names(feature)) {
        strsplit(as.character(feature$segment_colors)[1L], ",", fixed = TRUE)[[1L]][1L]
      } else NA_character_,
      stringsAsFactors = FALSE
    )
    out[[length(out) + 1L]] <- ggchord_common_candidate(
      accver, feature, "dna_near_exact", start, end,
      ggchord_common_strand(feature$directionality_label, observed), FALSE,
      segments, 1 - substitutions / nchar(reference), 1, metadata
    )
  }
  out
}

ggchord_genetic_code <- local({
  bases <- strsplit(c(
    "TTT F TTC F TTA L TTG L TCT S TCC S TCA S TCG S TAT Y TAC Y TAA * TAG *",
    "TGT C TGC C TGA * TGG W CTT L CTC L CTA L CTG L CCT P CCC P CCA P CCG P",
    "CAT H CAC H CAA Q CAG Q CGT R CGC R CGA R CGG R ATT I ATC I ATA I ATG M",
    "ACT T ACC T ACA T ACG T AAT N AAC N AAA K AAG K AGT S AGC S AGA R AGG R",
    "GTT V GTC V GTA V GTG V GCT A GCC A GCA A GCG A GAT D GAC D GAA E GAG E",
    "GGT G GGC G GGA G GGG G"
  ), "[[:space:]]+")
  tokens <- unlist(bases)
  stats::setNames(tokens[seq(2L, length(tokens), 2L)],
                  tokens[seq(1L, length(tokens), 2L)])
})

ggchord_translate_frame <- function(sequence, frame) {
  last <- nchar(sequence) - ((nchar(sequence) - frame + 1L) %% 3L)
  if (last < frame + 2L) return("")
  starts <- seq.int(frame, last - 2L, by = 3L)
  codons <- substring(sequence, starts, starts + 2L)
  aa <- unname(ggchord_genetic_code[codons])
  aa[is.na(aa)] <- "X"
  paste0(aa, collapse = "")
}

ggchord_normalize_protein <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) return("")
  toupper(gsub("[^A-Z*]", "", x))
}

ggchord_common_protein_frames <- function(target) {
  rc <- ggchord_reverse_complement(target)
  out <- list()
  for (strand in c("+", "-")) for (frame in 1:3) {
    dna <- if (strand == "+") target else rc
    out[[paste0(strand, frame)]] <- list(
      strand = strand, frame = frame,
      protein = ggchord_translate_frame(dna, frame)
    )
  }
  out
}

ggchord_common_protein_coordinates <- function(frame, aa_start, aa_length,
                                                seq_length) {
  rc_start <- frame$frame + (aa_start - 1L) * 3L
  rc_end <- rc_start + aa_length * 3L - 1L
  if (frame$strand == "+") c(start = rc_start, end = rc_end) else
    c(start = seq_length - rc_end + 1L, end = seq_length - rc_start + 1L)
}

ggchord_common_approx_hits <- function(reference, target,
                                       min_identity, min_coverage) {
  nr <- nchar(reference); nt <- nchar(target)
  if (nr < 40L || nt < 12L) return(data.frame())
  k <- min(8L, max(5L, floor(nr / 12L)))
  ref_positions <- unique(c(seq.int(1L, max(1L, nr - k + 1L),
    by = max(1L, floor((nr - k + 1L) / 12L))), nr - k + 1L))
  starts <- integer()
  for (rp in ref_positions) {
    seed <- substr(reference, rp, rp + k - 1L)
    hit <- gregexpr(seed, target, fixed = TRUE)[[1L]]
    if (!(length(hit) == 1L && hit[1L] < 0L)) starts <- c(starts, hit - rp + 1L)
  }
  starts <- sort(unique(starts))
  rows <- lapply(starts, function(s) {
    ref_lo <- max(1L, 2L - s)
    ref_hi <- min(nr, nt - s + 1L)
    overlap <- ref_hi - ref_lo + 1L
    if (overlap <= 0L) return(NULL)
    tar_lo <- s + ref_lo - 1L
    a <- strsplit(substr(reference, ref_lo, ref_hi), "", fixed = TRUE)[[1L]]
    b <- strsplit(substr(target, tar_lo, tar_lo + overlap - 1L), "", fixed = TRUE)[[1L]]
    identity <- mean(a == b)
    coverage <- overlap / nr
    if (identity < min_identity || coverage < min_coverage) return(NULL)
    data.frame(start = tar_lo, identity = identity, coverage = coverage,
      aa_length = overlap, ref_lo = ref_lo)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows) && length(starts) && nr <= 200L) {
    distance <- utils::adist(reference, target, partial = TRUE, counts = TRUE)
    counts <- attr(distance, "counts")
    offsets <- attr(distance, "offsets")
    aligned_reference <- nr - counts[1L, 1L, "del"]
    identity <- if (aligned_reference > 0L) {
      (aligned_reference - counts[1L, 1L, "sub"]) / aligned_reference
    } else 0
    coverage <- aligned_reference / nr
    if (identity >= min_identity && coverage >= min_coverage) {
      rows[[1L]] <- data.frame(
        start = offsets[1L, 1L, "first"], identity = identity,
        coverage = coverage,
        aa_length = offsets[1L, 1L, "last"] - offsets[1L, 1L, "first"] + 1L,
        ref_lo = 1L
      )
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(-out$coverage, -out$identity, out$start), , drop = FALSE]
}

ggchord_match_common_protein <- function(target, accver, feature, frames,
                                         approximate, min_identity,
                                         min_coverage, metadata) {
  reference_with_stop <- ggchord_normalize_protein(feature$reference_protein)
  if (!nzchar(reference_with_stop)) return(list())
  reference <- reference_with_stop
  reference <- sub("\\*$", "", reference)
  out <- list(); target_length <- nchar(target)
  for (frame in frames) {
    if (!approximate) {
      hit <- gregexpr(reference, frame$protein, fixed = TRUE)[[1L]]
      if (length(hit) == 1L && hit[1L] < 0L) next
      matches <- data.frame(
        start = hit, identity = 1, coverage = 1,
        aa_length = nchar(reference), ref_lo = 1L
      )
      method <- "protein_exact"
    } else {
      matches <- ggchord_common_approx_hits(
        reference, frame$protein, min_identity, min_coverage
      )
      method <- "protein_approx"
    }
    if (!nrow(matches)) next
    for (i in seq_len(nrow(matches))) {
      include_stop <- !approximate && endsWith(reference_with_stop, "*") &&
        isTRUE(as.logical(feature$hitsStopCodon %||% FALSE)) &&
        substr(frame$protein,
          matches$start[i] + matches$aa_length[i],
          matches$start[i] + matches$aa_length[i]) == "*"
      matched_aa_length <- matches$aa_length[i] + as.integer(include_stop)
      coords <- ggchord_common_protein_coordinates(
        frame, matches$start[i], matched_aa_length, target_length
      )
      seg <- data.frame(
        segment_index = 1L, segment_type = "standard",
        start = as.integer(coords["start"]), end = as.integer(coords["end"]),
        length_bp = as.integer(matched_aa_length * 3L), translated = TRUE,
        segment_name = NA_character_,
        color = if ("segment_colors" %in% names(feature)) {
          strsplit(as.character(feature$segment_colors)[1L], ",", fixed = TRUE)[[1L]][1L]
        } else NA_character_,
        stringsAsFactors = FALSE
      )
      out[[length(out) + 1L]] <- ggchord_common_candidate(
        accver, feature, method, coords["start"], coords["end"],
        ggchord_common_strand(feature$directionality_label, frame$strand),
        FALSE, seg, matches$identity[i], matches$coverage[i], metadata
      )
    }
  }
  out
}

ggchord_common_overlap <- function(a, b, length) {
  intervals <- function(start, end) {
    if (start <= end) matrix(c(start, end), nrow = 1L) else
      rbind(c(start, length), c(1, end))
  }
  aa <- intervals(a$start, a$end); bb <- intervals(b$start, b$end)
  overlap <- 0
  for (i in seq_len(nrow(aa))) for (j in seq_len(nrow(bb))) {
    overlap <- overlap + max(0, min(aa[i, 2], bb[j, 2]) -
      max(aa[i, 1], bb[j, 1]) + 1)
  }
  min(overlap / sum(aa[, 2] - aa[, 1] + 1),
      overlap / sum(bb[, 2] - bb[, 1] + 1))
}

ggchord_resolve_common_candidates <- function(candidates, lengths) {
  if (!nrow(candidates)) return(candidates)
  method_rank <- match(candidates$match_method,
    c("reference_exact", "dna_exact", "protein_exact", "dna_near_exact",
      "protein_approx"))
  span <- vapply(seq_len(nrow(candidates)), function(i) {
    if (candidates$start[i] <= candidates$end[i])
      candidates$end[i] - candidates$start[i] + 1 else
      lengths[[candidates$accver[i]]] - candidates$start[i] +
        candidates$end[i] + 1
  }, numeric(1))
  ord <- order(candidates$accver, candidates$type, -span,
    -candidates$coverage, -candidates$identity, method_rank,
    candidates$common_feature_id, candidates$start)
  keep <- logical(nrow(candidates))
  for (i in ord) {
    rivals <- which(keep & candidates$accver == candidates$accver[i] &
      candidates$type == candidates$type[i])
    duplicate <- any(vapply(rivals, function(j) {
      ggchord_common_overlap(
        candidates[i, ], candidates[j, ], lengths[[candidates$accver[i]]]
      ) >= .40
    }, logical(1)))
    if (!duplicate) keep[i] <- TRUE
  }
  candidates[keep, , drop = FALSE]
}

ggchord_empty_common_features <- function() {
  data.frame(
    accver = character(), start = integer(), end = integer(),
    strand = character(), cross_origin = logical(), type = character(),
    anno = character(), feature_color = character(),
    feature_label_colour = character(), feature_class = character(),
    feature_shape = character(), preferred_lane = integer(),
    .semantic_lane_hint = logical(),
    common_feature_id = character(),
    match_id = character(), match_method = character(), identity = numeric(),
    coverage = numeric(), confidence = character(), segment_count = integer(),
    segments = I(list()), database_id = character(),
    database_version = character(), source_sha256 = character(),
    stringsAsFactors = FALSE
  )
}

#' Find common biological features in DNA sequences
#'
#' Searches the built-in common-feature database (or a compatible custom
#' database) using curated exact-sequence annotations, exact DNA, and the
#' detection mode stored for each database record. Restriction sites are
#' intentionally outside this API; use
#' [find_restriction_sites()] for those.
#'
#' @param sequence DNA text, a named character vector, or a data frame with
#'   `accver` and `sequence`.
#' @param database Optional normalized database list or custom data frame.
#' @param types,features Optional feature-type and feature-name/ID filters.
#' @param circular Search across the sequence origin.
#' @param mode Matching mode. `"auto"` first uses curated annotations for an
#'   exactly known sequence, then follows each record's detection mode. Records
#'   marked `exactProteinMatch` use exact protein followed by near-exact DNA;
#'   other records use exact DNA. `"protein"` explicitly enables exact and
#'   approximate protein matching.
#' @param min_protein_identity,min_protein_coverage Approximate protein
#'   thresholds in `[0, 1]`.
#' @param resolve Return deterministic best annotations or every candidate.
#' @return A data frame directly usable by [geom_feature()] and
#'   [geom_feature_label_repel()]. `feature_class`, `feature_shape`, and
#'   `preferred_lane` provide normalized semantic layout hints; explicit user
#'   aesthetics and scales still take priority. Its `segments` list-column
#'   preserves the biological feature's segment structure.
#' @export
find_common_features <- function(
    sequence, database = NULL, types = NULL, features = NULL,
    circular = TRUE, mode = c("auto", "dna", "protein"),
    min_protein_identity = 0.90, min_protein_coverage = 0.80,
    resolve = c("best", "all")) {
  mode <- match.arg(mode); resolve <- match.arg(resolve)
  if (!is.logical(circular) || length(circular) != 1L || is.na(circular)) {
    ggchord_stop("find_common_features(): circular must be TRUE or FALSE")
  }
  thresholds <- c(min_protein_identity, min_protein_coverage)
  if (!is.numeric(thresholds) || any(!is.finite(thresholds)) ||
      any(thresholds < 0 | thresholds > 1)) {
    ggchord_stop("find_common_features(): protein thresholds must lie in [0, 1]")
  }
  input <- ggchord_normalize_sequence_input(sequence, "find_common_features()")
  db <- ggchord_common_database(database)
  feature_table <- ggchord_common_feature_columns(
    as.data.frame(db$features, stringsAsFactors = FALSE)
  )
  segment_table <- as.data.frame(db$segments, stringsAsFactors = FALSE)
  if (!all(c("common_feature_id", "name", "type") %in% names(feature_table)) ||
      !all(c("feature_id", "segment_index", "segment_type", "length_bp",
             "dna_sequence_top_strand") %in% names(segment_table))) {
    ggchord_stop("find_common_features(): malformed feature database")
  }
  if (!is.null(types)) {
    if (!is.character(types) || anyNA(types))
      ggchord_stop("find_common_features(): types must be character")
    feature_table <- feature_table[tolower(feature_table$type) %in%
      tolower(types), , drop = FALSE]
  }
  if (!is.null(features)) {
    if (!is.character(features) || anyNA(features))
      ggchord_stop("find_common_features(): features must be character")
    feature_table <- feature_table[
      tolower(feature_table$name) %in% tolower(features) |
        feature_table$common_feature_id %in% features, , drop = FALSE]
  }
  if (!nrow(feature_table)) return(ggchord_empty_common_features())
  split_segments <- split(segment_table,
    factor(segment_table$feature_id, levels = unique(segment_table$feature_id)))
  candidates <- list()
  for (s in seq_along(input$sequences)) {
    target <- input$sequences[s]; accver <- input$ids[s]
    if (mode == "auto") {
      reference_hits <- ggchord_common_reference_candidates(
        target, accver, db, types = types, features = features
      )
      if (length(reference_hits)) {
        candidates <- c(candidates, reference_hits)
        next
      }
    }
    frames <- NULL
    dna_matched <- character(); protein_exact_matched <- character()
    if (mode %in% c("auto", "dna")) {
      for (i in seq_len(nrow(feature_table))) {
        feature <- feature_table[i, , drop = FALSE]
        protein_mode <- identical(
          as.character(feature$detectionMode %||% NA_character_),
          "exactProteinMatch"
        )
        if (mode == "auto" && protein_mode) next
        seg <- split_segments[[as.character(feature$common_feature_id)]]
        if (is.null(seg)) next
        hit <- ggchord_match_common_dna(
          target, accver, feature, seg, circular, db$metadata
        )
        if (length(hit)) {
          candidates <- c(candidates, hit)
          dna_matched <- c(dna_matched, as.character(feature$common_feature_id))
        }
      }
    }
    if (mode %in% c("auto", "protein")) {
      frames <- ggchord_common_protein_frames(target)
      for (i in seq_len(nrow(feature_table))) {
        feature <- feature_table[i, , drop = FALSE]
        id <- as.character(feature$common_feature_id)
        protein_mode <- identical(
          as.character(feature$detectionMode %||% NA_character_),
          "exactProteinMatch"
        )
        if (mode == "auto" && !protein_mode) next
        if (mode == "auto" && id %in% dna_matched) next
        hit <- ggchord_match_common_protein(
          target, accver, feature, frames, FALSE,
          min_protein_identity, min_protein_coverage, db$metadata
        )
        if (length(hit)) {
          candidates <- c(candidates, hit)
          protein_exact_matched <- c(protein_exact_matched, id)
        }
      }
      protein_modes <- !is.na(feature_table$detectionMode) &
        feature_table$detectionMode == "exactProteinMatch"
      if (mode == "auto") {
        approximate_rows <- which(protein_modes)
      } else {
        translated <- if ("translated_any" %in% names(feature_table)) {
          !is.na(feature_table$translated_any) & feature_table$translated_any
        } else feature_table$type == "CDS"
        approximate_rows <- which(feature_table$type == "CDS" & translated)
      }
      for (i in approximate_rows) {
        feature <- feature_table[i, , drop = FALSE]
        id <- as.character(feature$common_feature_id)
        if (id %in% c(dna_matched, protein_exact_matched)) next
        hit <- if (mode == "auto") {
          ggchord_match_common_dna_near_exact(
            target, accver, feature, db$metadata
          )
        } else {
          ggchord_match_common_protein(
            target, accver, feature, frames, TRUE,
            min_protein_identity, min_protein_coverage, db$metadata
          )
        }
        if (length(hit)) candidates <- c(candidates, hit)
      }
    }
  }
  if (!length(candidates)) return(ggchord_empty_common_features())
  out <- ggchord_rbind_fill(candidates)
  if (resolve == "best") {
    lengths <- stats::setNames(nchar(input$sequences), input$ids)
    out <- ggchord_resolve_common_candidates(out, lengths)
  }
  out <- out[order(out$accver, out$start, out$end, out$type,
    out$common_feature_id, out$match_method), , drop = FALSE]
  rownames(out) <- NULL
  out
}

ggchord_expand_feature_segments <- function(data) {
  if (!"segments" %in% names(data) || !is.list(data$segments) || !nrow(data)) {
    return(data)
  }
  rows <- list()
  for (i in seq_len(nrow(data))) {
    segments <- data$segments[[i]]
    if (!is.data.frame(segments) || !nrow(segments)) {
      row <- data[i, , drop = FALSE]
      row$.biological_source_row <- i
      rows[[length(rows) + 1L]] <- row
      next
    }
    segments$.lo <- pmin(segments$start, segments$end)
    segments$.hi <- pmax(segments$start, segments$end)
    segments <- segments[order(segments$.lo, segments$.hi,
      segments$segment_index), , drop = FALSE]
    standard <- segments$segment_type != "gap"
    if (!any(standard)) next

    # A segmented database record describes one biological feature. Adjacent
    # standard segments therefore share one outline and one arrowhead; their
    # joins are retained as internal boundaries. Explicit gaps, coordinate
    # gaps, and colour changes start a new visible feature run.
    runs <- list()
    current <- integer()
    previous_hi <- NA_real_
    previous_colour <- NA_character_
    for (j in seq_len(nrow(segments))) {
      if (!standard[j]) {
        if (length(current)) runs[[length(runs) + 1L]] <- current
        current <- integer()
        previous_hi <- NA_real_
        previous_colour <- NA_character_
        next
      }
      colour <- if ("color" %in% names(segments)) {
        as.character(segments$color[j])
      } else NA_character_
      colour_changed <- length(current) && !is.na(colour) &&
        !is.na(previous_colour) && !identical(colour, previous_colour)
      disconnected <- length(current) && segments$.lo[j] > previous_hi + 1
      if (length(current) && (colour_changed || disconnected)) {
        runs[[length(runs) + 1L]] <- current
        current <- integer()
      }
      current <- c(current, j)
      previous_hi <- if (length(current) == 1L) segments$.hi[j] else
        max(previous_hi, segments$.hi[j])
      previous_colour <- colour
    }
    if (length(current)) runs[[length(runs) + 1L]] <- current

    for (j in seq_along(runs)) {
      run <- segments[runs[[j]], , drop = FALSE]
      row <- data[i, , drop = FALSE]
      row$start <- min(run$.lo)
      row$end <- max(run$.hi)
      row$.biological_source_row <- i
      row$.segment_index <- min(run$segment_index)
      boundaries <- sort(unique(run$.lo[-1L]))
      boundaries <- boundaries[boundaries > row$start & boundaries < row$end]
      row$.feature_boundaries <- I(list(as.numeric(boundaries)))
      boundary_styles <- rep(NA_character_, length(boundaries))
      if (length(boundaries) && "line_style" %in% names(run)) {
        style_rows <- match(boundaries, run$.lo)
        boundary_styles <- as.character(run$line_style[style_rows])
      }
      row$.feature_boundary_styles <- I(list(boundary_styles))
      if ("color" %in% names(run)) {
        run_colour <- run$color[!is.na(run$color) & nzchar(run$color)][1L]
        if (length(run_colour)) row$feature_color <- run_colour
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  out <- ggchord_rbind_fill(rows)
  rownames(out) <- NULL
  out
}
