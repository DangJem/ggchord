# Explicit helpers for dense ribbon data and deterministic sequence layout.

ggchord_ribbon_required_columns <- function() {
  c("qaccver", "saccver", "length", "pident",
    "qstart", "qend", "sstart", "send")
}

ggchord_check_ribbon_tables <- function(seq_data, ribbon_data, caller) {
  if (!is.data.frame(seq_data) ||
      !all(c("seq_id", "length") %in% names(seq_data))) {
    ggchord_stop(caller, ": seq_data must contain seq_id and length")
  }
  if (!is.data.frame(ribbon_data)) {
    ggchord_stop(caller, ": ribbon_data must be a data.frame")
  }
  required <- ggchord_ribbon_required_columns()
  missing <- setdiff(required, names(ribbon_data))
  if (length(missing)) {
    ggchord_stop(
      caller, ": ribbon_data is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  seq_ids <- as.character(seq_data$seq_id)
  seq_lengths <- seq_data$length
  if (anyNA(seq_ids) || any(!nzchar(seq_ids)) || anyDuplicated(seq_ids)) {
    ggchord_stop(caller, ": seq_data$seq_id must be unique and non-missing")
  }
  if (!is.numeric(seq_lengths) || any(!is.finite(seq_lengths)) ||
      any(seq_lengths <= 0)) {
    ggchord_stop(caller, ": seq_data$length must contain positive numbers")
  }

  numeric_columns <- c("length", "pident", "qstart", "qend", "sstart", "send")
  if (any(!vapply(ribbon_data[numeric_columns], is.numeric, logical(1)))) {
    ggchord_stop(caller, ": ribbon coordinate, length and pident columns must be numeric")
  }
  if (any(!is.finite(as.matrix(ribbon_data[numeric_columns])))) {
    ggchord_stop(caller, ": ribbon coordinate, length and pident values must be finite")
  }
  unknown <- setdiff(
    unique(c(as.character(ribbon_data$qaccver),
             as.character(ribbon_data$saccver))),
    seq_ids
  )
  if (length(unknown)) {
    ggchord_stop(
      caller, ": ribbon_data contains sequence IDs absent from seq_data: ",
      paste(unknown, collapse = ", ")
    )
  }

  list(
    seq_ids = seq_ids,
    lengths = stats::setNames(as.numeric(seq_lengths), seq_ids),
    required = required
  )
}

ggchord_ribbon_weights <- function(ribbon_data, weight) {
  switch(weight,
    count = rep(1, nrow(ribbon_data)),
    length = pmax(0, as.numeric(ribbon_data$length)),
    pident = pmax(0, as.numeric(ribbon_data$length)) *
      pmax(0, as.numeric(ribbon_data$pident)) / 100
  )
}

ggchord_typed_na <- function(x) {
  if (is.factor(x)) return(factor(NA_character_, levels = levels(x)))
  if (inherits(x, "Date")) return(as.Date(NA))
  if (inherits(x, "POSIXct")) return(as.POSIXct(NA_real_, origin = "1970-01-01"))
  x[NA_integer_][1]
}

#' Bundle dense alignment ribbons explicitly
#'
#' Aggregates nearby alignment ribbons on a normalized query/subject grid.
#' Bundling is always explicit: \code{geom_ribbon()} continues to draw one
#' ribbon per input row unless the returned data are supplied by the user.
#'
#' @param ribbon_data Alignment data in ggchord ribbon format.
#' @param seq_data Sequence data containing \code{seq_id} and \code{length}.
#' @param bins Positive integer number of normalized midpoint bins per
#'   sequence, default 80.
#' @param min_bundle Minimum number of rows required for aggregation, default
#'   2. Smaller groups remain as individual ribbons.
#' @param weight Bundle weight: \code{"length"}, \code{"count"}, or
#'   \code{"pident"}. Identity weight is aligned length multiplied by identity
#'   proportion.
#' @param group_by Optional character columns that must agree before ribbons
#'   can be bundled.
#'
#' @return A list with \code{data} and \code{report}. The data include
#'   \code{.bundle_n}, \code{.bundle_weight}, and \code{.bundle_density}; the
#'   report maps each output row back to its source rows.
#' @export
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' bundled <- bundle_ggchord_ribbons(ribbon_data_example, seq_data_example)
#' nrow(bundled$data)
#' bundled_plot <- ggchord(seq_data_example, bundled$data) +
#'   geom_seq() +
#'   geom_ribbon(ggplot2::aes(ribbon_alpha = .bundle_density)) +
#'   scale_ribbon_alpha_continuous()
bundle_ggchord_ribbons <- function(
    ribbon_data,
    seq_data,
    bins = 80L,
    min_bundle = 2L,
    weight = c("length", "count", "pident"),
    group_by = NULL) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  caller <- "bundle_ggchord_ribbons()"
  checked <- ggchord_check_ribbon_tables(seq_data, ribbon_data, caller)
  weight <- match.arg(weight)
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
      bins < 1 || bins != as.integer(bins)) {
    ggchord_stop(caller, ": bins must be one positive integer")
  }
  if (!is.numeric(min_bundle) || length(min_bundle) != 1L ||
      !is.finite(min_bundle) || min_bundle < 2 ||
      min_bundle != as.integer(min_bundle)) {
    ggchord_stop(caller, ": min_bundle must be an integer of at least 2")
  }
  bins <- as.integer(bins)
  min_bundle <- as.integer(min_bundle)
  if (!is.null(group_by)) {
    if (!is.character(group_by) || anyNA(group_by) ||
        any(!group_by %in% names(ribbon_data))) {
      ggchord_stop(caller, ": group_by must name columns in ribbon_data")
    }
    group_by <- unique(group_by)
  }
  reserved <- intersect(
    names(ribbon_data), c(".bundle_n", ".bundle_weight", ".bundle_density")
  )
  if (length(reserved)) {
    ggchord_stop(
      caller, ": input already contains reserved bundle column(s): ",
      paste(reserved, collapse = ", ")
    )
  }

  input <- as.data.frame(ribbon_data, stringsAsFactors = FALSE)
  n_input <- nrow(input)
  if (n_input == 0L) {
    input$.bundle_n <- integer(0)
    input$.bundle_weight <- numeric(0)
    input$.bundle_density <- numeric(0)
    return(list(
      data = input,
      report = data.frame(
        output_row = integer(0), source_rows = character(0),
        n_input = integer(0), bundled = logical(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  lengths <- checked$lengths
  qid <- as.character(input$qaccver)
  sid <- as.character(input$saccver)
  qmid <- (input$qstart + input$qend) / 2
  smid <- (input$sstart + input$send) / 2
  qbin <- pmin(bins, pmax(1L, floor((qmid - 1) / lengths[qid] * bins) + 1L))
  sbin <- pmin(bins, pmax(1L, floor((smid - 1) / lengths[sid] * bins) + 1L))
  qdir <- sign(input$qend - input$qstart)
  sdir <- sign(input$send - input$sstart)
  direction <- ifelse(qdir == sdir, "same", "opposite")
  row_weight <- ggchord_ribbon_weights(input, weight)

  key_parts <- list(qid, sid, direction, qbin, sbin)
  if (length(group_by)) {
    key_parts <- c(key_parts, lapply(input[group_by], as.character))
  }
  keys <- do.call(paste, c(key_parts, sep = "\r"))
  groups <- split(seq_len(n_input), factor(keys, levels = unique(keys)))

  required <- checked$required
  extra <- setdiff(names(input), required)
  out_rows <- list()
  source_map <- list()

  append_single <- function(i) {
    row <- input[i, , drop = FALSE]
    row$.bundle_n <- 1L
    row$.bundle_weight <- row_weight[i]
    out_rows[[length(out_rows) + 1L]] <<- row
    source_map[[length(source_map) + 1L]] <<- i
  }

  append_bundle <- function(idx) {
    row <- input[idx[1], , drop = FALSE]
    local_weight <- row_weight[idx]
    mean_weight <- pmax(as.numeric(input$length[idx]), .Machine$double.eps)
    qlo <- min(input$qstart[idx], input$qend[idx])
    qhi <- max(input$qstart[idx], input$qend[idx])
    slo <- min(input$sstart[idx], input$send[idx])
    shi <- max(input$sstart[idx], input$send[idx])
    row$qstart <- if (qdir[idx[1]] >= 0) qlo else qhi
    row$qend <- if (qdir[idx[1]] >= 0) qhi else qlo
    row$sstart <- if (sdir[idx[1]] >= 0) slo else shi
    row$send <- if (sdir[idx[1]] >= 0) shi else slo
    row$length <- round(stats::weighted.mean(input$length[idx], mean_weight))
    row$pident <- stats::weighted.mean(input$pident[idx], mean_weight)

    for (nm in extra) {
      values <- input[[nm]][idx]
      equal <- if (length(values) <= 1L) TRUE else {
        all(vapply(seq_along(values)[-1L], function(j) {
          identical(values[j], values[1])
        }, logical(1)))
      }
      if (!equal) row[[nm]] <- ggchord_typed_na(input[[nm]])
    }
    row$.bundle_n <- length(idx)
    row$.bundle_weight <- sum(local_weight)
    out_rows[[length(out_rows) + 1L]] <<- row
    source_map[[length(source_map) + 1L]] <<- idx
  }

  for (idx in groups) {
    if (length(idx) < min_bundle) {
      for (i in idx) append_single(i)
    } else {
      append_bundle(idx)
    }
  }

  first_source <- vapply(source_map, min, integer(1))
  ord <- order(first_source)
  out_rows <- out_rows[ord]
  source_map <- source_map[ord]
  out <- ggchord_rbind_fill(out_rows)

  pair <- paste(out$qaccver, out$saccver, sep = "\r")
  pair_max <- stats::ave(out$.bundle_weight, pair, FUN = max)
  out$.bundle_density <- ifelse(pair_max > 0, out$.bundle_weight / pair_max, 0)
  attr(out, "source_rows") <- unlist(source_map, use.names = FALSE)
  report <- data.frame(
    output_row = seq_len(nrow(out)),
    source_rows = vapply(source_map, paste, collapse = ",", character(1)),
    n_input = lengths(source_map),
    bundled = lengths(source_map) >= min_bundle,
    stringsAsFactors = FALSE
  )

  list(data = out, report = report)
}

ggchord_resolve_orientation <- function(seq_orientation, seq_ids, caller) {
  if (is.null(seq_orientation)) {
    return(stats::setNames(rep(1, length(seq_ids)), seq_ids))
  }
  if (!is.numeric(seq_orientation) || anyNA(seq_orientation)) {
    ggchord_stop(caller, ": seq_orientation must contain only 1 and -1")
  }
  if (length(seq_orientation) == 1L) {
    out <- rep(seq_orientation, length(seq_ids))
  } else if (!is.null(names(seq_orientation))) {
    if (!all(seq_ids %in% names(seq_orientation))) {
      ggchord_stop(caller, ": named seq_orientation must include every sequence")
    }
    out <- seq_orientation[seq_ids]
  } else if (length(seq_orientation) == length(seq_ids)) {
    out <- seq_orientation
  } else {
    ggchord_stop(caller, ": seq_orientation must have length 1 or match seq_data")
  }
  if (any(!out %in% c(-1, 1))) {
    ggchord_stop(caller, ": seq_orientation must contain only 1 and -1")
  }
  stats::setNames(as.numeric(out), seq_ids)
}

ggchord_optimizer_rows <- function(ribbon_data, lengths, weights, approximate) {
  qid <- as.character(ribbon_data$qaccver)
  sid <- as.character(ribbon_data$saccver)
  qfrac <- ((ribbon_data$qstart + ribbon_data$qend) / 2 - 1) / lengths[qid]
  sfrac <- ((ribbon_data$sstart + ribbon_data$send) / 2 - 1) / lengths[sid]
  qfrac <- pmin(1, pmax(0, qfrac))
  sfrac <- pmin(1, pmax(0, sfrac))
  scored <- data.frame(qid = qid, sid = sid, qfrac = qfrac, sfrac = sfrac,
                       weight = weights, stringsAsFactors = FALSE)
  if (!isTRUE(approximate) || nrow(scored) == 0L) return(scored)

  qbin <- pmin(64L, floor(scored$qfrac * 64L) + 1L)
  sbin <- pmin(64L, floor(scored$sfrac * 64L) + 1L)
  key <- paste(scored$qid, scored$sid, qbin, sbin, sep = "\r")
  groups <- split(seq_len(nrow(scored)), factor(key, levels = unique(key)))
  do.call(rbind, lapply(groups, function(idx) {
    w <- scored$weight[idx]
    if (sum(w) <= 0) w <- rep(1, length(idx))
    data.frame(
      qid = scored$qid[idx[1]], sid = scored$sid[idx[1]],
      qfrac = stats::weighted.mean(scored$qfrac[idx], w),
      sfrac = stats::weighted.mean(scored$sfrac[idx], w),
      weight = sum(scored$weight[idx]), stringsAsFactors = FALSE
    )
  }))
}

ggchord_layout_score <- function(scored, seq_order, orientation, lengths) {
  if (nrow(scored) == 0L) return(0)
  ordered_lengths <- lengths[seq_order]
  starts <- c(0, utils::head(cumsum(ordered_lengths), -1)) /
    sum(ordered_lengths)
  names(starts) <- seq_order
  widths <- ordered_lengths / sum(ordered_lengths)
  names(widths) <- seq_order
  qf <- ifelse(orientation[scored$qid] == 1,
                scored$qfrac, 1 - scored$qfrac)
  sf <- ifelse(orientation[scored$sid] == 1,
                scored$sfrac, 1 - scored$sfrac)
  p1 <- starts[scored$qid] + qf * widths[scored$qid]
  p2 <- starts[scored$sid] + sf * widths[scored$sid]
  lo <- pmin(p1, p2)
  hi <- pmax(p1, p2)
  w <- pmax(0, scored$weight)
  if (sum(w) <= 0) w <- rep(1, length(w))
  w <- w / sum(w)

  crossing <- 0
  if (length(lo) > 1L) {
    for (i in seq_len(length(lo) - 1L)) {
      j <- seq.int(i + 1L, length(lo))
      inside_lo <- lo[j] > lo[i] & lo[j] < hi[i]
      inside_hi <- hi[j] > lo[i] & hi[j] < hi[i]
      crossed <- xor(inside_lo, inside_hi)
      if (any(crossed)) crossing <- crossing + sum(w[i] * w[j[crossed]])
    }
  }
  span <- pmin(hi - lo, 1 - (hi - lo))
  crossing + 0.05 * sum(w * span)
}

#' Optimize sequence order and orientation for ribbon readability
#'
#' Searches deterministic sequence-order and orientation changes that reduce
#' weighted ribbon centre-line crossings and long chord spans. It never
#' modifies the supplied data and never accepts a worse arrangement. Inputs
#' above 2,000 ribbons use a deterministic 64-bin approximation and record
#' \code{approximate = TRUE} in the returned report.
#'
#' @param seq_data Sequence data containing \code{seq_id} and \code{length}.
#' @param ribbon_data Alignment data in ggchord ribbon format.
#' @param seq_order Optional initial complete sequence order. The default uses
#'   the row order of \code{seq_data}.
#' @param seq_orientation Optional initial orientation (1 or -1), scalar,
#'   complete vector, or named vector.
#' @param optimize Any of \code{"order"} and \code{"orientation"}.
#' @param weight Objective weight: \code{"length"}, \code{"count"}, or
#'   \code{"pident"}.
#' @param max_iter Positive integer maximum number of improving search rounds.
#'
#' @return An object of class \code{ggchord_layout_optimization} containing
#'   \code{seq_order}, \code{seq_orientation}, before/after scores and a report.
#' @export
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' optimized <- optimize_ggchord_layout(seq_data_example, ribbon_data_example)
#' optimized$seq_order
#' optimized_plot <- ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq(
#'     seq_order = optimized$seq_order,
#'     seq_orientation = optimized$seq_orientation
#'   ) +
#'   geom_ribbon()
optimize_ggchord_layout <- function(
    seq_data,
    ribbon_data,
    seq_order = NULL,
    seq_orientation = NULL,
    optimize = c("order", "orientation"),
    weight = c("length", "count", "pident"),
    max_iter = 50L) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  caller <- "optimize_ggchord_layout()"
  checked <- ggchord_check_ribbon_tables(seq_data, ribbon_data, caller)
  weight <- match.arg(weight)
  allowed <- c("order", "orientation")
  if (!is.character(optimize) || anyNA(optimize) ||
      any(!optimize %in% allowed)) {
    ggchord_stop(caller, ": optimize may contain only 'order' and 'orientation'")
  }
  optimize <- unique(optimize)
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !is.finite(max_iter) ||
      max_iter < 1 || max_iter != as.integer(max_iter)) {
    ggchord_stop(caller, ": max_iter must be one positive integer")
  }
  max_iter <- as.integer(max_iter)

  original_order <- checked$seq_ids
  if (is.null(seq_order)) seq_order <- original_order
  if (!is.character(seq_order) || anyNA(seq_order) ||
      length(seq_order) != length(original_order) ||
      !setequal(seq_order, original_order) || anyDuplicated(seq_order)) {
    ggchord_stop(caller, ": seq_order must contain every sequence exactly once")
  }
  orientation <- ggchord_resolve_orientation(
    seq_orientation, original_order, caller
  )
  weights <- ggchord_ribbon_weights(ribbon_data, weight)
  approximate <- nrow(ribbon_data) > 2000L
  scored <- ggchord_optimizer_rows(
    ribbon_data, checked$lengths, weights, approximate
  )
  current_order <- seq_order
  current_score <- ggchord_layout_score(
    scored, current_order, orientation, checked$lengths
  )
  before_score <- current_score
  iterations <- 0L
  tolerance <- sqrt(.Machine$double.eps)

  for (iteration in seq_len(max_iter)) {
    candidates <- list()
    if ("order" %in% optimize && length(current_order) > 2L) {
      for (i in seq.int(2L, length(current_order) - 1L)) {
        candidate <- current_order
        candidate[c(i, i + 1L)] <- candidate[c(i + 1L, i)]
        candidates[[length(candidates) + 1L]] <- list(
          order = candidate, orientation = orientation
        )
      }
    }
    if ("orientation" %in% optimize) {
      for (id in original_order) {
        candidate_orientation <- orientation
        candidate_orientation[id] <- -candidate_orientation[id]
        candidates[[length(candidates) + 1L]] <- list(
          order = current_order, orientation = candidate_orientation
        )
      }
    }
    if (!length(candidates)) break

    scores <- vapply(candidates, function(candidate) {
      ggchord_layout_score(
        scored, candidate$order, candidate$orientation, checked$lengths
      )
    }, numeric(1))
    best <- which.min(scores)
    if (!is.finite(scores[best]) || scores[best] >= current_score - tolerance) {
      break
    }
    current_order <- candidates[[best]]$order
    orientation <- candidates[[best]]$orientation
    current_score <- scores[best]
    iterations <- iteration
  }

  report <- list(
    approximate = approximate,
    n_ribbons = nrow(ribbon_data),
    n_scored = nrow(scored),
    iterations = iterations,
    weight = weight,
    optimize = optimize,
    improved = current_score < before_score - tolerance
  )
  out <- list(
    seq_order = current_order,
    seq_orientation = orientation[current_order],
    score_before = before_score,
    score_after = current_score,
    report = report
  )
  class(out) <- c("ggchord_layout_optimization", "list")
  out
}
