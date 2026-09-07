# Feature placement is resolved in the sequence-local frame before polygon
# generation.  The ggplot2 Position objects are deliberately identity at the
# Cartesian stage: their fields are consumed by ggchord's layout engine.

PositionGgchordFeature <- ggplot2::ggproto(
  "PositionGgchordFeature", ggplot2::PositionIdentity,
  ggchord_feature_position = TRUE,
  placement = "identity",
  offset = 0
)

#' Place features on opposite sides according to strand
#'
#' Positive signed offsets point away from the plot centre and negative
#' offsets point inward.  A scalar is expanded to `+ = offset` and
#' `- = -offset`; explicitly strand-named values are already signed and are
#' not flipped again.
#'
#' @param offset A finite scalar or a flexible per-sequence/per-strand
#'   specification.
#' @return A ggplot2 Position object.
#' @export
position_strand <- function(offset = 0.1) {
  ggplot2::ggproto(
    NULL, PositionGgchordFeature,
    placement = "strand", offset = offset
  )
}

#' Place both feature strands on one shared band
#'
#' @param offset Signed local-normal offset. Negative values place the shared
#'   band inside a circular sequence and positive values place it outside.
#' @return A ggplot2 Position object.
#' @export
position_plasmid <- function(offset = -0.1) {
  ggplot2::ggproto(
    NULL, PositionGgchordFeature,
    placement = "plasmid", offset = offset
  )
}

ggchord_as_feature_position <- function(position, caller = "feature layer") {
  if (is.character(position)) {
    if (length(position) != 1L || is.na(position) ||
        !position %in% c("identity", "strand", "plasmid")) {
      ggchord_stop(
        caller, ": position must be 'identity', 'strand', 'plasmid', ",
        "position_strand(), position_plasmid(), or position_feature_stack()"
      )
    }
    return(switch(
      position,
      identity = ggplot2::position_identity(),
      strand = position_strand(),
      plasmid = position_plasmid()
    ))
  }
  if (inherits(position, "PositionIdentity") ||
      isTRUE(position$ggchord_feature_position) ||
      isTRUE(position$ggchord_feature_stack)) return(position)
  ggchord_stop(
    caller, ": only identity, strand, plasmid, and feature-stack Positions ",
    "operate in the sequence-local frame"
  )
}

ggchord_is_identity_position <- function(position) {
  inherits(position, "PositionIdentity") &&
    !isTRUE(position$ggchord_feature_position) &&
    !isTRUE(position$ggchord_feature_stack)
}

ggchord_validate_offset_value <- function(x, caller) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    ggchord_stop(caller, ": offsets must be finite numeric scalars")
  }
  as.numeric(x)
}

# Resolve the same sequence/list forms accepted by process_gene_param(), while
# retaining whether a strand value was explicit (and therefore already signed).
ggchord_resolve_position_offset <- function(offset, seqs, placement,
                                             caller = "position") {
  n <- length(seqs)
  expand_scalar <- function(x) {
    x <- ggchord_validate_offset_value(x, caller)
    if (identical(placement, "strand")) c("+" = x, "-" = -x)
    else c("+" = x, "-" = x)
  }
  explicit_strands <- function(x) {
    if (is.null(names(x)) || !length(x) ||
        any(!names(x) %in% c("+", "-")) || anyDuplicated(names(x))) {
      ggchord_stop(caller, ": strand offsets must be uniquely named '+' and/or '-'")
    }
    if (identical(placement, "plasmid")) {
      ggchord_stop(caller, ": position_plasmid() accepts shared offsets, not strand-specific offsets")
    }
    out <- c("+" = 0, "-" = 0)
    for (st in names(x)) out[st] <- ggchord_validate_offset_value(x[[st]], caller)
    out
  }
  resolve_element <- function(x) {
    if (length(x) == 1L && is.null(names(x))) return(expand_scalar(x))
    explicit_strands(x)
  }

  if (!is.list(offset)) {
    nm <- names(offset)
    if (length(offset) == 1L && is.null(nm)) {
      val <- expand_scalar(offset)
      return(stats::setNames(rep(list(val), n), seqs))
    }
    if (!is.null(nm) && all(nm %in% c("+", "-"))) {
      val <- explicit_strands(offset)
      return(stats::setNames(rep(list(val), n), seqs))
    }
    if (!is.null(nm)) {
      unknown <- setdiff(nm, seqs)
      if (length(unknown)) ggchord_stop(caller, ": unknown accver: ", paste(unknown, collapse = ", "))
      out <- stats::setNames(rep(list(expand_scalar(0)), n), seqs)
      for (id in nm) out[[id]] <- expand_scalar(offset[[id]])
      return(out)
    }
    if (length(offset) != n) {
      ggchord_stop(caller, ": unnamed offsets must have length 1 or match the sequence count")
    }
    return(stats::setNames(lapply(offset, expand_scalar), seqs))
  }

  nm <- names(offset)
  if (length(offset) == 1L && (is.null(nm) || !nzchar(nm[1L]))) {
    val <- resolve_element(offset[[1L]])
    return(stats::setNames(rep(list(val), n), seqs))
  }
  if (!is.null(nm) && length(nm) && all(nm %in% seqs)) {
    out <- stats::setNames(rep(list(expand_scalar(0)), n), seqs)
    for (id in nm) out[[id]] <- resolve_element(offset[[id]])
    return(out)
  }
  idx_names <- as.character(seq_along(seqs))
  if (!is.null(nm) && length(nm) && all(nm %in% idx_names)) {
    out <- stats::setNames(rep(list(expand_scalar(0)), n), seqs)
    for (i in seq_along(nm)) out[[seqs[as.integer(nm[i])]]] <- resolve_element(offset[[i]])
    return(out)
  }
  if (is.null(nm) || all(!nzchar(nm))) {
    if (length(offset) != n) ggchord_stop(caller, ": unnamed offset lists must match the sequence count")
    return(stats::setNames(lapply(offset, resolve_element), seqs))
  }
  ggchord_stop(caller, ": invalid flexible offset specification")
}

ggchord_resolve_legacy_offset <- function(offset, seqs, name) {
  old <- process_gene_param(offset, seqs, name, 0, FALSE)
  lapply(old, function(x) c("+" = unname(x["+"]), "-" = -unname(x["-"])))
}

ggchord_allocate_feature_lanes <- function(data, base, direction, spacing,
                                            circular = FALSE, lengths = NULL) {
  lane <- integer(nrow(data))
  band <- paste(
    data$accver,
    format(base, digits = 15, scientific = FALSE, trim = TRUE),
    direction,
    sep = "\r"
  )
  groups <- split(seq_len(nrow(data)), factor(band, levels = unique(band)))
  for (idx in groups) {
    sid <- as.character(data$accver[idx[1L]])
    len <- if (!is.null(lengths)) unname(lengths[sid]) else NA_real_
    intervals <- lapply(idx, function(i) {
      a <- as.numeric(data$start[i]); b <- as.numeric(data$end[i])
      if (isTRUE(circular) && is.finite(len) && a > b) {
        rbind(c(a, len), c(1, b))
      } else matrix(c(min(a, b), max(a, b)), nrow = 1L)
    })
    ord <- order(vapply(intervals, function(x) min(x[, 1]), numeric(1)),
                 vapply(intervals, function(x) max(x[, 2]), numeric(1)), idx)
    lane_intervals <- list()
    overlaps <- function(a, occupied) {
      any(vapply(occupied, function(b) any(
        outer(seq_len(nrow(a)), seq_len(nrow(b)), Vectorize(function(i, j) {
          a[i, 1] <= b[j, 2] && b[j, 1] <= a[i, 2]
        }))
      ), logical(1)))
    }
    for (local in ord) {
      chosen <- 1L
      while (chosen <= length(lane_intervals) &&
             overlaps(intervals[[local]], lane_intervals[[chosen]])) chosen <- chosen + 1L
      if (chosen > length(lane_intervals)) lane_intervals[[chosen]] <- list()
      lane_intervals[[chosen]][[length(lane_intervals[[chosen]]) + 1L]] <- intervals[[local]]
      lane[idx[local]] <- chosen - 1L
    }
  }
  lane
}

ggchord_apply_feature_position <- function(data, position, seqs, lengths,
                                            circular = FALSE,
                                            legacy_offset = NULL,
                                            legacy_name = "gene_offset") {
  if (is.null(data) || !nrow(data)) return(data)
  ggchord_require_columns(
    data, c("accver", "start", "end", "strand"), legacy_name
  )
  unknown <- setdiff(unique(as.character(data$accver)), seqs)
  if (length(unknown)) {
    ggchord_stop(legacy_name, ": unknown accver: ", paste(unknown, collapse = ", "))
  }
  if (anyNA(data$strand) || any(!as.character(data$strand) %in% c("+", "-"))) {
    ggchord_stop(legacy_name, ": strand must contain only '+' or '-'")
  }
  if (!is.numeric(data$start) || !is.numeric(data$end) ||
      any(!is.finite(data$start)) || any(!is.finite(data$end))) {
    ggchord_stop(legacy_name, ": start and end must be finite numeric values")
  }
  position <- position %||% ggplot2::position_identity()
  if (!is.null(legacy_offset)) {
    resolved <- ggchord_resolve_legacy_offset(legacy_offset, seqs, legacy_name)
    position_name <- "legacy_strand"
  } else if (isTRUE(position$ggchord_feature_stack)) {
    base_position <- position$base_position
    if (is.null(base_position)) {
      side <- position$side %||% "strand"
      resolved <- switch(
        side,
        strand = ggchord_resolve_position_offset(0.1, seqs, "strand"),
        outside = ggchord_resolve_position_offset(-0.1, seqs, "plasmid"),
        inside = ggchord_resolve_position_offset(0.1, seqs, "plasmid")
      )
      position_name <- "feature_stack_legacy"
    } else {
      resolved <- ggchord_resolve_position_offset(
        base_position$offset %||% 0, seqs,
        base_position$placement %||% "plasmid", "position_feature_stack(base_position)"
      )
      position_name <- paste0("feature_stack_", base_position$placement %||% "identity")
    }
  } else if (isTRUE(position$ggchord_feature_position)) {
    resolved <- ggchord_resolve_position_offset(
      position$offset, seqs, position$placement,
      paste0("position_", position$placement, "()")
    )
    position_name <- position$placement
  } else {
    resolved <- stats::setNames(rep(list(c("+" = 0, "-" = 0)), length(seqs)), seqs)
    position_name <- "identity"
  }

  out <- as.data.frame(data, stringsAsFactors = FALSE)
  base <- vapply(seq_len(nrow(out)), function(i) {
    unname(resolved[[as.character(out$accver[i])]][as.character(out$strand[i])])
  }, numeric(1))
  lane <- integer(nrow(out)); lane_offset <- numeric(nrow(out))
  if (isTRUE(position$ggchord_feature_stack)) {
    direction <- ifelse(base < 0, -1, 1)
    lane <- ggchord_allocate_feature_lanes(
      out, base, direction, position$spacing %||% 0.08,
      circular = circular, lengths = lengths
    )
    lane_offset <- direction * lane * (position$spacing %||% 0.08)
  }
  out$.position_name <- position_name
  out$.position_base_offset <- base
  out$.feature_stack_lane <- lane
  out$.position_lane_offset <- lane_offset
  out$.normal_offset <- base + lane_offset
  out
}
