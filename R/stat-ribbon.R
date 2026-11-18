# Formal ribbon statistics. Computation happens before shared chord geometry;
# the ggproto objects then expose the computed columns to after_stat().

StatGgchordRibbon <- ggplot2::ggproto(
  "StatGgchordRibbon", ggplot2::Stat,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(
    .ggchord_bundle_n = NA_real_,
    .ggchord_bundle_weight = NA_real_,
    .ggchord_density = NA_real_
  ),
  compute_panel = function(data, scales) {
    data$bundle_n <- data$.ggchord_bundle_n
    data$bundle_weight <- data$.ggchord_bundle_weight
    data$density <- data$.ggchord_density
    data
  }
)

StatGgchordRibbonBundle <- ggplot2::ggproto(
  "StatGgchordRibbonBundle", StatGgchordRibbon,
  source_type = "bundle",
  source_parameters = c("bins", "min_bundle", "weight", "group_by"),
  compute_source = function(ribbon_data, seq_data, ...) {
    result <- bundle_ggchord_ribbons(ribbon_data, seq_data, ...)
    result$data$bundle_n <- result$data$.bundle_n
    result$data$bundle_weight <- result$data$.bundle_weight
    result$data$density <- result$data$.bundle_density
    result
  }
)

StatGgchordRibbonDensity <- ggplot2::ggproto(
  "StatGgchordRibbonDensity", StatGgchordRibbon,
  source_type = "density",
  source_parameters = c("bins", "weight", "group_by"),
  compute_source = function(ribbon_data, seq_data, bins = 80L, weight = "length", group_by = NULL) {
    ggchord_ribbon_density(ribbon_data, seq_data, bins, match.arg(weight, c("length", "count", "pident")),
      group_by, "stat_ribbon_density()")
  }
)

ggchord_ribbon_density <- function(ribbon_data, seq_data, bins, weight,
                                   group_by, caller) {
  checked <- ggchord_check_ribbon_tables(seq_data, ribbon_data, caller)
  if (!is.numeric(bins) || length(bins) != 1L || !is.finite(bins) ||
      bins < 1 || bins != as.integer(bins)) {
    ggchord_stop(caller, ": bins must be one positive integer")
  }
  bins <- as.integer(bins)
  if (!is.null(group_by) &&
      (!is.character(group_by) || anyNA(group_by) ||
       any(!group_by %in% names(ribbon_data)))) {
    ggchord_stop(caller, ": group_by must name columns in ribbon_data")
  }
  input <- as.data.frame(ribbon_data, stringsAsFactors = FALSE)
  if (nrow(input) == 0L) {
    input$.bundle_n <- integer(0)
    input$.bundle_weight <- numeric(0)
    input$.bundle_density <- numeric(0)
    input$bundle_n <- integer(0)
    input$bundle_weight <- numeric(0)
    input$density <- numeric(0)
    return(list(data = input, report = data.frame()))
  }
  lengths <- checked$lengths
  qid <- as.character(input$qaccver)
  sid <- as.character(input$saccver)
  qmid <- (input$qstart + input$qend) / 2
  smid <- (input$sstart + input$send) / 2
  qbin <- pmin(bins, pmax(1L, floor((qmid - 1) / lengths[qid] * bins) + 1L))
  sbin <- pmin(bins, pmax(1L, floor((smid - 1) / lengths[sid] * bins) + 1L))
  direction <- ggchord_ribbon_direction(input)
  parts <- list(qid, sid, direction, qbin, sbin)
  if (length(group_by)) parts <- c(parts, lapply(input[group_by], as.character))
  key <- do.call(paste, c(parts, sep = "\r"))
  row_weight <- ggchord_ribbon_weights(input, weight)
  group_n <- stats::ave(rep(1L, nrow(input)), key, FUN = sum)
  group_weight <- stats::ave(row_weight, key, FUN = sum)
  pair <- paste(qid, sid, sep = "\r")
  pair_max <- stats::ave(group_weight, pair, FUN = max)
  input$.bundle_n <- as.integer(group_n)
  input$.bundle_weight <- as.numeric(group_weight)
  input$.bundle_density <- ifelse(pair_max > 0, group_weight / pair_max, 0)
  input$bundle_n <- input$.bundle_n
  input$bundle_weight <- input$.bundle_weight
  input$density <- input$.bundle_density
  list(
    data = input,
    report = data.frame(
      source_row = seq_len(nrow(input)),
      bin_count = input$.bundle_n,
      bin_weight = input$.bundle_weight,
      density = input$.bundle_density,
      stringsAsFactors = FALSE
    )
  )
}

ggchord_ribbon_stat <- function(name, dots) {
  stats <- list(ribbon_bundle = StatGgchordRibbonBundle, ribbon_density = StatGgchordRibbonDensity)
  if (!is.character(name) || length(name) != 1L || !name %in% c("identity", names(stats)))
    ggchord_stop("Choose ribbon stat from: identity, ", paste(names(stats), collapse = ", "))
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) || anyDuplicated(names(dots))))
    ggchord_stop("Supply each ribbon geometry or statistic parameter once, by name")
  stat <- stats[[name]]
  controls <- unique(unlist(lapply(stats, function(s) s$source_parameters)))
  arguments <- dots[intersect(names(dots), controls)]
  if (length(setdiff(names(arguments), stat$source_parameters)))
    ggchord_stop("Unsupported parameter for ribbon stat '", name, "': ",
      paste(setdiff(names(arguments), stat$source_parameters), collapse = ", "))
  list(stat = stat, arguments = arguments, geom_parameters = dots[!names(dots) %in% controls])
}

ggchord_apply_ribbon_stat <- function(lyr, specification) {
  # ggplot2 deliberately drops input columns that are not mapped before a Stat
  # runs. Carry the three precomputed values through private aesthetics, then
  # expose their public names from compute_panel() for after_stat().
  lyr$mapping[[".ggchord_bundle_n"]] <- as.name("bundle_n")
  lyr$mapping[[".ggchord_bundle_weight"]] <- as.name("bundle_weight")
  lyr$mapping[[".ggchord_density"]] <- as.name("density")
  lyr$stat <- specification$stat
  lyr$ggchord_params$ribbon_stat <- list(type = specification$stat$source_type,
    compute = specification$stat$compute_source, arguments = specification$arguments)
  lyr
}

#' Bundle ribbons as a ggchord stat
#'
#' A layer form of [bundle_ggchord_ribbons()] that computes bundles before the
#' shared chord geometry and exposes `after_stat(bundle_n)`,
#' `after_stat(bundle_weight)`, and `after_stat(density)`.
#'
#' @param mapping,data Standard layer mapping and data arguments.
#' @param bins,min_bundle,weight,group_by Passed to
#'   [bundle_ggchord_ribbons()].
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Ribbon geometry and appearance arguments passed to
#'   [geom_link_ribbon()].
#' @return A ggchord ribbon layer.
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
stat_ribbon_bundle <- function(
    mapping = NULL, data = NULL,
    bins = 80L, min_bundle = 2L,
    weight = c("length", "count", "pident"),
    group_by = NULL, position = "identity",
    show.legend = TRUE, inherit.aes = TRUE, ...) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_link_ribbon(mapping = mapping, data = data, stat = "ribbon_bundle",
    bins = bins, min_bundle = min_bundle, weight = match.arg(weight), group_by = group_by,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes, ...)
}

#' Compute local ribbon density as a ggchord stat
#'
#' Keeps one output ribbon per input row while measuring density in normalized
#' query/subject midpoint bins. The computed variables are the same as
#' [stat_ribbon_bundle()], so density can be mapped with
#' `aes(ribbon_alpha = after_stat(density))`.
#'
#' @inheritParams stat_ribbon_bundle
#' @param bins Positive integer number of normalized midpoint bins per
#'   sequence.
#' @param weight Density weight: `"length"`, `"count"`, or `"pident"`.
#' @param group_by Optional columns that define independent density groups.
#' @return A ggchord ribbon layer.
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
stat_ribbon_density <- function(
    mapping = NULL, data = NULL,
    bins = 80L,
    weight = c("length", "count", "pident"),
    group_by = NULL, position = "identity",
    show.legend = TRUE, inherit.aes = TRUE, ...) {
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_link_ribbon(mapping = mapping, data = data, stat = "ribbon_density",
    bins = bins, weight = match.arg(weight), group_by = group_by,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes, ...)
}
