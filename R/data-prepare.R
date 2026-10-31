# Explicit preparation composes existing helpers without changing drawing.

ggchord_prepare_options <- function(x, allowed, caller) {
  if (!is.list(x) || (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x)) || any(!names(x) %in% allowed)))) {
    ggchord_stop(caller, ": expected a named list; allowed entries: ", paste(allowed, collapse = ", "))
  }
  x
}

ggchord_prepare_collection <- function(x, mapping, converter, caller) {
  if (is.null(x)) x <- list()
  if (!is.list(x) || is.data.frame(x) || (length(x) &&
      (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))))) {
    ggchord_stop(caller, ": collections must be named lists of input tables")
  }
  if (is.null(mapping)) mapping <- list()
  ggchord_prepare_options(mapping, names(x), caller)
  out <- lapply(names(x), function(name) converter(x[[name]], mapping = mapping[[name]]))
  stats::setNames(out, names(x))
}

ggchord_prepare_feature_issues <- function(features, seq_data) {
  out <- list()
  add <- function(table, category, rows, column, message, severity) {
    if (length(rows)) out[[length(out) + 1L]] <<- data.frame(
      table = table, category = category, row = rows, column = column,
      message = message, severity = severity, stringsAsFactors = FALSE)
  }
  for (name in names(features)) {
    data <- features[[name]]
    index <- match(data$accver, seq_data$accver)
    add(name, "unknown_id", which(is.na(index)), "accver", "Unknown sequence ID", "error")
    bad <- data$start > seq_data$length[index] | data$end > seq_data$length[index]
    add(name, "out_of_range", which(bad), "start/end", "Interval exceeds sequence length", "error")
    if (!"feature_id" %in% names(data)) {
      missing <- seq_len(nrow(data))
    } else missing <- which(is.na(data$feature_id) | !nzchar(as.character(data$feature_id)))
    add(name, "missing_feature_id", missing, "feature_id",
      "No stable source feature ID; ID-based joins require explicit IDs", "warning")
  }
  if (length(out)) do.call(rbind, out) else data.frame(table = character(),
    category = character(), row = integer(), column = character(),
    message = character(), severity = character())
}

ggchord_prepare_direction <- function(data, caller) {
  if (is.null(data) || !"direction" %in% names(data)) return(invisible(NULL))
  endpoint <- sign(data$qend - data$qstart) * sign(data$send - data$sstart)
  direction <- as.character(data$direction)
  conflict <- (direction %in% c("reverse", "opposite") & endpoint > 0) |
    (direction %in% "same" & endpoint < 0)
  if (any(conflict)) ggchord_stop(caller, ": direction contradicts or cannot be represented by endpoint order; ",
    "supply original directional endpoints rather than previously sorted links")
  invisible(NULL)
}

#' Prepare named ggchord tables with an explicit audit trail
#'
#' Converts inputs using the existing `as_ggchord_*()` adapters and returns
#' named feature/link collections alongside the constructor's three tables.
#' Collections are passed explicitly to layers, for example
#' `geom_feature(data = result$features$sites)`; preparation never adds layers.
#' Import files first with the appropriate reader. No file format, alignment,
#' sequence length, genome membership or biological identity is guessed.
#'
#' By default this only converts and validates. Optional operations run in the
#' fixed order `clean`, `filter`, `deduplicate`, `merge`, `bundle`, `optimize`.
#' Each entry in `steps` is a named list of arguments for the corresponding
#' existing helper. An empty list requests its defaults. `clean` affects only
#' the constructor tables; the other operations affect `ribbon_data` only.
#' Preparation sets cleaning's `reversed_interval = "keep"` by default to
#' preserve directional endpoints. Explicit `"sort"` cannot be followed by
#' deduplication, merging, bundling or optimization, which use those endpoints.
#' Those operations also reject previously sorted input whose `direction`
#' contradicts the endpoints; supply the original directional link table.
#' Geometry, named feature/link tables and the original inputs are never
#' silently filtered, merged or reoriented. Optimization returns a suggestion.
#' Expensive pairwise duplicate diagnostics are disabled in validation here;
#' request `steps = list(deduplicate = list(...))` to deduplicate explicitly.
#'
#' An optional `genome` column in `seq_data` groups contigs as metadata.
#' `accver` must remain globally unique; contigs are never concatenated.
#' Stable feature IDs and parent ID lists are preserved when supplied, not
#' manufactured from row numbers. Repeated IDs can describe feature segments;
#' uniqueness for a specific join is checked by the consumer, such as
#' [find_ggchord_loci()].
#'
#' @param seq_data Sequence input accepted by [as_ggchord_seq()].
#' @param ribbon_data Optional paired intervals accepted by [as_ggchord_link()].
#' @param gene_data Optional directional gene input accepted by
#'   [as_ggchord_feature()]. Nondirectional intervals belong in `features`.
#' @param features,links Named lists of feature or paired-interval inputs.
#' @param mapping Named list with optional `seq`, `ribbon`, `gene` column
#'   mappings, and `features`/`links` lists of mappings named like the respective
#'   collections. Each column mapping is a named character vector.
#' @param steps Named list of explicit helper argument lists. Supported names:
#'   `clean`, `filter`, `deduplicate`, `merge`, `bundle`, `optimize`.
#'   Data arguments are supplied by preparation and cannot be overridden.
#' @return A list with `seq_data`, `ribbon_data`, `gene_data`, `features`,
#'   `links` and `report`. `report$valid` describes final validity, with
#'   `validation` before/after, `feature_issues`, per-link validations,
#'   `conversions`, `steps`, effective `step_arguments`, a ribbon `source_rows` list (original input rows
#'   contributing to each output row), `pair_counts`, and optional
#'   `layout_suggestion`. Structural conversion errors stop; biological
#'   reference/bounds errors remain in the report and are never repaired by
#'   default. Helper stage reports number rows in that stage's input.
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' seqs <- data.frame(accver = c("A", "B"), length = c(100, 100),
#'   genome = c("genome1", "genome2"))
#' sites <- data.frame(accver = "A", start = 10, end = 20,
#'   feature_id = "site1")
#' rb <- data.frame(qaccver = "A", saccver = "B", qstart = 1,
#'   qend = 20, sstart = 50, send = 31)
#' prepared <- prepare_ggchord_data(seqs, rb, features = list(sites = sites),
#'   steps = list(deduplicate = list(keep = "first")))
#' prepared$report$valid
#' ggchord(prepared$seq_data, prepared$ribbon_data) + geom_seq() +
#'   geom_link_ribbon() + geom_feature(data = prepared$features$sites) +
#'   coord_collinear()
prepare_ggchord_data <- function(seq_data, ribbon_data = NULL, gene_data = NULL,
    features = list(), links = list(), mapping = list(), steps = list()) {
  caller <- "prepare_ggchord_data()"
  mapping <- ggchord_prepare_options(mapping, c("seq", "ribbon", "gene", "features", "links"), caller)
  stage_order <- c("clean", "filter", "deduplicate", "merge", "bundle", "optimize")
  steps <- ggchord_prepare_options(steps, stage_order, caller)
  helpers <- list(clean = clean_ggchord_data, filter = filter_ggchord_ribbons,
    deduplicate = deduplicate_ggchord_ribbons, merge = merge_ggchord_ribbons,
    bundle = bundle_ggchord_ribbons, optimize = optimize_ggchord_layout)
  for (stage in names(steps)) {
    allowed <- setdiff(names(formals(helpers[[stage]])), c("seq_data", "ribbon_data", "gene_data"))
    ggchord_prepare_options(steps[[stage]], allowed, paste0(caller, " steps$", stage))
    if (stage != "clean" && is.null(ribbon_data)) ggchord_stop(
      caller, ": step '", stage, "' requires ribbon_data")
  }
  if ("clean" %in% names(steps)) {
    if (is.null(steps$clean$reversed_interval)) steps$clean$reversed_interval <- "keep"
    steps$clean$reversed_interval <- match.arg(steps$clean$reversed_interval,
      c("sort", "drop", "error", "keep"))
    if (identical(steps$clean$reversed_interval, "sort") &&
        any(c("deduplicate", "merge", "bundle", "optimize") %in% names(steps))) {
      ggchord_stop(caller, ": cleaning with reversed_interval = 'sort' loses endpoint direction ",
        "needed by later ribbon steps; use reversed_interval = 'keep'")
    }
  }
  seq_data <- as_ggchord_seq(seq_data, mapping = mapping$seq)
  if ("genome" %in% names(seq_data) &&
      (anyNA(seq_data$genome) || any(!nzchar(as.character(seq_data$genome))))) {
    ggchord_stop(caller, ": supplied genome membership must be non-missing and non-empty")
  }
  if (!is.null(ribbon_data)) ribbon_data <- as_ggchord_link(ribbon_data, mapping$ribbon)
  if (!is.null(gene_data)) gene_data <- as_ggchord_feature(gene_data, mapping$gene)
  features <- ggchord_prepare_collection(features, mapping$features, as_ggchord_feature, caller)
  links <- ggchord_prepare_collection(links, mapping$links, as_ggchord_link, caller)
  # Cleaning can overwrite a stored direction, so check its input as well.
  if ("clean" %in% names(steps) &&
      any(c("deduplicate", "merge", "bundle", "optimize") %in% names(steps))) {
    ggchord_prepare_direction(ribbon_data, caller)
  }
  conversions <- list(seq = attr(seq_data, "ggchord_conversion"),
    ribbon = attr(ribbon_data, "ggchord_conversion"), gene = attr(gene_data, "ggchord_conversion"),
    features = lapply(features, attr, "ggchord_conversion"),
    links = lapply(links, attr, "ggchord_conversion"))
  validate <- function() validate_ggchord_data(seq_data, ribbon_data, gene_data,
    check_duplicates = FALSE)
  before <- validate()
  source_rows <- if (is.null(ribbon_data)) list() else as.list(seq_len(nrow(ribbon_data)))
  stage_reports <- list()
  suggestion <- NULL
  for (stage in stage_order[stage_order %in% names(steps)]) {
    args <- if (stage == "clean") list(seq_data = seq_data, ribbon_data = ribbon_data, gene_data = gene_data) else
      if (stage %in% c("bundle", "optimize")) list(seq_data = seq_data, ribbon_data = ribbon_data) else
        list(ribbon_data = ribbon_data)
    # Optional transformations must not consume malformed numeric or reference
    # inputs. Only clean has explicit policies for such inputs.
    if (stage != "clean" && !validate_ggchord_data(seq_data, ribbon_data,
        check_duplicates = FALSE)$valid) ggchord_stop(
      caller, ": invalid ribbon input before '", stage, "'; inspect validation or request explicit cleaning")
    if (stage %in% c("deduplicate", "merge", "bundle", "optimize"))
      ggchord_prepare_direction(ribbon_data, paste0(caller, " before ", stage))
    result <- do.call(helpers[[stage]], c(args, steps[[stage]]))
    if (stage == "optimize") {
      suggestion <- result
      stage_reports[[stage]] <- result$report
      next
    }
    if (stage == "clean") {
      seq_data <- result$seq_data; gene_data <- result$gene_data
      ribbon_data <- result$ribbon_data
      if (!is.null(ribbon_data)) source_rows <- source_rows[ribbon_data$.source_row]
    } else {
      old_source <- source_rows
      if (stage %in% c("filter", "deduplicate")) {
        source_rows <- old_source[attr(result$data, "source_rows")]
      } else {
        column <- if (stage == "merge") "from_rows" else "source_rows"
        source_rows <- lapply(result$report[[column]], function(rows) {
          index <- as.integer(strsplit(rows, ",", fixed = TRUE)[[1]])
          unique(unlist(old_source[index], use.names = FALSE))
        })
      }
      ribbon_data <- result$data
    }
    stage_reports[[stage]] <- result$report
  }
  after <- if (!length(setdiff(names(steps), "optimize"))) before else validate()
  feature_issues <- ggchord_prepare_feature_issues(features, seq_data)
  link_validation <- lapply(links, function(x) validate_ggchord_data(seq_data, x,
    check_duplicates = FALSE))
  pair_counts <- data.frame(qaccver = character(), saccver = character(), n = integer())
  if (!is.null(ribbon_data) && nrow(ribbon_data)) {
    pair_counts <- stats::aggregate(list(n = rep(1L, nrow(ribbon_data))),
      ribbon_data[c("qaccver", "saccver")], sum)
  }
  list(seq_data = seq_data, ribbon_data = ribbon_data, gene_data = gene_data,
    features = features, links = links, report = list(
      valid = after$valid && !any(feature_issues$severity == "error") &&
        all(vapply(link_validation, function(x) x$valid, logical(1))),
      validation = list(before = before, after = after), feature_issues = feature_issues,
      link_validation = link_validation, conversions = conversions, steps = stage_reports,
      step_arguments = steps,
      source_rows = source_rows, pair_counts = pair_counts, layout_suggestion = suggestion))
}
