# The collection organizes tables without assigning any plot-layer defaults.
ggchord_table_collection <- function(x, name, caller) {
  if (is.data.frame(x)) return(list(main = x))
  if (!is.list(x) || (length(x) && (is.null(names(x)) ||
      anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))) ||
      !all(vapply(x, is.data.frame, logical(1)))) {
    ggchord_stop(caller, ": ", name, " must be a data frame or a list of uniquely named data frames")
  }
  x
}

ggchord_import_collection <- function(sequences, features = list(), links = list(),
                                       tracks = list(), alignments = list(), report) {
  collections <- list(features = features, links = links, tracks = tracks, alignments = alignments)
  report <- c(list(operation = "import"), report)
  structure(c(list(sequences = sequences), collections, list(report = report)),
    class = "ggchord_import")
}

#' Organize a set of biological data tables
#'
#' A data collection keeps one complete reference dictionary and independently
#' named observation tables. Choose each table explicitly when drawing.
#'
#' @param sequences Complete reference sequences.
#' @param features Annotation data frame or named list of annotation tables.
#' @param links Point or interval relationship table or named list of such tables.
#' @param tracks Point or interval track table or named list of such tables.
#' @param alignments Character alignment table or named list of alignment tables.
#' @return A ggchord_data collection with a reference and four named table sets.
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
#' d <- ggchord_data(seqs,
#'   features = list(genes = genes), links = list(blast = links))
#' ggchord(d$sequences) + geom_seq() +
#'   geom_gene(data = d$features$genes) +
#'   geom_link_ribbon(data = d$links$blast)
ggchord_data <- function(sequences, features = list(), links = list(),
                          tracks = list(), alignments = list()) {
  ggchord_local_console_errors()
  caller <- "ggchord_data()"
  sequences <- as_ggchord_seq(sequences)
  collections <- list(features = features, links = links, tracks = tracks, alignments = alignments)
  for (role in names(collections)) {
    collections[[role]] <- ggchord_table_collection(collections[[role]], role, caller)
  }
  x <- structure(c(list(sequences = sequences), collections), class = "ggchord_data")
  # Keep coordinates available for explicit cleaning after collection.
  validate_ggchord_data(x, action = "error", checks = c("schema", "units"))
  ggchord_attach_data_report(x, "ggchord_data",
    inputs = list(sequences = data_report(sequences), tables = lapply(collections, names)),
    lineage = lapply(collections, function(tables) lapply(tables, function(table) {
      data_report(table, sections = "lineage")$lineage
    })))
}

ggchord_is_data_collection <- function(x) {
  inherits(x, "ggchord_data") || inherits(x, "ggchord_import")
}

ggchord_attach_data_report <- function(x, operation, inputs = NULL, lineage = NULL,
                                       issues = NULL, changes = NULL, metadata = NULL) {
  attr(x, "ggchord_data_report") <- list(operation = operation, inputs = inputs,
    lineage = lineage, issues = issues, changes = changes, metadata = metadata)
  x
}

#' Read the provenance and changes of a data operation
#'
#' @param x A ggchord data table or collection.
#' @param sections Named report sections to select. NULL returns all sections.
#' @return A report with section-level availability, operation, inputs,
#'   lineage, issues, changes and metadata.
#' @export
#' @examples
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' selected <- filter_ggchord_links(links, min_pident=95)
#' data_report(selected)
#' data_report(selected, sections=c("changes", "lineage"))
data_report <- function(x, sections = NULL) {
  ggchord_local_console_errors()
  report_names <- c("operation", "inputs", "lineage", "issues", "changes", "metadata")
  result <- attr(x, "ggchord_data_report")
  if (is.null(result)) {
    provenance <- attributes(x)[intersect(names(attributes(x)),
      c("ggchord_import", "ggchord_conversion", "ggchord_projection", "ggchord_focus",
        "ggchord_link_report", "ggchord_translation"))]
    report <- if (is.list(x) && !is.data.frame(x)) x$report else NULL
    if (is.null(report) && length(provenance)) {
      latest <- provenance[[length(provenance)]]
      lineage <- if (is.data.frame(x)) {
        fields <- intersect(c(".source_file", ".source_input", "source_input", ".source_row", ".source_line",
          ".source_record", ".source_feature", "source_rows", "source_row", "source_block",
          "source_segment", "source_members"), names(x))
        if (length(fields)) x[fields]
      }
      report <- list(operation = names(provenance)[length(provenance)], lineage = lineage,
        inputs = latest$files %||% latest$source_class, metadata = provenance)
    }
    if (!is.null(report)) {
      metadata <- report[setdiff(names(report), report_names)]
      result <- report[intersect(names(report), report_names)]
      if (length(metadata)) result$metadata <- c(result$metadata, metadata)
    }
  }
  requested <- sections %||% report_names
  if (!is.character(requested) || anyNA(requested) || anyDuplicated(requested) ||
      any(!requested %in% report_names)) {
    ggchord_stop("data_report(): sections must name distinct report sections")
  }
  # Source row indices and segments follow the current annotation rows.
  location <- attr(x, "ggchord_annotation_location", exact = TRUE)
  if (is.data.frame(x) && !is.null(location) &&
      all(c("accver", ggchord_annotation_coordinates(location)) %in% names(x)) &&
      any(c("lineage", "metadata") %in% requested)) {
    identity <- ggchord_annotation_identity(x, location, caller = "data_report()")
    if ("lineage" %in% requested) result$lineage <- ggchord_annotation_lineage(x, identity)
    if (location == "interval" && "metadata" %in% requested)
      result$metadata$segments <- ggchord_annotation_segment_index(x, identity)
  }
  available <- vapply(setNames(requested, requested), function(section) {
    if (is.null(result[[section]])) "unavailable" else "available"
  }, character(1))
  values <- setNames(lapply(requested, function(section) result[[section]]), requested)
  c(values, list(section_status = available))
}
