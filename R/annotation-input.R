# Constructor annotations share one source model; rendering keeps its own rows.
ggchord_feature_source_palette <- function(data) {
  if (!is.data.frame(data) || !all(c("anno", "feature_color") %in% names(data))) return(NULL)
  usable <- !is.na(data$anno) & nzchar(as.character(data$anno)) &
    !is.na(data$feature_color) & nzchar(as.character(data$feature_color))
  colours <- stats::setNames(as.character(data$feature_color[usable]), as.character(data$anno[usable]))
  colours[!duplicated(names(colours), fromLast = TRUE)]
}

ggchord_normalize_feature_input <- function(data) {
  if (is.null(data)) return(NULL)
  if (inherits(data, "GRanges")) data <- as_ggchord_feature(data)
  data <- ggchord_normalize_accver(data)
  if (!is.data.frame(data)) return(data) # structured validation reports this
  if (!"strand" %in% names(data)) {
    data$strand <- rep(".", nrow(data))
    direction_column <- intersect(c("directionality", "direction"), names(data))[1L]
    if (!is.na(direction_column)) {
      map <- c(forward = "+", reverse = "-", bidirectional = "+/-",
        none = ".", nondirectional = ".")
      data$strand <- unname(map[tolower(as.character(data[[direction_column]]))])
    }
  }
  data
}

ggchord_annotation_fallback <- function(data) {
  unified <- !is.null(data$feature_data)
  input <- if (unified) data$feature_data else data$gene_data
  if (is.null(input)) return(NULL)
  attr(input, "ggchord_annotation_input") <- unified
  attr(input, "ggchord_source") <- list(
    table = if (unified) "feature" else "gene", rows = seq_len(nrow(input)))
  input
}

ggchord_annotation_gene_rows <- function(data) {
  if (!nrow(data)) return(integer())
  field <- intersect(c("feature_class", "type", "role"), names(data))[1L]
  if (is.na(field)) ggchord_stop(
    "geom_gene(): feature_data needs feature_class, type or role to select gene/CDS; ",
    "supply an explicit gene layer data table, or use geom_feature()")
  key <- tolower(trimws(as.character(data[[field]])))
  which(!is.na(key) & key %in% c("gene", "cds"))
}

ggchord_source_registry <- function(data) {
  parts <- lapply(names(data), function(table) {
    input <- data[[table]]
    if (!is.data.frame(input) || !nrow(input)) return(NULL)
    n <- nrow(input)
    field <- intersect(c("feature_class", "type", "role"), names(input))[1L]
    data.frame(source_id = paste0(table, ":", seq_len(n)), source_table = table,
      source_row = seq_len(n),
      accver = if ("accver" %in% names(input)) as.character(input$accver) else rep(NA_character_, n),
      feature_class = if (is.na(field)) rep(NA_character_, n) else as.character(input[[field]]),
      stringsAsFactors = FALSE)
  })
  ggchord_rbind_fill(parts) %||% data.frame()
}

# Reuse interval checks while retaining the generic feature strand vocabulary.
validate_feature_data <- function(data, seq_data, col, check_coordinates, check_duplicates) {
  if (is.null(data)) return(col)
  local <- validate_gene_data(data, seq_data, new_validation_collector(),
    check_coordinates, check_duplicates)
  local$errors <- Filter(function(x) !identical(x$category[1L], "invalid_strand"), local$errors)
  local$warnings <- Filter(function(x) !x$category[1L] %in% c("missing_anno", "overlapping"), local$warnings)
  if (is.data.frame(data) && "strand" %in% names(data)) {
    invalid <- which(is.na(data$strand) | !data$strand %in% c("+", "-", ".", "+/-"))
    if (length(invalid)) local <- add_validation_issue(local, "gene", "invalid_strand",
      invalid, "strand", "feature_data$strand must be '+', '-', '.' or '+/-'", "error")
  }
  for (kind in c("errors", "warnings")) {
    local[[kind]] <- lapply(local[[kind]], function(x) {
      x$table <- "feature"
      x$message <- gsub("gene_data", "feature_data", x$message, fixed = TRUE)
      x$message <- gsub("gene features", "features", x$message, fixed = TRUE)
      x$message <- gsub(" (clean_ggchord_data(out_of_range = 'clip') clamps it)",
        "; correct or explicitly clip the source intervals before plotting", x$message, fixed = TRUE)
      x
    })
    col[[kind]] <- c(col[[kind]], local[[kind]])
  }
  col
}

# Row-level source sidecars never overwrite user metadata columns.
ggchord_annotation_source_maps <- function(layout) {
  source_maps <- list()
  for (id in names(layout$layer_inputs)) for (component in names(layout$layer_inputs[[id]])) {
    if (!component %in% c("gene_poly", "gene_text", "gene_text_repel", "gene_label_segment", "gene_label_repel")) next
    input <- layout$layer_inputs[[id]][[component]]
    if (!is.data.frame(input) || !nrow(input)) next
    source <- attr(input, "ggchord_source")
    source_table <- source$table %||% id
    source_rows <- source$rows %||% seq_len(nrow(input))
    if (length(source_rows) != nrow(input)) ggchord_stop("Annotation source rows do not match resolved input")
    field <- intersect(c("feature_class", "type", "role"), names(input))[1L]
    source_maps[[length(source_maps) + 1L]] <- data.frame(
      layer_id = id, component = component, input_row = seq_len(nrow(input)),
      source_id = paste0(source_table, ":", source_rows), source_table = source_table,
      source_row = source_rows,
      feature_class = if (is.na(field)) rep(NA_character_, nrow(input)) else as.character(input[[field]]),
      source_scope = if (is.null(source)) "layer" else "constructor", stringsAsFactors = FALSE)
  }
  ggchord_rbind_fill(source_maps) %||% data.frame()
}

ggchord_annotation_duplicate_sources <- function(sources) {
  if (!nrow(sources)) return(data.frame())
  rows <- sources[sources$component == "gene_poly" & sources$source_scope == "constructor", , drop = FALSE]
  parts <- lapply(split(rows, rows$source_id), function(x) {
    ids <- unique(x$layer_id)
    if (length(ids) < 2L) return(NULL)
    data.frame(source_id = x$source_id[1L], n_geometry_layers = length(ids),
      layer_ids = paste(ids, collapse = ","), reason = "same_constructor_annotation_in_multiple_geometry_layers")
  })
  ggchord_rbind_fill(parts) %||% data.frame()
}
