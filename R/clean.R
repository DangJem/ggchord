#' Apply explicit repairs to observation collections
#' @param x Data or import collection.
#' @param unknown_id Whether to report or drop unknown references.
#' @param out_of_range Whether to report, drop or clip out-of-range observations.
#' @param reversed_interval Whether to report or normalize reversed endpoints.
#' @param invalid_pident Whether to report, drop or clear invalid identity values.
#' @param empty_label Whether to keep, drop or replace empty label values.
#' @param replacement_label Text used with empty_label = "replace".
#' @return A collection of the input class, with a change report.
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
#' bundle <- ggchord_data(seqs, features=list(genes=genes))
#' cleaned <- clean_ggchord_data(bundle, out_of_range="clip")
#' data_report(cleaned)$changes
clean_ggchord_data <- function(x, unknown_id = "error", out_of_range = "error",
    reversed_interval = "error", invalid_pident = "error", empty_label = "keep",
    replacement_label = "unannotated") {
  caller <- "clean_ggchord_data()"
  if (!ggchord_is_data_collection(x)) ggchord_stop(caller, ": expected a data or import collection")
  unknown_id <- ggchord_match_arg(unknown_id, c("error", "drop"))
  out_of_range <- ggchord_match_arg(out_of_range, c("error", "drop", "clip"))
  reversed_interval <- ggchord_match_arg(reversed_interval, c("error", "normalize"))
  invalid_pident <- ggchord_match_arg(invalid_pident, c("error", "drop", "set_na"))
  empty_label <- ggchord_match_arg(empty_label, c("keep", "drop", "replace"))
  if (!is.character(replacement_label) || length(replacement_label) != 1L || is.na(replacement_label))
    ggchord_stop(caller, ": replacement_label must be one character value")
  references <- as_ggchord_seq(x$sequences)
  collections <- c(features = "feature", links = "link", tracks = "track")
  changes <- list()
  lineage <- list()
  out <- x
  for (collection in names(collections)) {
    for (name in names(x[[collection]])) {
      table <- x[[collection]][[name]]
      path <- paste(collection, name, sep = "/")
      role <- collections[[collection]]
      validate_ggchord_data(table, role = role, sequences = references, action = "error",
        checks = c("schema", "units"))
      live <- rep(TRUE, nrow(table))
      original <- table$.source_row %||% seq_len(nrow(table))
      record <- function(index, column, before, after, action, reason) {
        if (!length(index)) return(invisible(NULL))
        changes[[length(changes) + 1L]] <<- data.frame(collection = path,
          source_row = original[index], column = column, before = as.character(before),
          after = as.character(after), action = action, reason = reason)
        invisible(NULL)
      }
      fail <- function(index, issue) {
        ggchord_stop(caller, ": ", path, " row ", original[index[1L]], " has ", issue)
      }
      for (id in names(ggchord_data_endpoints(table, role))) {
        columns <- ggchord_data_endpoints(table, role)[[id]]
        at <- match(table[[id]], references$accver)
        bad <- which(live & is.na(at))
        if (length(bad)) {
          if (unknown_id == "error") fail(bad, "an unknown reference")
          record(bad, id, table[[id]][bad], NA_character_, "drop", "unknown_id")
          live[bad] <- FALSE
        }
        for (column in columns) {
          values <- table[[column]]
          if (!is.numeric(values)) ggchord_stop(caller, ": ", path, "$", column, " must be numeric")
          bad <- which(live & (!is.finite(values) | values != floor(values)))
          if (length(bad)) fail(bad, "nonfinite or nonintegral coordinates")
        }
        if (length(columns) == 2L) {
          a <- table[[columns[1L]]]
          b <- table[[columns[2L]]]
          reversed <- which(live & a > b)
          if (length(reversed)) {
            if (reversed_interval == "error") fail(reversed, "reversed endpoints")
            circular <- references$topology[at[reversed]] == "circular"
            if (role != "link" && any(circular))
              fail(reversed[circular], "ambiguous circular endpoints; supply ordered segments")
            record(reversed, columns[1L], a[reversed], b[reversed], "normalize", "reversed_interval")
            record(reversed, columns[2L], b[reversed], a[reversed], "normalize", "reversed_interval")
            for (column in columns) {
              source <- paste0("source_", column)
              if (!source %in% names(table)) table[[source]] <- table[[column]]
            }
            table[[columns[1L]]][reversed] <- b[reversed]
            table[[columns[2L]]][reversed] <- a[reversed]
          }
        }
        start <- table[[columns[1L]]]
        end <- table[[tail(columns, 1L)]]
        limit <- references$length[at]
        outside <- which(live & (start < 1 | end > limit))
        if (!length(outside)) next
        if (out_of_range == "error") fail(outside, "out-of-range endpoints")
        no_overlap <- outside[start[outside] > limit[outside] | end[outside] < 1]
        dropped <- if (out_of_range == "drop" || length(columns) == 1L) outside else no_overlap
        record(dropped, id, table[[id]][dropped], NA_character_, "drop", "out_of_range")
        live[dropped] <- FALSE
        clipped <- setdiff(outside, dropped)
        for (column in columns) {
          before <- table[[column]][clipped]
          adjusted <- pmax(1, pmin(before, limit[clipped]))
          different <- before != adjusted
          record(clipped[different], column, before[different], adjusted[different], "clip", "out_of_range")
          if (any(different)) {
            source <- paste0("source_", column)
            if (!source %in% names(table)) table[[source]] <- table[[column]]
          }
          table[[column]][clipped] <- adjusted
        }
        if (length(clipped) && role == "link") {
          for (column in intersect(c("pident", "length"), names(table))) {
            source <- paste0("source_", column)
            if (!source %in% names(table)) table[[source]] <- table[[column]]
            table[[column]][clipped] <- NA
          }
          if (!"coordinate_status" %in% names(table)) table$coordinate_status <- rep("source", nrow(table))
          table$coordinate_status[clipped] <- "estimated"
        }
        if (length(clipped) && role == "feature" && "segments" %in% names(table)) {
          for (index in clipped) {
            parts <- table$segments[[index]]
            if (!is.data.frame(parts)) next
            parts$start <- pmax(parts$start, 1)
            parts$end <- pmin(parts$end, limit[index])
            parts <- parts[parts$start <= parts$end, , drop = FALSE]
            if (!nrow(parts)) {
              record(index, "segments", "outside", NA_character_, "drop", "out_of_range")
              live[index] <- FALSE
            } else {
              table$segments[[index]] <- parts
            }
          }
        }
      }
      if ("pident" %in% names(table)) {
        value <- table$pident
        if (!is.numeric(value)) ggchord_stop(caller, ": ", path, "$pident must be numeric")
        bad <- which(live & !is.na(value) & (!is.finite(value) | value < 0 | value > 100))
        if (length(bad)) {
          if (invalid_pident == "error") fail(bad, "invalid identity percentages")
          record(bad, "pident", value[bad], NA_real_, invalid_pident, "invalid_pident")
          if (invalid_pident == "drop") live[bad] <- FALSE else {
            if (!"source_pident" %in% names(table)) table$source_pident <- value
            table$pident[bad] <- NA_real_
          }
        }
      }
      if ("label" %in% names(table) && empty_label != "keep") {
        label <- as.character(table$label)
        bad <- which(live & (is.na(label) | !nzchar(label)))
        if (length(bad)) {
          replacement <- if (empty_label == "replace") replacement_label else NA_character_
          record(bad, "label", label[bad], replacement, empty_label, "empty_label")
          if (empty_label == "drop") live[bad] <- FALSE else {
            label[bad] <- replacement
            table$label <- label
          }
        }
      }
      out[[collection]][[name]] <- table[live, , drop = FALSE]
      lineage[[path]] <- original[live]
    }
  }
  validate_ggchord_data(out, action = "error")
  empty <- data.frame(collection = character(), source_row = integer(), column = character(),
    before = character(), after = character(), action = character(), reason = character())
  ggchord_attach_data_report(out, "clean_ggchord_data", inputs = data_report(x),
    lineage = lineage, changes = ggchord_rbind_fill(changes) %||% empty)
}
