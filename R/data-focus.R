ggchord_explicit_focus_loci <- function(sequences, loci, expand, caller) {
  if (!is.data.frame(loci)) ggchord_stop(caller, ": loci must be a data frame")
  ggchord_require_columns(loci, c("accver", "start", "end", "locus_id"), caller)
  ggchord_conversion_coordinates(loci, c("start", "end"), caller)
  if (any(loci$start > loci$end)) ggchord_stop(caller, ": locus endpoints must be ordered")
  out <- ggchord_focus_validate_loci(sequences, loci, expand, caller)
  if (any(out$.focus_start != floor(out$.focus_start) | out$.focus_end != floor(out$.focus_end)))
    ggchord_stop(caller, ": expanded locus boundaries must be integral source positions")
  out
}

ggchord_focus_link_table <- function(data, sequences, loci, boundary, caller) {
  point <- all(c("qpos", "spos") %in% names(data)) &&
    !any(ggchord_ribbon_required_columns()[3:6] %in% names(data))
  canonical <- as_ggchord_link(data, kind = if (point) "point" else "interval")
  validate_ggchord_data(canonical, role = "link", sequences = sequences, action = "error")
  if (point) {
    rows <- list()
    original <- ggchord_focus_source(sequences)
    for (i in seq_len(nrow(canonical))) {
      q <- loci$locus_id[loci$accver == canonical$qaccver[i] &
        loci$start <= canonical$qpos[i] & loci$end >= canonical$qpos[i]]
      s <- loci$locus_id[loci$accver == canonical$saccver[i] &
        loci$start <= canonical$spos[i] & loci$end >= canonical$spos[i]]
      for (qid in q) for (sid in s) {
        row <- canonical[i, , drop = FALSE]
        for (side in c("q", "s")) {
          out_id <- if (side == "q") qid else sid
          id <- paste0(side, "accver")
          pos <- paste0(side, "pos")
          at <- match(row[[id]], sequences$accver)
          row[[paste0("source_", id)]] <- row[[paste0("source_", id)]] %||% original$accver[at]
          row[[paste0("source_", pos)]] <- row[[paste0("source_", pos)]] %||%
            (row[[pos]] + original$origin[at] - 1)
          row[[id]] <- out_id
          row[[pos]] <- row[[pos]] - loci$start[match(out_id, loci$locus_id)] + 1
        }
        rows[[length(rows) + 1L]] <- row
      }
    }
    return(ggchord_rbind_fill(rows) %||% canonical[0, , drop = FALSE])
  }
  # The clipping routine traverses intervals in alignment order. Recreate
  # that order transiently from evidence, then restore ordered coordinates.
  direction <- ggchord_link_orientation(canonical, caller)
  q_negative <- canonical$qstrand %in% "-"
  s_negative <- canonical$sstrand %in% "-"
  q_negative <- q_negative | (!canonical$qstrand %in% c("+", "-") &
    canonical$sstrand %in% c("+", "-") & ((direction == "opposite") != s_negative))
  s_negative <- s_negative | (!canonical$sstrand %in% c("+", "-") &
    ((direction == "opposite") != q_negative))
  for (side in c("q", "s")) {
    start <- paste0(side, "start")
    end <- paste0(side, "end")
    reverse <- if (side == "q") q_negative else s_negative
    saved <- canonical[[start]][reverse]
    canonical[[start]][reverse] <- canonical[[end]][reverse]
    canonical[[end]][reverse] <- saved
  }
  selected <- ggchord_focus_tables(sequences, ribbon_data = canonical, loci = loci,
    boundary = boundary)$ribbon_data
  clipped <- selected$.alignment_clipped %in% TRUE
  for (column in intersect(c("pident", "length"), names(selected))) {
    source <- if (column == "length") ".source_alignment_length" else ".source_pident"
    selected[[paste0("source_", column)]] <- selected[[source]]
    selected[[column]][clipped] <- NA
  }
  selected$coordinate_status <- ifelse(clipped, "estimated", "source")
  for (side in c("q", "s")) for (field in c("accver", "start", "end")) {
    prior <- paste0(".source_", side, field)
    if (!is.null(selected[[prior]])) selected[[paste0("source_", side, field)]] <- selected[[prior]]
  }
  ggchord_order_link_endpoints(selected, caller)
}

ggchord_focus_alignment_table <- function(data, sequences, loci) {
  rows <- list()
  root <- ggchord_focus_source(sequences)
  if (!"source_scope" %in% names(data) || !"source_accver" %in% names(data))
    return(data[0, , drop = FALSE])
  for (j in seq_len(nrow(loci))) {
    at <- match(loci$accver[j], sequences$accver)
    lower <- loci$start[j] + root$origin[at] - 1
    upper <- loci$end[j] + root$origin[at] - 1
    candidates <- which(data$source_scope %in% "source" & data$source_accver == root$accver[at])
    context <- intersect(c("source_input", ".source_input", "source_block", "accver"), names(data))
    groups <- split(candidates, interaction(data[candidates, context, drop = FALSE], drop = TRUE))
    for (rows_for_source in groups) {
      cells <- data[rows_for_source, , drop = FALSE]
      members <- which(!cells$is_gap & cells$source_position >= lower & cells$source_position <= upper)
      if (!length(members)) next
      keep <- seq(min(members), max(members))
      keep <- keep[cells$is_gap[keep] | keep %in% members]
      cells <- cells[keep, , drop = FALSE]
      cells$alignment_row <- cells$accver
      cells$accver <- loci$locus_id[j]
      cells$local_position <- cells$source_position - lower + 1
      rows[[length(rows) + 1L]] <- cells
    }
  }
  ggchord_rbind_fill(rows) %||% data[0, , drop = FALSE]
}

#' Crop explicit source windows together with their observation tables
#' @param x A data collection with a complete reference.
#' @param loci Window table with accver, start, end and unique locus_id.
#' @param expand Optional padding in source units.
#' @param boundary Trim overlaps or retain only complete observations.
#' @return A collection of the same class as x in local coordinates.
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
#' bundle <- ggchord_data(seqs, features=list(genes=genes), links=list(matches=links))
#' loci <- data.frame(accver=c("A","B"), start=c(50,100), end=c(700,700),
#'                    locus_id=c("A_window","B_window"))
#' local <- focus_ggchord_data(bundle, loci)
#' ggchord(local$sequences) + geom_gene(data=local$features$genes)
focus_ggchord_data <- function(x, loci, expand = 0, boundary = "trim") {
  caller <- "focus_ggchord_data()"
  if (!ggchord_is_data_collection(x)) ggchord_stop(caller, ": expected a data or import collection")
  boundary <- match.arg(boundary, c("trim", "drop"))
  sequences <- as_ggchord_seq(x$sequences)
  windows <- ggchord_explicit_focus_loci(sequences, loci, expand, caller)
  normalized <- loci
  normalized$start <- windows$.focus_start
  normalized$end <- windows$.focus_end
  out <- x
  base <- ggchord_focus_tables(sequences, loci = normalized)
  out$sequences <- base$seq_data
  out$sequences$source_accver <- out$sequences$.source_accver
  out$sequences$source_origin <- out$sequences$.source_origin
  out$sequences$source_length <- out$sequences$.source_length
  cropped <- normalized$start > 1 | normalized$end < sequences$length[match(normalized$accver, sequences$accver)]
  out$sequences$topology[cropped] <- "linear"
  out$sequences$circular[cropped] <- FALSE
  lineage <- list()
  for (name in names(out$features)) {
    features <- as_ggchord_feature(out$features[[name]])
    validate_ggchord_data(features, role = "feature", sequences = sequences, action = "error")
    scoped <- ggchord_focus_tables(sequences, feature_data = features,
      loci = normalized, boundary = boundary)
    table <- scoped$feature_data
    table$source_accver <- table$.source_accver
    table$source_start <- table$.source_start
    table$source_end <- table$.source_end
    out$features[[name]] <- table
    lineage[[paste("features", name, sep = "/")]] <- scoped$report
  }
  for (name in names(out$links)) {
    table <- ggchord_focus_link_table(out$links[[name]], sequences, normalized, boundary, caller)
    lineage[[paste("links", name, sep = "/")]] <- list(source_rows = table$.source_row,
      input_rows = nrow(out$links[[name]]), output_rows = nrow(table), estimated = table$coordinate_status)
    out$links[[name]] <- table
  }
  for (name in names(out$tracks)) {
    table <- out$tracks[[name]]
    type <- if (any(c("start", "end") %in% names(table))) "interval" else "point"
    out$tracks[[name]] <- focus_ggchord_track(table, sequences, normalized, type, boundary)
    lineage[[paste("tracks", name, sep = "/")]] <- data_report(out$tracks[[name]])$lineage
  }
  for (name in names(out$alignments)) {
    table <- ggchord_focus_alignment_table(out$alignments[[name]], sequences, normalized)
    lineage[[paste("alignments", name, sep = "/")]] <- list(source_rows = table$.source_row,
      input_rows = nrow(out$alignments[[name]]), output_rows = nrow(table))
    out$alignments[[name]] <- table
  }
  ggchord_attach_data_report(out, "focus_ggchord_data", inputs = data_report(x),
    lineage = lineage, metadata = list(windows = windows, boundary = boundary))
}

#' Crop data and arrange the resulting local references in rows
#' @param x A data collection with a complete reference.
#' @param loci Window table, optionally containing row_id.
#' @param boundary Trim overlaps or retain only complete observations.
#' @param orientation Display orientation for the local references.
#' @param segment_gap Physical gap between segments in the same row.
#' @return A list with local data, a sequence layout, and the crop report.
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
#' loci <- data.frame(accver=c("A","A"), start=c(1,400), end=c(300,700),
#'                    locus_id=c("left","right"), row_id="row1")
#' seg <- segment_ggchord_data(bundle, loci)
#' ggchord(seg$data$sequences) + geom_gene(data=seg$data$features$genes) +
#'   seg$layout + coord_collinear()
segment_ggchord_data <- function(x, loci, boundary = "trim",
                                  orientation = 1, segment_gap = 2) {
  data <- focus_ggchord_data(x, loci, boundary = boundary)
  sequences <- data$sequences
  row_id <- loci$row_id %||% loci$locus_id
  windows <- data.frame(accver = sequences$accver, instance_id = sequences$accver,
    start = 1, end = sequences$length, row_id = row_id)
  layout <- layout_sequences(order = sequences$accver, orientation = orientation,
    segment_gap = segment_gap, windows = windows)
  list(data = data, layout = layout, report = list(focus = data_report(data), windows = windows))
}
