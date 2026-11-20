#' Biological data and source coordinates
#'
#' Sequence tables carry stable reference IDs, complete lengths, coordinate
#' units and biological topology. Reference lengths are supplied by the
#' source or declared explicitly. Features and interval links use ordered,
#' one-based closed coordinates. Biological strands and relative alignment
#' direction remain separate from endpoints and display orientation.
#'
#' [ggchord_data()] optionally organizes a reference and named observation
#' tables. Each drawing layer uses its own observation table. [data_report()]
#' retrieves source information and operation results, including the members
#' of aggregated relationships and the original positions of focused windows.
#'
#' [validate_ggchord_data()] reports input issues and [clean_ggchord_data()]
#' applies explicitly chosen repairs. [focus_ggchord_data()] builds new local
#' data coordinates. A layout display window leaves the reference dictionary
#' and source statistics in their original coordinate system.
#'
#' Character alignments use real alignment columns and nullable source
#' positions. A gap consumes a column without occupying a source residue.
#' Restriction cut boundaries are interbase coordinates, from zero to the
#' complete reference length. Source coordinate metadata records unit and
#' origin without relying on the selected plotting coordinate system.
#'
#' @name ggchord_data_contract
NULL

# Resolve the native coordinate domain without inferring it from alphabet or IDs.
ggchord_sequence_unit <- function(seq_data, caller = "ggchord") {
  if (!is.data.frame(seq_data) || !"sequence_unit" %in% names(seq_data)) return("bp")
  units <- as.character(seq_data$sequence_unit)
  if (!length(units) || anyNA(units) || any(!units %in% c("bp", "aa"))) {
    ggchord_stop(caller, ": sequence_unit must be 'bp' or 'aa'")
  }
  if (length(unique(units)) != 1L) ggchord_stop(caller,
    ": bp and aa cannot share one native coordinate domain; project protein intervals to CDS coordinates or use separate plots")
  units[1L]
}

# Apply explicit interval declarations to the sequence coordinate domain.
ggchord_check_interval_unit <- function(data, unit, caller) {
  if (is.data.frame(data) && nrow(data) && "sequence_unit" %in% names(data) &&
      !identical(ggchord_sequence_unit(data, caller), unit)) {
    ggchord_stop(caller, ": interval sequence_unit (", ggchord_sequence_unit(data, caller),
      ") disagrees with the sequence table (", unit, ")")
  }
  invisible(NULL)
}
