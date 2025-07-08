#' Example gene annotation data
#'
#' Gene annotation data for ggchord demonstration (short genes have been filtered out)
#'
#' @format A data frame containing the following columns:
#' \itemize{
#'   \item seq_id: Sequence ID
#'   \item start: Gene start position
#'   \item end: Gene end position
#'   \item strand: Strand direction (+/-)
#'   \item anno: Gene annotation category
#' }
"gene_data_example"

#' Example sequence data
#'
#' Sequence length data for ggchord demonstration
#'
#' @format A data frame containing columns: seq_id, length
"seq_data_example"

#' Example alignment data
#'
#' BLAST alignment data for ggchord demonstration (length >= 100)
#'
#' @format A data frame containing standard BLAST columns (qaccver, saccver, pident, etc.)
"ribbon_data_example"
