#' Example gene annotation data
#'
#' A compact plotting fixture generated deterministically from
#' \code{examples/gene_track.tsv}. Eight features of varied visible lengths are
#' distributed across each sequence and favour informative annotations.
#' Because the source records are
#' strongly strand-biased, the demonstration \code{strand} alternates between
#' \code{"+"} and \code{"-"}; the biological source value remains available as
#' \code{source_strand}.
#'
#' @format A data frame containing the following columns:
#' \itemize{
#'   \item accver: Sequence accession/version or an unchanged custom identifier
#'   \item start: Gene start position
#'   \item end: Gene end position
#'   \item strand: Strand direction (+/-)
#'   \item anno: Gene annotation category
#'   \item source_strand: Strand direction in the source annotation file
#' }
"gene_data_example"

#' Example sequence data
#'
#' Sequence length data for ggchord demonstration
#'
#' @format A data frame containing columns: accver, length
"seq_data_example"

#' Example alignment data
#'
#' A compact alignment fixture generated deterministically from the BLAST
#' files under \code{examples/blastn}. Alignments shorter than 300 bases are
#' omitted, and dense sequence pairs retain at most three spatially distributed
#' representatives so the default chord remains readable. Sparse pairs are
#' preserved unchanged.
#'
#' @format A data frame containing standard alignment columns (qaccver, saccver, pident, etc.)
"ribbon_data_example"

#' Minimal single-genome sequence fixture
#'
#' @format A one-row data frame with `accver` and `length`.
#' @source `examples/single-genome/minimal.fasta`.
#' @name single_genome_example
NULL

#' Minimal single-genome gene fixture
#'
#' @format A data frame with `accver`, `start`, `end`, `strand`, and `anno`.
#' @source `examples/single-genome/features.tsv`.
#' @name single_gene_example
NULL

#' Restriction sites for the minimal single-genome fixture
#'
#' @format A data frame returned by [find_restriction_sites()].
#' @source Generated from `examples/single-genome/minimal.fasta` by
#'   `data-raw/generate_example_data.R`.
#' @name restriction_site_example
NULL

#' Circular plasmid sequence examples
#'
#' Thirteen one-row data frames generated from every reference FASTA record
#' under `examples/plasmid/`. The file stem is used unchanged as both `accver`
#' and the display `label`; the DNA sequence is copied without biological
#' modification. Feature and restriction annotations are generated at run time
#' with [find_common_features()] and [find_restriction_sites()].
#'
#' @format A one-row data frame with `accver`, display `label`, `length`, and
#'   the DNA `sequence`.
#' @source `examples/plasmid/`.
#' @name plasmid_example_pUC19
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pUC19c
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pBR322
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pBluescript_II_SK_plus
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pCAMBIA1300
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pcDNA3_1_plus
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pDONR221
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pEarleyGate_201
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pET_28a_plus
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pETDuet_1
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pSB1C3
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pSpCas9_BB_2A_GFP_PX458
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pTRE_Tight_BI
NULL

#' @rdname plasmid_example_pUC19
#' @name plasmid_example_pTRIPZ
NULL
