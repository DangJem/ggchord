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

#' Circular plasmid sequence fixtures
#'
#' Three unchanged FASTA records used by the v0.13 circular-map examples.
#'
#' @format A three-row data frame with `accver`, display `label`, `length`, and
#'   the DNA `sequence`.
#' @source `examples/plasmid/`.
"plasmid_sequence_example"

#' Compact plasmid feature fixtures
#'
#' Core pBR322 and pUC19c annotations used to demonstrate gene and general
#' feature layers. This compact table is not a complete replacement for the
#' source records' feature tables.
#'
#' @format A data frame with `accver`, `start`, `end`, `strand`, `type`,
#'   `anno`, and `source`.
#' @source GenBank J01749.1 and L09137.2.
"plasmid_feature_example"

#' Common restriction sites in the plasmid examples
#'
#' Matches for a small independently specified panel of common motifs. It is
#' generated from `plasmid_sequence_example`; the complete REBASE database is
#' not bundled pending confirmation of redistribution terms.
#'
#' @format A data frame returned by [find_restriction_sites()].
#' @source `examples/plasmid/` and `data-raw/generate_example_data.R`.
"plasmid_restriction_example"
