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
#' The `plasmid_example_pUC19c` name is an alias
#' of `plasmid_example_pUC19`, including its `accver`. It must not be confused
#' with an independently imported GenBank pUC19c record such as L09137.2.
#'
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


#' Four complete protein sequences for comparison
#'
#' Three EMALE proteins and one Mavirus protein selected from public complete
#' FASTA records. Their six pairwise comparisons are in
#' [protein_ribbon_data_example]. Lengths describe complete proteins.
#'
#' @format A data frame with four rows and columns:
#' \describe{
#'   \item{accver}{Unchanged protein identifier.}
#'   \item{length}{Full length in amino acids.}
#'   \item{sequence_unit}{The declared coordinate unit, `"aa"`.}
#'   \item{sequence}{Complete protein residues retained from the source FASTA.}
#' }
#' @source Public EMALE/Mavirus FASTA records. Selection, source hashes,
#'   BLASTp version and arguments are in `examples/protein/sources.json`.
#'   Rebuild with `data-raw/prepare-bundled-proteins.py` and
#'   `data-raw/generate_protein_data.R`.
#' @seealso [protein_ribbon_data_example], [ggchord]
"protein_seq_data_example"

#' Six pairwise BLASTp alignments between the example proteins
#'
#' Observed BLASTp alignments between every unordered pair of the four
#' complete proteins in [protein_seq_data_example]. The raw all-against-all
#' output retains self hits and both query/subject directions. This teaching
#' table retains one direction per pair, ordered by the sequence table.
#'
#' @format A data frame with six rows containing the twelve standard BLAST
#'   fields, physical `.source_line` provenance and `sequence_unit = "aa"`.
#'   Endpoints are one-based inclusive amino-acid positions.
#' @source `examples/protein/all-vs-all-blastp.tsv`, generated with NCBI
#'   BLASTp 2.17.0+ from unchanged complete protein records. The command,
#'   source hashes and sequence selection are recorded in
#'   `examples/protein/sources.json`.
#' @seealso [protein_seq_data_example], [geom_link_ribbon]
"protein_ribbon_data_example"

#' Complete genome comparison inputs and flagship display subsets
#'
#' Unchanged sequence, GFF3 and BLAST records for four genomes, together with
#' the explicit subsets used by the website flagships. Rebuild with
#' `data-raw/generate_website_data.R`; source checksums are recorded in
#' `data-raw/website-data-sources.tsv`. No intervals or strands are invented.
#' Full source tables remain separate from the display selections.
#'
#' @format A named list: `sequences` includes complete DNA strings; `seq_data`
#'   contains accession and length; `annotations`, `alignments` and
#'   `source_alignments` retain complete imported records (the latter adds
#'   source-file and HSP identifiers). `genes` contains all CDS. `links`,
#'   `features`, `chord_features`, `linear_features`, `adjacent_links` and
#'   `cross_links` are the existing deterministic display selections, in their
#'   original order. Coordinates are one-based inclusive, lengths are in bp,
#'   and pident is a percentage. Displayed negative subject endpoints may be
#'   descending. Import provenance is retained on the imported tables.
#' @source Repository `examples/fasta/`, `examples/gff3/` and
#'   `examples/blastn-standard/`; their source documentation accompanies the
#'   original files. The generator records the exact filtering rules.
#' @examples
#' data(genome_comparison_example)
#' head(genome_comparison_example$seq_data)
"genome_comparison_example"

#' Source-backed alignment and coding sequence teaching inputs
#'
#' Prepared inputs for sequence-detail, Dotplot and Fusion examples. File
#' readers are demonstrated separately; these tables require no checkout or
#' working-directory configuration. Rebuild with
#' `data-raw/generate_website_data.R`.
#'
#' @format A named list: `plasmid_source` retains the complete pUC19 input;
#'   `paf` is the complete 11-HSP EMALE pair; `dna_alignment`
#'   and `protein_alignment` are long-form observed alignment cells;
#'   `aligned_strings` contains the original gapped DNA strings;
#'   `coding_genome`, `coding_features` and `coding_proteins` retain the genomic,
#'   annotation and protein inputs for the CDS translation example.
#'   `protein_hits` and `source_coding_features` retain the full protein
#'   projection inputs with source-file provenance. DNA
#'   coordinates are in bp, protein coordinates in residues. No alignment gaps
#'   or coding relationships are inferred to improve a plot.
#' @source Repository `examples/alignment/` and `examples/format-inputs/`,
#'   with upstream provenance in their README files and checksums in
#'   `data-raw/website-data-sources.tsv`.
#' @examples
#' data(sequence_detail_example)
#' head(sequence_detail_example$paf)
"sequence_detail_example"
