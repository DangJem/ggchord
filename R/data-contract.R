#' Data workflow and biological coordinate contracts
#'
#' ggchord separates reading, adaptation, validation, explicit cleaning and
#' preparation. None of these steps runs BLAST, infers a missing full sequence
#' length from features, verifies homology experimentally, or adds plot layers.
#'
#' \describe{
#'   \item{Read}{The six \code{read_*()} functions parse local files. Both
#'     \code{file} and \code{files} accept path vectors and glob patterns. When
#'     both are supplied, \code{file} precedes \code{files}; matches expand in
#'     pattern order, and each normalized path is read once, retaining its first
#'     supplied spelling. Missing, empty, unmatched or directory paths fail.
#'     Text input must be UTF-8 (an initial BOM is accepted by table readers).
#'     Convert other encodings explicitly before reading. Import attributes
#'     record format, files and output counts. Tabular readers and FASTA lengths
#'     retain physical \code{.source_line}; optional \code{.source_file} gives
#'     row-level file identity. GenBank instead retains record/feature indices
#'     and source text in its report. Malformed field counts stop, rather than
#'     shifting fields between rows. Numeric conversion problems in BLAST/GFF3
#'     warn and retain missing values for inspection and validation.}
#'   \item{Adapt}{\code{as_ggchord_*()} accepts documented tables or optional
#'     interval/sequence objects. Required IDs and positive integral coordinates
#'     are checked immediately; zero-based tables must be converted explicitly.
#'     Metadata name collisions receive prefixes and a column map. Required
#'     coordinate NA is an error, while unrelated metadata NA survives. Seqinfo,
#'     object metadata, names and upstream provenance are retained in the
#'     conversion attribute. Optional packages are never installed automatically.}
#'   \item{Validate}{\code{validate_ggchord_data()} reports all collected issues
#'     without changing input. \code{strict = TRUE} stops after collection if
#'     errors exist. Logical argument errors stop immediately. This validator
#'     describes constructor seq/ribbon/gene tables; generic unknown-strand
#'     features belong in preparation's named feature collections. Coordinate
#'     checks do not prove an alignment or a biological annotation.}
#'   \item{Clean}{\code{clean_ggchord_data()} explicitly applies selected
#'     policies and reports original rows and old/new values. Its historical
#'     \code{reversed_interval = "sort"} is retained; use \code{"keep"} for
#'     directional links or circular-origin features. Sorting cannot recover
#'     their original meaning later. Unknown-ID \code{"keep"} is an audit
#'     policy, not a guarantee that the result is drawable. NULL means an absent
#'     optional table, and an empty table remains a table.}
#'   \item{Prepare}{\code{prepare_ggchord_data()} adapts first, validates, then
#'     runs only requested steps in the documented order. Its cleaning default
#'     is \code{reversed_interval = "keep"}, because downstream operations use
#'     endpoint directions. Adaptation errors therefore occur before cleaning;
#'     clean/drop malformed IDs separately when needed. Named features/links
#'     remain independent collections. Reports distinguish stage-local rows from
#'     composed original ribbon source sets.}
#' }
#'
#' @section Coordinates and direction:
#' Genomic intervals are one-based inclusive: \code{start == end} is one base.
#' The explicit adapters enforce integral coordinates. Legacy constructor
#' validation, cleaning and raw focus accept finite numeric positions; their
#' success alone does not certify integral genomic coordinates. Prefer integral
#' windows for biological data and use adapters before these legacy entries
#' when strict integer checks are required.
#' BED/PAF half-open \code{[a,b)} becomes \code{[a+1,b]}; zero-width records
#' cannot become nonempty intervals and are rejected. Link endpoint order can
#' encode relative direction; it is not a circular-origin join. Feature
#' \code{start > end} is context-dependent: circular plotting can split an origin
#' join, while linear window workflows require explicit ordered segments and an
#' envelope. A raw focus locus uses min/max of its bounds, not an origin wrap.
#' \code{strand} is biological annotation, \code{orientation} is display choice,
#' and link \code{direction} is relative. Single-base link endpoints have no
#' order: supplied strand or
#' relative direction can provide evidence, otherwise direction remains unknown.
#' Merge and bundle retain these rows independently instead of manufacturing
#' aggregate endpoint directions. Unknown orientation is kept separate from
#' known direction in density and deduplication.
#' Unknown strands remain unknown in
#' adapters; GFF3's explicit legacy \code{unstranded = "plus"} is an exception,
#' with \code{source_strand} retained. Prefer \code{unstranded = "keep"} and
#' generic feature layers when direction is unknown.
#'
#' Plotting uses a continuous display axis, not a base-cell raster: feature/axis
#' chord anchors use pos/length, point-link anchors use (pos-1)/length. These
#' historical graphical offsets do not change stored coordinates. Single-base
#' features can use a point or widened glyph. No new endpoint convention is
#' imposed on plots in this data-workflow revision.
#'
#' @section Selection, projection and aggregation:
#' Focus crops all selected tables and composes earliest source coordinates.
#' Overlapping windows can copy a row. Ribbon clipping interpolates an endpoint
#' correspondence, rounds to bp and marks estimated length; identity is retained,
#' not recalculated. Segment preparation creates separate windows and explicit
#' gaps; it never concatenates their biological sequences. Position mapping is
#' invertible within the selected linear window and rejects segmented gaps.
#' CDS projection uses observed complete codons and explicit phase/partial
#' policies; it does not translate with the stored genetic-code ID. Protein
#' identity/length remain source-prefixed metadata when links are projected.
#'
#' Filtering retains full rows. Deduplication selects representatives and
#' reports discarded rows; it does not combine their statistics. Merge creates
#' envelopes and reports source sets, lengths and identities. Bundle groups by
#' pair, relative endpoint direction and midpoint bins; both endpoints may be
#' reversed together in one group. Envelope endpoints follow the first row's
#' absolute direction. Bundle length is a rounded length-weighted mean of source
#' lengths, identity is a length-weighted mean, and conflicting extra metadata
#' become typed NA. These summaries are not new alignments. Supply
#' \code{group_by} for metadata that must remain separate. Density retains every
#' input row and measures bin weight relative to that pair's maximum bin; it
#' is not coverage or statistical probability. Optimization suggests display
#' order/orientation without changing biological strand or query/subject roles.
#'
#' @section Search helpers:
#' Cluster links use declared membership, with explicit paralog and segment
#' policies, not sequence alignment. Primer exact and three-prime suffix searches
#' report zero mismatches; a partial suffix is not full-oligo validation. Wallace
#' temperature is a simple composition estimate, not an experimental prediction.
#' Restriction results are motif matches and source-backed cut offsets; missing
#' cuts remain unknown. Methylation status uses supplied/reference context and a
#' limited sensitivity catalogue, not a complete digestion model. Linear context
#' never wraps sequence ends. Common-feature confidence is a deterministic
#' heuristic rank, not a probability or experimental confirmation. Its six-frame
#' protein search uses the standard codon table without alternative starts,
#' splicing or circular-origin protein matches. Source \code{geneticCode}
#' metadata is not interpreted as a selectable translation table. Use exact DNA
#' mode or an external code-aware analysis for nonstandard translation.
#'
#' @name ggchord_data_contract
#' @seealso \code{\link{prepare_ggchord_data}}, \code{\link{read_genbank}},
#'   \code{\link{as_ggchord_feature}}, \code{\link{project_ggchord_links}}
NULL
