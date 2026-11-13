{{code:courseblockaz}}

{{prepare:setup}}

Extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains examples/. Each section names its optional dependencies.

## Variant records and structural breakpoints {#variants}

Install optional **VariantAnnotation** for VCF and additionally **StructuralVariantAnnotation** for explicit structural decoding. These examples read the dependencies' unchanged source-backed examples: the chr22 1000 Genomes subset and a real caller VCF distributed with StructuralVariantAnnotation. Data stays in those packages; no third-party parser implementation or input is copied into ggchord. `check-variant-import.R` records exact package versions and source hashes.

{{code:courseblockba}}

Each placed record remains one feature row with an ALT list, including multiallelic calls. Features show the reference allele span, including the padding base, rather than an inferred symbolic SV footprint. Complete INFO, sample/genotype arrays and headers remain in `report$records[[1]]$parsed`. Inspect actual caller fields there; neither missing support nor Identity becomes zero or a fabricated count. Unknown full lengths remain NA; provide reliable reference lengths before plotting. Telomeric/out-of-bounds spans stay in `unplaced_features` and the source, without clipping.

Structural decoding is opt-in, keeps confidence intervals, disables absent-mate inference, and creates one link per reciprocal pair with both original record indices. Warnings identify records not paired by the decoder, while all records remain in the report. `.plottable = FALSE` marks breakends outside known reference bounds; affected pairs remain in `unplaced_links`. Exact single-base breakpoints suit `geom_link_line()` after explicitly selecting their positions; uncertainty intervals suit interval geometry. `q_breakend_orientation` and `s_breakend_orientation` describe the derivative sequence, not relative homology strand. Never manufacture alignment Identity to colour these links.

For modern BCF (`.bcf` or `.bcf.gz`), supply `python` naming an existing interpreter with optional **pysam**: `read_vcf("calls.bcf", python = "/path/to/environment/bin/python")`. Only BCF invokes it, transcodes to temporary VCF via HTSlib, and records decoder/version/messages. Virtual-environment interpreter paths retain their environment. Original record order, values and genotypes are preserved; transcoded text is not an original physical BCF line. No installation occurs. This reader loads selected local files into memory; select or subset large datasets upstream.

## Genomic alignment blocks: MAF and chain {#genome-alignments}

Use an existing Python environment with optional **Biopython Bio.Align** MAF/chain support (directly verified with 1.85). `ALIGNMENT_PYTHON` can select that interpreter for this lesson. The source checkout/input bundle must contain `examples/genome-alignments/`. Its two unchanged fixed-version upstream files have hashes, permission attribution and a verification script; no aligner is run. An unavailable optional parser skips the example, not a successful import.

{{code:courseblockbb}}

MAF retains the actual gapped characters and full declared source lengths. Negative components map to decreasing **forward source** coordinates. First-component anchor relationships, source blocks/files and original coordinate matrices remain explicit; a block is not automatically expanded to every possible pair. Gap cells have missing source positions, and full lengths alone do not independently verify symbols. Filter one file/block before `map_alignment_position()` or a character panel because columns restart per block. Original quality, context and empty-component evidence stays in the source/report. An empty first anchor never silently becomes another anchor.

The chain example explicitly selects source block 1, the source-backed NR_046654.1/chr3 relationship; upstream `_modified` format controls remain in the unchanged input and are not claimed as new biological observations. The three aligned exon blocks retain intervening gaps. Chain provides no mismatches or character cells, and a chain score is not Identity. Query maps to qaccver, target to saccver; original IDs remain separately traceable. Declare `sequence_unit` explicitly when its default bp is inappropriate.

For two assemblies both named `chr1`, declare identities, for example `query_ids = c(chr1 = "assemblyB.chr1")` and `subject_ids = c(chr1 = "assemblyA.chr1")`. Defaults preserve IDs and reject ambiguous cross-role collisions; they never auto-prefix or merge assemblies. Full-length conflicts fail, and shared final IDs require explicit consistent declaration on both sides. This is source identity, separate from display instances. No pident, support count or complete source string is inferred from either format's endpoints.

## Existing delta alignments {#delta}

`read_delta()` requires an existing **show-coords** program; `DELTA_DECODER` can select its executable for this lesson. It reads existing NUCMER/PROMER results, never invokes an aligner or opens the historical FASTA paths in the file preamble. The unchanged upstream Drosophila example, pinned version, original permission notice and hashes are in `examples/genome-alignments/`. An absent executable skips this optional example; that is not successful decoding.

{{code:courseblockbc}}

This explicitly selected query remains the literal string `"3210101"`. Query/reference map to qaccver/saccver without shifting the original one-based directional endpoints. Every native output row is checked against its actual source header/endpoints before error counts, signed delta scripts and physical lines are attached. Full lengths come from source declarations, not the last alignment. Cross-role source-name collisions follow the explicit ID maps of `read_chain()`; no assembly identity is guessed.

`pident` is the tool's reported percentage with its printed precision, not a value reconstructed from two endpoints. In PROMER it is an **AA alignment statistic**, while interval coordinates and full sequence lengths remain **bp**. Native similarity/stop percentages and frames remain separate. qspan/sspan are genomic span widths, not alignment-column counts; a BLAST `length`, read support and sequence characters are not manufactured. Original text and native output/messages remain in `report`. If the executable has no reliable version option, its version remains unknown in the reader report; record your prepared environment's provenance separately. The directed audit records the actual fixed tool source tag and executable checksum. No third-party executable or implementation is redistributed in ggchord.

The source bundle also includes `arabidopsis-first-pair.delta`, an unchanged prefix containing the original preamble and all 39 records of the first Arabidopsis sequence-pair block. `arabidopsis-source.json` records the pinned original compressed file, decompressed hash, complete-block rule and unchanged product hash; the author's MIT permission notice is retained. This checks real nucleotide alignments independently of the translated example, without newly running alignment or certifying the remaining source blocks. Use the same existing decoder:

{{code:courseblockbd}}

## Existing protein signatures {#protein-signatures}

These installed-package inputs contain existing Pfam annotations cached from the InterPro API for two human proteins. The original JSON bytes, retrieval date, URLs, CC0 attribution and hashes are in `inst/extdata/signature-sources.json`. `data-raw/prepare-signature-inputs.py` explicitly re-encodes the same observed fragments as TSV for a coordinate round trip. It does not run InterProScan: absent MD5, status, run date and integrated description remain `-`. The TSV is a derived format example, not a new analysis result.

{{code:courseblockbe}}

Full lengths are declared AA lengths, not maximum feature ends. The API returns lowercase IDs (`p00517`, `p51587`); these remain unchanged, as do accession versions and repeated model matches. Model/signature IDs are not unique feature-instance IDs, so `feature_id` stays missing. Any source identity reconciliation with another file must be explicit.

Each cached JSON location is one feature with actual `segments`; its start/end only bound the footprint. Existing Feature geometry leaves coordinate gaps empty and retains one source match. The representative flag is metadata, never an automatic filter. Original entry/protein/location indices and per-file `.source_match` stay traceable; use `.source_input` when combining files. TSV rows lack JSON grouping evidence and are not automatically grouped by model ID. Scores remain a native list-column: JSON numbers, TSV strings and missing values retain their source type. Do not compare different libraries' score/e-value systems as a universal confidence or invent identity/read support. For a numeric scale, first select the intended library and establish its actual score meaning explicitly.

JSON needs optional **jsonlite** and reads only supplied local caches. Remaining-page URLs and declared counts stay in the report; nothing is fetched at draw time. Quality/context metadata, original fields and complete source objects remain available for auditing. These location annotations provide no complete protein characters or independent character verification. The [official format guide](https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html) and [score explanation](https://interpro-documentation.readthedocs.io/en/latest/faq.html#why-are-there-no-e-values-associated-with-interpro-entries) describe the underlying conventions.

## Existing feature relationships {#feature-relations}

`read_feature_relations()` reads existing pairs or group memberships offline. The installed plant subsets preserve complete original MCScanX data lines from a pinned source commit. Permission, original URLs, hashes and physical subset indices are recorded in `inst/extdata/feature-relation-sources.json`; `data-raw/prepare-feature-relations.py` reproduces them without running an analysis.

{{code:courseblockbf}}

The caller declares ascending one-based-inclusive contiguous anchors; the file suffix does not establish the coordinate convention. The companion position tables have no full chromosome lengths, so this example inspects relationships and does not invent `seq_data` to draw a whole chromosome. Supply independently sourced full lengths before plotting `ribbon_data` as interval relationships. These are not basewise HSPs. Map a score only after establishing its meaning; the third input column remains an opaque string, not Identity or e-value.

`pairs` accepts two IDs and an optional score. `named_groups` retains the first field as a native character group ID, including numeric-looking IDs; `groups` treats every field as a member and only records a file-local group index. Memberships never expand into all pairwise connections. Duplicates, singleton groups, member order, raw fields and physical lines stay visible.

Optional `feature_data` matches exact unique `feature_id` values. Case, versions and prefixes are never edited; the pair example and group example have different original ID namespaces and companion tables. Missing IDs fail by default; `unmatched = "report"` preserves original rows and flags incomplete mappings. Ambiguous duplicate feature IDs always fail. Without positions the reader still returns all relationships, with no mapped links. Annotation strand is retained separately, without assigning alignment direction; segmented annotations require an explicit contiguous anchor table. No full lengths, Identity, support or source characters are inferred.

## Graph paths: exploratory subsets {#graph-paths}

GFA graph edges describe adjacency/overlap, rather than automatically providing HSPs between whole genomes. Existing explicit path steps can be shown through ordinary sequence/feature tables when their coordinate construction is justified. Node IDs can recur within one path; retain the path, source-step index, occurrence and traversal orientation. Do not pair every occurrence or equate node/edge counts with orthology, Identity or read support.

A source-backed exploratory example uses ODGI's unchanged human LPA graph at a pinned commit, decoded by optional **gfapy**. All 5,195 original edges have zero overlap; the prototype validates every original path transition and node length, then projects 202,806 steps across 13 supplied paths. An explicit display subset keeps the first 100 steps of each of two named paths, with original traversal orientation and IDs. Complete supplied path lengths remain separate from chromosome lengths. No pair links are generated. Code, source URL/hash, original MIT notice, output tables and actual image are in `tools/long-term-audit/OPTIONAL-GRAPH-BRIDGE-v0.27.0.md` and its linked prototype artifacts in the source checkout; `data-raw/prepare-graph-prototype.py` retrieves only the two pinned originals on explicit invocation.

This is not a generic `read_gfa()` or graph Coord. Nonzero/unknown overlaps, missing source characters, GFA W/J/GFA2 and native ODGI binaries are outside this prototype. Use mature upstream tools to export an explicit path subset and retain their log/source mapping; do not infer missing coordinate or pair evidence. A documented upstream export route is distinct from a binary decoder actually validated here.
