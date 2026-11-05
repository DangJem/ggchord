These helpers inspect inputs, select windows and consume results without inferring biological relationships for the analyst. Retain source tables and choose tools by a declared question.

## Prepare this page {#input}

Independently load ggchord and optional dplyr, using bundled complete sequences, observed HSPs and real pUC19 GC windows.

{{code:comparescaledata}}
{{code:trackdata}}


## A catalogue is not a search result {#catalog}

`primer_catalog` returns oligonucleotide records; aliases=TRUE returns source aliases. Catalogue direction is descriptive; find_primer_bindings determines binding direction on target DNA. Multiple names for one sequence do not represent distinct oligonucleotides. Continue with [primers](primers.html) to search and draw.

{{code:comparecatalog}}
{{console:comparecatalog}}


## Target-size diagnostics and explicit ordering {#diagnostics}

`diagnose_ggchord` inspects one PNG size and returns measurement, text and guide hints without mutating the plot. `export_ggchord_layout` returns layer tables and original_data=TRUE retains source rows. `get_chord_layout` exposes the complete public layout. `optimize_ggchord_layout` returns proposed order/directions and before/after scores; it neither edits inputs nor applies changes automatically. Readability optimization is not biological correctness.

{{code:comparediagnostics}}
{{console:comparediagnostics}}


## Tracks and windows share source coordinates {#track-focus}

Focus pUC19 to 500–1500 bp and trim the same observed GC windows as intervals. `focus_ggchord_track` does not recompute averages, invent boundary samples or bridge omitted windows. Return to actual DNA and declare new windows/denominators if finer GC summaries are needed.

{{code:comparetrackfocus}}
{{console:comparetrackfocus}}


## Where to learn each helper family {#tools}

| Task | Input, helpers and output | Main lesson |
| --- | --- | --- |
| Read files | read_fasta_lengths/read_gff3/read_blast/read_bed/read_paf/read_genbank: source files to tables and provenance | [Formats](formats.html), [alignments](alignments.html) |
| Adapt supported objects | as_ggchord_seq/as_ggchord_feature/as_ggchord_link: explicit IDs/full lengths to frames | [Adapters](formats.html#adapters) |
| Validate and repair | validate_ggchord_data/clean_ggchord_data/prepare_ggchord_data: issue versus repair reports | [Workflow](data.html#validation) |
| Genome/contig keys | validate_ggchord_hierarchy consumes explicit membership | [Hierarchy](data.html#hierarchy) |
| Prepare links | filter/deduplicate/merge/bundle_ggchord_ribbons retains membership/estimates/reports | [Preparation](ribbons.html) |
| Select targets/windows | find_ggchord_loci/focus_ggchord_data/segment_ggchord_data/map_ggchord_position: explicit selection/reversible mapping | [Real target](feature-geometry.html#locus-selection), [windows](views.html) |
| Homology and coding projection | link_ggchord_clusters/build_ggchord_cds/project_ggchord_features/project_ggchord_links require membership/CDS/phase evidence | [Homology](homology.html), [projection](projection.html) |
| Search annotations | primer_catalog/find_common_features/find_primer_bindings/find_restriction_sites/filter_restriction_sites separate catalogue, matches and display set | [Recognition](recognition.html), [primers](primers.html), [restriction](restriction-sites.html) |
| Preview and consume | view_ggchord/diagnose_ggchord/export_ggchord_layout/get_chord_layout/optimize_ggchord_layout inspect target dimensions and named results | [Output](output.html) |




## Exercises and replacement input {#adapt}

Change only diagnostic width and inspect measurement hints. Focus the same track to another real source interval and inspect source metadata. Retain unmodified source tables and read reports after each step. Consume named layout fields and never treat build-local layer_id as a permanent biological identifier.
