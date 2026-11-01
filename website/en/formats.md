## Adapt an existing table {#adapters}
Run the [shared setup](index.html#quick-start). These examples rename columns of
bundled source-backed sequences, links and genes; they do not invent alignments.

{{code:adapters}}

Adapters return ordinary frames with `ggchord_conversion` provenance. IDs remain
unchanged. Supply reliable full lengths rather than using maximum feature ends.
See [alignment input](alignments.html) for minimal fields and coordinate units.

## BED and PAF: convert coordinates once {#zero-based}
The BED example writes known pUC19 intervals as minimal BED3 and reads them back.
PAF comes directly from real minimap2 output in the public tutorial, with original
optional tags/values. The display retains mapq ≥30 and receives only six endpoint roles.

{{code:formats}}
{{figure:paf}}

The reader converts BED zero-based half-open intervals to one-based inclusive
coordinates; do not add one again. BED3 has no strand, so it returns unknown
strand; use a strand-bearing format when needed. Negative PAF relative direction
becomes decreasing subject endpoints. These are valid reverse hits, not values
to sort indiscriminately. Quality/identity use producer fields, not new alignments.

## GenBank: retain record identity {#genbank}
{{code:genbank}}

This is the original NCBI L09137.2 record, named pUC19c; do not rename it to the
bundled pUC19. Results include `seq_data`, `features$annotation` and `report`,
retaining segmented locations and record/feature provenance. The selected Python
needs Biopython already installed; no automatic installation occurs.
`unsupported = "report"` records unsupported locations rather than guessing bounds.

## Optional interval objects {#ranges}
{{code:optionalranges}}
{{results:optional}}

An originally unknown-strand pUC19 annotation remains unknown as `*`. Known full
lengths populate Seqinfo. `GInteractions` can become links but still needs explicit
sequence/length context. Missing dependencies are marked skipped, not validated
by source presence. See the [data contract](reference/ggchord_data_contract.html).
SAM/BAM CIGAR is not parsed directly: use an alignment-aware tool before table
adaptation; renaming columns alone cannot resolve splicing.
