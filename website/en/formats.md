## Start with canonical roles {#adapters}
Run the [shared setup](index.html#quick-start). A renamed column is an explicit
role mapping, not a request to rewrite accession versions or infer full lengths.

{{code:adapters}}

Inspect the ordinary returned tables and `ggchord_conversion` attributes.
`as_ggchord_seq()` requires known full lengths; `as_ggchord_feature()` preserves
unknown strands; `as_ggchord_link()` preserves anchors and available metadata.
The same adapters are used by supported direct S4 inputs to layers.

## BED and PAF: convert once {#zero-based}
These files are deliberately synthetic. They demonstrate format conversion,
not a measured alignment or annotation of a real sequence.

{{code:formats}}

The BED interval `[0,100)` becomes one-based inclusive `1..100`. PAF negative
relative strand yields decreasing subject endpoints; mapq 255 is unavailable,
retained separately from the converted missing value. PAF identity is derived
from producer-reported counts, not recomputed from bases. Keep import attributes
and source-line columns; do not subtract or add one again before drawing.

## GenBank: records, segments and provenance {#genbank}
{{code:genbank}}

This optional example runs only if the selected Python already has Biopython.
The function installs nothing. Its result contains `seq_data`, named
`features$annotation` and `report`; joined locations retain real segments and
source record/feature indices. Unknown lengths are not inferred from feature
ends. Use `unsupported = "report"` to retain unsupported locations in a report,
not to flatten them into guessed local intervals.

## Optional interval objects {#ranges}
{{code:optionalranges}}
{{results:optional}}

The unknown strand `*` remains unknown. Known seqlengths come from Seqinfo,
not the largest observed interval. InteractionSet `GInteractions` additionally
adapts to links; keep a common region/length context and explicit biological
membership. Missing optional packages are reported as skipped in the example
record, not as validated support. See the complete [data contract](reference/ggchord_data_contract.html).
