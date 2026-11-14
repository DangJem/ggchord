A single-sequence circular map starts with two columns. Feature, primer and restriction searches additionally need real DNA; existing annotations can be drawn directly without running recognition again.

<a id="tables"></a>

## Tables and optional information {#tables}
| Input | Required content | Optional content |
| --- | --- | --- |
| Backbone | accver, length; exactly one sequence | label for display; sequence for searches |
| Features | accver, start, end and declared directionality | type, anno/label, feature_id and shape categories |
| Primers | Computed binding interval, strand/name roles required by the geom | Recognition and provenance fields |
| Restriction sites | Position/name roles required by the geom | Enzyme catalogue, methylation assumption and filter provenance |

Backbones do not require DNA. Prefer search-result tables for the latter layers and consult their help for detailed fields. Unknown-strand annotations belong in generic feature layers.

<a id="inspect"></a>

## Inspect bundled real records {#inspect}
{{code:courseblockai}}

The compact names retain ggchord's early BLAST workflow history. `accver` need not be a database accession. Full DNA and length must agree. Recognition is computational matching, not experimental validation.

<a id="import"></a>

## Import your own plasmid {#import}
`read_fasta_lengths()` reads IDs and full lengths only. Searches additionally need actual DNA, for example read with optional Biostrings `readDNAStringSet()` and `as.character()`, retaining matching IDs. GFF3/BED/GenBank or GRanges annotations follow the [format guide](formats.md), retaining provenance. Handle circular-origin and segmented records according to the reader contract; do not invent one continuous interval by merging every segment.

<a id="path"></a>

## A progressive case sequence {#path}

Build pUC19 step by step in the [layer section](circular-start.md#input), then complete the [pUC19 and pETDuet-1 cases](plasmids.md). Consult [search assumptions](recognition.md) and the [13 reference maps](plasmid-catalog.md) as needed. The bundled pUC19c alias refers to the same example.
