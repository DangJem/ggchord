## Tables {#tables}
| Input | Required fields | Meaning |
| --- | --- | --- |
| Sequence | `accver`, `length` | Unique ID; positive full length in bp |
| Interval link | `qaccver`, `saccver`, `qstart`, `qend`, `sstart`, `send` | Query and subject anchors |
| Gene | `accver`, `start`, `end`, `strand` | Known `+`/`-`; optional `anno` |
| Feature | `accver`, `start`, `end`, `strand` or declared direction | May retain unknown strand; optional `type`, `label`, `feature_id` |

Identifiers may be arbitrary strings. Accession versions are preserved exactly.
`pident` and alignment `length` are optional unless a mapping, filter or statistic
needs them. A feature interval does not imply a sequence's full length.

## Coordinate contract {#coordinates}
Genomic coordinates are **one-based, inclusive bp**. BED/PAF readers convert their
source conventions; do not convert twice. Gene endpoints are ordered with
strand stored separately; link endpoint order may encode direction. A one-bp
link without direction evidence stays unknown. Unknown strand is never silently
changed to `+`; use `geom_feature()` instead of the strict gene layer.

Amino-acid coordinates require [explicit CDS projection](projection.html).
Native Cartesian layout units, genomic bp, angles in degrees, physical mm and
track value units are different domains. See the full
[data contract](reference/ggchord_data_contract.html).

## Import original records {#import}
Run this repository-only example from its root. Original files are not rewritten.

{{code:imports}}

`read_genbank()` retains record/feature provenance and segmented locations.
`read_bed()` and `read_paf()` retain their conversion metadata. `as_ggchord_seq()`,
`as_ggchord_feature()` and `as_ggchord_link()` adapt supported S4 objects; provide
full lengths and IDs explicitly when the object lacks them. `file` selects one
file; `files` supports an explicitly expanded collection. Missing sequence or
homology data must be supplied by the analyst, not invented by the plotting layer.

## Validate, clean, prepare {#validation}
After the [shared setup](index.html#quick-start):

{{code:validation}}

Validation reports issues; cleaning makes explicit repairs with a report.
Preparation combines workflow steps without inferring biological membership.
Inspect `validation`, `cleaned$report` and `prepared$report` before plotting.
`source_row` in geometry refers to the resolved layer input, while `.source_line`
means a physical file line; GenBank has separate record/feature provenance.
Aggregates may have `source_rows` and `.length_estimated` flags. Preserve these
when filtering or merging links.
