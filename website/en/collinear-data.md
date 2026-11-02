Collinearity extends ggchord's original chord grammar into a linear view. The design goal is to retain shared data and layers, switching coordinates wherever the contracts permit.

## The same three tables {#tables}
| Input | Required columns | Linear interpretation |
| --- | --- | --- |
| seq_data | accver, length | One sequence per row with a shared bp x scale |
| ribbon_data | qaccver, saccver, qstart, qend, sstart, send | An observed interval correspondence between rows |
| gene_data | accver, start, end, strand | Genes separated above/below rows by a Position |

Optional fields and the historical BLAST naming follow [1.1](chord-data.html). Identity needs pident; labels need anno. Source fields preserve coordinates in local windows; display gaps do not represent invented DNA lengths.

## Obtain the inputs {#import}
Reuse the FASTA/GFF3/BLAST workflow in [1.2](chord-import.html), or import PAF, GRanges and GInteractions through [adapters](formats.html). Protein alignments require [explicit CDS projection](projection.html); aa endpoints must not be presented as bp.

## Current conversion boundaries {#contract}
The [next page](collinear-build.html) demonstrates replacing coord alone. Explicit circular radius, curvature, ring and radial sequence-label placement are rejected by the current linear contract rather than silently ignored. Keep coordinate-specific styling in the appropriate showcase. Perfect conversion of every possible composition remains a design goal, not a current guarantee.

Linear views default to adjacent displayed rows. Set `adjacent_only = FALSE` in shared layers to retain the chord link set. Choose manual sequence directions or coordinate-owned automatic orientation, not both at once.
