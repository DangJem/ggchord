## Start with source-backed CDS {#cds}
Run the [identical-CDS comparison](homology.html#membership) first. It uses real
bundled plasmid annotations and DNA. Complete CDS use explicit phase 0 and
translation table 11.

{{code:projection}}

This example selects residues 2–20 for display; it does not claim a detected
protein domain. Projection returns genomic intervals with strand/segment context.
Actual domains require exact parent CDS IDs and amino-acid ranges; annotation
names or maximum feature endpoints cannot supply missing evidence.

## Segments, phase and partial regions {#boundaries}
`build_ggchord_cds()` constructs a segment map; `project_ggchord_features()` then
splits amino-acid intervals through it. Crossing a boundary can return several
rows; multiplying coordinates by three is insufficient. Inspect policies/reports
for partial CDS, phase and truncation. DNA matching, interval projection and
actual translation/alignment are distinct operations.

## BLASTp query spans and two-ended links {#links}
{{code:linkprojection}}

This ungapped case relies on [verified identical coding DNA](homology.html).
Real BLASTp HSPs may contain gaps: query spans can still be projected, but paired
links need the subject CDS and justified ungapped correspondence. See
[real BLASTp import/projection](alignments.html#protein). This version does not
reinterpret aa as bp or add native protein-coordinate plotting.
