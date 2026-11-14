<div class="task-intro"><p>Project known protein intervals to genomic positions through real CDS, retaining phase, strand and segmentation explicitly.</p></div>

## Prepare this page {#prepare}
Use the [verified AmpR input](homology.html) or save and run the preparation script below.

{{recipe:ampr}}

## Start with source-backed CDS {#cds}
Run the [identical-CDS comparison](homology.html#membership) first. It uses real
bundled plasmid annotations and DNA. Complete CDS use explicit phase 0 and
translation table 11.

{{code:projection}}

{{console:projection}}

This example selects residues 2–20 for display; it does not claim a detected
protein domain. Projection returns genomic intervals with strand/segment context.
Actual domains require exact parent CDS IDs and amino-acid ranges; annotation
names or maximum feature endpoints cannot supply missing evidence.

## Segments, phase and partial regions {#boundaries}
`build_ggchord_cds()` constructs a segment map; `project_ggchord_features()` then
projects amino-acid intervals through it. Each input feature returns one row;
its genomic pieces are retained in the `segments` list-column. Two-ended link
projection can instead split into several paired rows. Multiplying coordinates
by three is insufficient. Inspect policies/reports
for partial CDS, phase and truncation. DNA matching, interval projection and
actual translation/alignment are distinct operations.

## BLASTp query spans and two-ended links {#links}
{{code:linkprojection}}

{{console:linkprojection}}

This ungapped case relies on [verified identical coding DNA](homology.html).
Real BLASTp HSPs may contain gaps: query spans can still be projected, but paired
links need the subject CDS and justified ungapped correspondence. See
[real BLASTp import/projection](protein-alignments.html#projection). Genomic plots use the projected bp coordinates.

## Result checks {#check-result}
Check parent CDS IDs, projection reports and returned segments. The selected residue range is a display example, not a detected domain.

<div class="next-step">Continue: <a href="alignments.html">Import nucleotide and protein alignments →</a></div>
