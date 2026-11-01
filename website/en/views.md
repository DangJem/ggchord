## Focus without losing origins {#focus}
After the [shared setup](index.html#quick-start):

{{code:focus}}
{{figure:focus}}

Focus relocates retained positions to local one-based coordinates and retains
source origins, truncation flags and a report. `boundary = "trim"` can clip
links by a linear display projection; it does not realign sequences. A scaled
alignment length is flagged as estimated. `coordinates = "source"` restores
source axis labels, including after repeated focusing.

## Independent windows and local reversal {#segments}
{{code:segment}}

Each window has a reversible source mapping; omitted gaps have no mapping.
`map_ggchord_position()` maps explicit source positions into a plot for ordinary
ggplot2 annotations. Supply the window ID when windows share a source sequence.
Do not interpret display orientation as a change to the original strand.

## Explicit biological relationships {#homology}
`find_ggchord_loci()` selects declared feature/name queries. `link_ggchord_clusters()`
requires an explicit membership table and sequence pairs; equal gene names do
not establish homology. `filter_ggchord_ribbons()`, `deduplicate_ggchord_ribbons()`
and `merge_ggchord_ribbons()` prepare existing alignments, preserving source
membership and estimated-length flags. Bundle and density statistics summarize
provided links rather than creating new sequence alignments.

## Map an annotation back to source coordinates {#mapping}
{{code:mappingposition}}

Run the [focus example](views.html#focus) first. The helper returns Cartesian x/y
for a source position 100 bp inside the first selected window. It does not make all ordinary ggplot2 coordinates
become genomic bp. In segmented views choose the window ID, and check positions
in omitted gaps rather than extrapolating across them. Use [explicit homologs](homology.html)
when the locus selection comes from declared feature/cluster evidence.
