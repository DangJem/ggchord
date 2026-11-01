## Derive membership from actual coding sequences {#membership}
After the [shared setup](index.html#quick-start), select AmpR annotations from
bundled pUC19 and pETDuet-1, extract actual DNA and correct for strand. Group only
identical oriented CDS sequences; sharing the label AmpR is not itself evidence.

{{code:homologs}}
{{figure:homologs}}

This strict example demonstrates identical CDS, not a requirement that all
homologs be identical. General memberships should come from reliable external
analysis. Inspect `homologs$report$comparisons` and `$members`; unknown members
and ambiguous paralogs are reported. `paralogs = "all"` requests bounded Cartesian
products explicitly. A segmented feature's outer span is an anchor, not an exon
alignment. `find_ggchord_loci()` selects exact IDs and extends biologically.

## Project an evidenced protein interval to both CDS {#projection}
{{code:linkprojection}}
{{results:projection}}

The first 20 residues of identical CDS justify `ungapped = TRUE`. Equal spans in
an ordinary BLASTp HSP do not. Projection splits at CDS segment boundaries and
retains protein metadata without inventing nucleotide identity. For a real
BLASTp query-span example see [alignment input](alignments.html#protein); see
[coordinate projection](projection.html) for the feature workflow.
