## Membership is supplied evidence {#membership}
After the [shared setup](index.html#quick-start), this declared synthetic cluster
connects one feature on each sequence. It is a teaching assertion, not an
inference from the repeated gene annotation text.

{{code:homologs}}
{{figure:homologs}}

Inspect `homologs$report$comparisons` and `$members`. Unknown members are omitted
and reported; ambiguous paralogs default to omission. `paralogs = "all"` explicitly
requests bounded Cartesian products. Repeated feature IDs require an explicit
segment policy, and an outer span is an anchor, not an exon alignment.
`find_ggchord_loci()` selects exact feature IDs and expands in biological direction;
unknown direction requires symmetric extensions. Continue with [focus](views.html).

## Project a declared protein correspondence {#projection}
{{code:linkprojection}}
{{results:projection}}

Equal protein endpoint spans are not evidence of ungapped alignment. The explicit
`ungapped = TRUE` assertion is required. The helper splits where either CDS crosses
a coding-segment boundary and preserves source protein metadata. It does not emit
protein identity as nucleotide identity or invent nucleotide alignment length.
For domains and partial/phase policies, see [CDS projection](projection.html).
