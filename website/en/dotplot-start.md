{{code:courseblockam}}

## Plot original HSP endpoints {#dotplot}

The source-backed EMALE PAF pair has 11 observed HSPs. Its declared sequence
lengths supply both axes. `read_paf()` converts zero-based half-open intervals
to the internal one-based convention and reverses target endpoints for minus
alignments. The files and their provenance are in the [input download](formats.html).
The following code needs the optional `dplyr` package for data preparation.

{{code:dotplotdata}}

Each line runs from `(qstart, sstart)` to `(qend, send)`. `geom_seq()` is not
required. Use `link_*` mappings and scales for colour, alpha, width and type.
The default Identity domain is 0–100 with ticks 0/50/80/90/95/100; an explicit
user scale overrides it.

{{code:dotplotview}}

{{figure:viewdotplot}}

The legacy constructor may warn about descending endpoints in these valid reverse HSPs. Keep the normalized ordering; sorting each interval would erase the displayed inversion.

`x_seqs` and `y_seqs` independently select/order accessions. `gap` is display
spacing in native bp/aa, not omitted source sequence. Axis directions accept
1/-1 or complete named vectors. Display reversal changes the slope without
changing the biological source table. Each axis starts at source position 1;
concatenated segments retain their own accession and local tick labels.

{{code:dotplotreverse}}

Use `map_ggchord_position(..., axis = "x"/"y")` for native/global conversion,
with `inverse = TRUE` to go back. Display gaps return NA. Focused source
origins use `source = TRUE`; default axis labels stay window-relative.
Unknown directions are warned, omitted and reported; contradictory direction
metadata is rejected. Repeated declared HSP IDs require `allow_repeated=TRUE`
in all participating layers. Unidentified duplicate/self rows remain distinct.
Without cg/cs/CIGAR or another alignment trace, these are HSP segments, not
base-level mismatches or gap paths.

## Boundaries and next steps {#support}

Backbone gene/feature/track/ribbon layers, facets and Fusion endpoints are unsupported. To inspect gaps or mismatches, provide actual aligned characters; see [the alignment workflow](aligned-sequences.html). Exercise: trace four endpoints of one segment and explain slope and colour separately.
