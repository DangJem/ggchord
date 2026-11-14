Start with observed HSPs from one sequence pair, then extend to a four-sequence matrix. The minimal inputs are complete sequence lengths and four HSP endpoints; annotations are optional.

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

## Multi-sequence HSP matrix {#matrix}

These inputs come from original repository FASTA and standard BLAST tables. All 32 HSPs across four sequences are retained; self-alignments, reciprocal records and missing HSPs are not invented. The 266 CDS belong to the shared source preparation and are not submitted to the Dotplot segment layer. An empty area means no supplied HSP, not proven absence of homology.

{{prepare:dotplotmatrixdata}}

{{code:flagshipdot}}

{{manualfigure:flagdotplot}}

## Inspect source direction {#direction}

Relative endpoint direction comes from the source record. Reversing a display axis changes the visible slope while leaving the source direction unchanged. Count the supplied directions per sequence pair:

```r
hsp_direction <- flagship_hsp |>
  mutate(relative_direction = if_else(
    (qend - qstart) * (send - sstart) < 0, "reverse", "forward"
  )) |>
  count(qaccver, saccver, relative_direction, name = "HSP")
hsp_direction
```

## Query an exact position {#mapping}

Map the first HSP query start onto the concatenated x axis. Inverse lookup uses `inverse = TRUE` and an explicit axis; display gaps return NA. The four endpoints remain unchanged.

```r
map_ggchord_position(
  flag_dotplot, flagship_hsp$qaccver[1], flagship_hsp$qstart[1], axis = "x"
)
```

## Boundaries and next steps {#support}

Backbone gene/feature/track/ribbon layers, facets and Fusion endpoints are unsupported. To inspect gaps or mismatches, provide actual aligned characters; see [the alignment workflow](aligned-sequences.html). Example: trace four endpoints of one segment and explain slope and colour separately.
