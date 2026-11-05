Identify binding intervals on actual DNA before choosing arc and label placement. Computational matches are not wet-lab validation, and catalogue direction hints do not determine the binding strand on the target.

## Catalogue, search and actual input {#input}

Independently load pUC19 and retain the first three exact universal-catalogue hits for clear display. Results retain primer name, oligonucleotide, one-based interval, strand, match method and uniqueness. Select after searching; absence from this display subset is not absence from the full result.

{{code:comparefinishdnabase}}
{{code:comparefinishprimer}}
{{code:comparefinishprimerrows}}
{{output:comparefinishprimerrows}}


## Primer arc placement and width {#primer-arcs}

Default between follows the backbone centre line. Outside/inside use arc_offset (default 0.045) along the local normal. The between/outside/inside examples change only side. The wider-arc example returns to between and changes only arc_linewidth from 0.026 to 0.05. This local geometry width is not bp or text size.

{{code:compareprimerbetween}}
{{figure:compareprimerbetween}}
{{code:compareprimeroutside}}
{{figure:compareprimeroutside}}

Relative to outside placement, change only arc_offset from 0.045 to 0.09:

{{code:compareprimeroffset}}
{{figure:compareprimeroffset}}
{{code:compareprimerinside}}
{{figure:compareprimerinside}}
{{code:compareprimerwidth}}
{{figure:compareprimerwidth}}


## Primer names, range text and side {#primer-labels}

Labels default to outside with show_location=TRUE; leaders originate at the directional arc’s 3′ tip. Compare baseline with only location text hidden, or only side changed to inside. Hiding text does not shorten binding intervals. Interior text can overlap other features; default max_overlaps=Inf does not guarantee collision-free mixed layers.

{{code:compareprimerlabel}}
{{figure:compareprimerlabel}}
{{code:compareprimernolocation}}
{{figure:compareprimernolocation}}
{{code:compareprimerlabelinside}}
{{figure:compareprimerlabelinside}}


## Search versus drawing choices {#parameters}

`primer_catalog(set="all", aliases=FALSE)` returns deduplicated oligonucleotides; universal selects common records and aliases=TRUE shows source names for shared sequences. `find_primer_bindings()` computes matches on DNA; circular searches include origin-spanning matches, while linear targets need circular=FALSE. Exact three_prime suffix matching is not a general mismatch thermodynamic model. Wallace temperature is a composition-based heuristic.

Primer geometry defaults: colour=#A020F0, side=between, arc_offset=0.045, arc_linewidth=0.026. Text defaults: side=outside, show_location=TRUE, max_overlaps=Inf. Both require the same target sequence; renaming accver does not transfer bindings to another DNA. Use single circular mode; aa inputs reject nucleotide-only layers.



## Practice and your primers {#adapt}

Change only arc_offset, then return to defaults and hide location text. Provide actual oligonucleotides and names, inspect match fields before selecting displayed hits. See the [parameter atlas](parameter-atlas.html) for complete current defaults and [restriction sites](restriction-sites.html) or the [full plasmid workflow](plasmids.html) next.
