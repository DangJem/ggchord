<div class="task-intro"><p>Build a readable plasmid map from bundled pUC19 DNA, with separate annotation, primer and restriction layers. Keep display choices separate from complete results.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed.

{{prepare:featurepalette,docstyle}}

The preparation above declares the palette and shared appearance.

## Use a source-backed sequence {#input}
The bundled `plasmid_example_pUC19` is copied from the repository's original
`examples/plasmid/pUC19.fna`; bases and sequence identity are not edited for display.
The workflow uses the package's source-backed recognition catalogues. It needs no
network request or experimental measurements fabricated for the figure.

{{code:plasmidinput}}

Inspect `plasmid_features`, `plasmid_primers` and `display_sites` to distinguish complete search results from display choices.

## Compose independent layers {#draw}
Add the backbone, annotations, labels, primers and restriction sites in that order. Features and their labels share one Position; each layer reads its own data.

{{code:plasmid}}
{{figure:plasmid}}
{{results:plasmid}}

Mixed exterior annotations share root-ordered placement, including sparse restriction labels beside primers. This example uses the dedicated plasmid theme and default feature fitting: internal text is tried before adjacent/external placement. Inspect at the final export dimensions instead of hiding crossings by shrinking text or forcing offsets.


The original annotation table is retained in `plasmid_feature_source`. Dedicated
primer layers replace duplicate `primer_bind` annotations in the displayed feature
table; only geometry, type and annotation text are passed to feature layers.

## Read the annotation layers {#layers}
The double backbone establishes one circular sequence. General features use a
plasmid preset and explicit stacking Position; labels use the same Position.
Primer arcs indicate binding ranges, restriction ticks indicate cleavage/display
positions, and the centre label identifies the sequence and length. A cropped or
filtered display does not change the full-sequence restriction count.

## Choose the displayed sites deliberately {#sites}
`find_restriction_sites()` returns all matches. This figure explicitly selects
EcoRI, BamHI and HindIII so that a compact restriction cluster remains readable;
`plasmid_sites` retains the complete search result. The `reference` filter uses a saved
reference display set only when the complete supplied sequence matches that
reference. Other tasks can select `unique`, `six_plus` or explicit enzyme names.
Inspect `display_sites` and the site's methylation/source columns before drawing.
The reference set is a presentation choice, not a claim about an unknown sample.

## Understand catalogue and label limits {#limits}
`primer_names` exposes catalogue entries. This example searches the universal catalogue exactly; reference primer
annotations and exact oligonucleotide searches have distinct provenance. Near-exact protein feature
recognition has the limitations in [recognition](recognition.html); it does not
validate biological function or a primer experiment. Text can move, shrink within
its explicit minimum, or hide. Keep the geometry and inspect exported label
reasons at the final size rather than assuming every annotation must have text.

## Adapt the workflow {#adapt}
For your own named DNA sequence, replace the input and inspect search results
before adding geoms. For externally curated annotations, pass your feature table
instead of running recognition. Use [stacking and shapes](annotations.html) for
geometry, [labels](labels.html) for fitting, and [output](output.html) for devices.

`geom_primer_label_repel(side = "inside")` keeps callouts inside the backbone,
even alongside exterior restriction and feature labels. This is corrected in the
unreleased 0.16.6.9000 checkout; v0.16.5 can move those requested inside labels
outside. The default outside placement and existing source coordinates are unchanged.

Inside restriction labels also retain their inner geometry and now report the
inner side in the exported annotation registry; earlier versions mislabelled
those diagnostic rows as exterior resources.

The full inside scene can still place long primer text over other feature regions.
The default `max_overlaps = Inf` does not certify zero overlap across feature
layers. Prefer the default outside field when space is tight, or use
`show_location = FALSE` to shorten the display; inspect the result at its final size.

## Check your result {#check-result}
Start with the centre name and length, then read positions, category colours, arrow directions and outer labels. Check enzyme names in display_sites; selecting a display subset leaves plasmid_sites intact.

<div class="next-step">Continue: <a href="labels.html">Labels and local placement →</a></div>

## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

```r
ggsave("pUC19-map.png", p_plasmid, width = 10, height = 9, units = "in", dpi = 144)
```

## Complete example: a dense pETDuet-1 map {#complex}

Use the pUC19 example above to learn the layers, then move to the richer pETDuet-1 map. Inputs come from the bundled sequence, feature and reference-primer catalogues. The restriction subset retains unique six-base-or-longer recognition sites from commercial_all; it is not every cut made by every enzyme.

{{code:showcaseplasmiddata}}

{{code:showcaseplasmid}}

{{figure:showcaseplasmid}}

Inputs contain 20 features, five reference primers and 58 grouped cut positions. The feature preset retains plum coding regions, mint resistance genes, yellow origins and blue regulatory features; primers are purple, consistently across this site. Dense labels depend on export size; open the full image to read them.

The explicit viewport keeps the dense map at the documented export size. `label_fit = "shrink"` allows short-feature labels to fit down to 2.5 mm; labels that do not fit remain callouts. Preview and save at the same 16.5 × 8 in size.

```r
ggsave("pETDuet-1-map.png", p_show_plasmid, width = 16.5, height = 8, units = "in", dpi = 144)
```
