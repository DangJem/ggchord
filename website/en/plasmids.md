Run the [shared setup](index.html#quick-start) for the explicit palette first.

## Use a source-backed sequence {#input}
The bundled `plasmid_example_pUC19` is copied from the repository's original
`examples/plasmid/pUC19.fna`; bases and sequence identity are not edited for display.
The workflow uses the package's source-backed recognition catalogues. It needs no
network request or experimental measurements fabricated for the figure.

{{code:plasmid}}
{{figure:plasmid}}
{{results:plasmid}}

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
