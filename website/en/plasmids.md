## Research question and sequence identity {#question}

The pUC19 case first asks how complete actual DNA, identified annotations, primer bindings and restriction sites become separate layers. The dense pETDuet-1 case asks how to preserve hierarchy at a fixed export size. Full lengths come from their own original FASTA records, not maximum annotation ends; renaming pUC19c does not make it pUC19.

Retain separate feature catalogue methods/provenance, primer catalogue versus exact matches, restriction database and methylation background. Curated reference annotations are not experimental validation for an unknown sample. Saved reference sets require exact complete-sequence matching; inspect search results for other targets.
## Use a source-backed sequence {#input}
The bundled `plasmid_example_pUC19` is copied from the repository's original
`examples/plasmid/pUC19.fna`; bases and sequence identity are not edited for display.
The workflow uses the package's source-backed recognition catalogues. It needs no
network request or experimental measurements fabricated for the figure.

{{code:plasmidinput,frontplasmidsummary}}

{{output:frontplasmidsummary}}

Inspect `plasmid_features`, `plasmid_primers` and `display_sites` to distinguish complete search results from display choices.
## Compose independent layers {#draw}
Add the backbone, annotations, labels, primers and restriction sites in that order. Features and their labels share one Position; each layer reads its own data.

{{code:plasmid}}
{{figure:plasmid}}
{{results:plasmid}}

Start with the displayed cut-site count. The layer receives only three selected sites. The full search result remains in `plasmid_sites` for looking up sites outside the figure.

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
This example searches the universal catalogue exactly; reference primer
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

`geom_primer_label_repel(side = "inside")` places primer callouts inside the backbone. Interior restriction labels appear inside the backbone and are recorded there in the exported annotation registry.

The full inside scene can still place long primer text over other feature regions.
The default `max_overlaps = Inf` does not certify zero overlap across feature
layers. Prefer the default outside field when space is tight, or use
`show_location = FALSE` to shorten the display; inspect the result at its final size.
## Check your result {#check-result}
Start with the centre name and length, then read positions, category colours, arrow directions and outer labels. Check enzyme names in display_sites; selecting a display subset leaves plasmid_sites intact.
## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

{{code:courseblocku}}
## Evidence checks before saving {#evidence-check}

Preview with view_ggchord at the declared width/height and save with ggsave at matching dimensions. Check IDs, source strands, endpoint directions, selection reports, legend units and hidden-label reasons. A successful PNG does not certify font portability across PDF/SVG viewers. Keep complete source tables, display subsets, parameters, version and input hashes. Download the [asset/provenance record](../assets/ASSET-SOURCES.md).
## Source rights and reproducibility {#source-rights}

The package code is MIT-licensed; source database records and reference catalogues retain their original terms. The package license does not relicense third-party data. Keep accession identities and provenance with exported figures. The [asset/source record](../assets/ASSET-SOURCES.md), input archive hashes and displayed preparation script document this example; algorithmic matches remain hypotheses for experimental verification.
