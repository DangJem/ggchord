<div class="task-intro"><p>Combine original FASTA, GFF3 and standard BLAST tables into a traceable comparison. Inspect inputs, choose a display subset, then compare chord and full-length views.</p></div>

## Prepare this page {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) so your working directory contains `examples/`.

## Read and inspect original records {#read}
{{code:realimport}}
{{console:realimport}}

FASTA supplies IDs and lengths, GFF3 supplies annotation ranges and strands, and BLAST supplies observed alignments. The empty validation summary means no issues were reported for these inputs.

## Make the display selection explicit {#selection}
{{code:realprepare}}
{{console:realprepare}}

Keep alignments with at least 80% identity and an alignment length of at least 500 bp. For the chord view, show up to six of the longest CDS per sequence so the links remain readable.

{{code:realfeatures}}
{{output:realfeatures}}
{{results:real}}

## Build a chord comparison {#chord}
{{code:realchord}}
{{figure:realchord}}

Read sequence labels first, then source-strand feature arrows and interval links.
Ribbons encode existing aligned intervals and mapped identity; their crossing is
not itself evidence of a rearrangement. The selected features remain visibly
separate from the much wider alignment spans. Use the full input tables when
scientific interpretation needs annotations outside this display subset.

## Inspect all annotations in a linear view {#linear}
{{code:reallinear}}
{{figure:reallinear}}

This view draws all imported annotations. Ribbons connect adjacent displayed
sequences by default; retaining every comparison is an explicit layer decision.
The bp scale bar belongs to collinear coordinates. See [local views](views.html)
for cropping, source-coordinate axes and independent windows, and
[output](output.html) to preview and save at the gallery's stated dimensions.

## Keep the evidence with the figure {#provenance}
The input archive contains unmodified FASTA/GFF3/BLAST files, a SHA256 manifest
and the shared example script. The plot manifest records the code hash, PNG device,
width, height and dpi. The accession names identify source records; a drawing is
not a replacement for the original annotations or alignment evidence.

## Select fields before passing data {#minimal-fields}
Backbones use `accver` and `length`. Ribbons use IDs and endpoints on both sequences. This example also filters on `pident` and `length`. See [alignment inputs](alignments.html).

## Check your result {#check-result}
Check that IDs agree across the tables, inspect filtering reports and verify that guides encode identity and strand. Use the overview to locate regions, then focus a neighbourhood to study a CDS.

<div class="next-step">Continue: <a href="homology.html">Explicit homologs and coding-coordinate links →</a></div>

## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

```r
ggsave("phage-comparison.png", p_real_chord, width = 10, height = 8, units = "in", dpi = 144)
```
