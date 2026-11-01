<div class="task-intro"><p>Combine original FASTA, GFF3 and standard BLAST tables into a traceable comparison. Inspect inputs, choose a display subset, then compare chord and full-length views.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed. For file-based examples, first extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains `examples/`.

{{prepare:docstyle}}

This complete workflow uses four sequence records already present in the source
repository. It requires the matching development checkout or the tutorial input
bundle, not new alignments or external gene prediction. [Download the inputs](../assets/tutorial-inputs.zip)
and unpack them so `examples/` is under the working directory. Record the package
version and [build manifest](../build-manifest.json) with the exported figure.

## 1. Read and inspect original records {#read}
{{code:realimport}}

The FASTA headers supply exact accession versions and lengths. GFF3 contributes
CDS intervals and their original strands. BLAST output contributes existing
alignments; the reader does not run BLAST. Inspect `real_validation` before
continuing: decreasing link endpoints can be legitimate reverse-direction hits.
Never turn unknown annotation strands into `+` to make a gene shortcut accept them.

## 2. Make the display selection explicit {#selection}
{{code:realprepare}}
{{results:real}}

The link display retains identity ≥80% and reported alignment length ≥500 bp.
The chord view selects the six longest imported CDS per sequence, labelling at most two selected CDS with informative annotations (excluding
hypothetical protein). These are display rules, not biological findings or a claim
that omitted genes are absent. Raw `feature_real` and `link_real` remain intact;
`real_prepared$report` and `real_filtered$report` retain preparation/filter evidence.
No display strands are recoded. This differs from the small packaged teaching
fixture whose documented generation includes display-strand recoding.

## 3. Build a chord comparison {#chord}
{{code:realchord}}
{{figure:realchord}}

Read sequence labels first, then source-strand feature arrows and interval links.
Ribbons encode existing aligned intervals and mapped identity; their crossing is
not itself evidence of a rearrangement. The selected features remain visibly
separate from the much wider alignment spans. Use the full input tables when
scientific interpretation needs annotations outside this display subset.

## 4. Inspect all annotations in a linear view {#linear}
{{code:reallinear}}
{{figure:reallinear}}

This view draws all imported annotations. Ribbons connect adjacent displayed
sequences by default; retaining every comparison is an explicit layer decision.
The bp scale bar belongs to collinear coordinates. See [local views](views.html)
for cropping, source-coordinate axes and independent windows, and
[output](output.html) to preview and save at the gallery's stated dimensions.

## 5. Keep the evidence with the figure {#provenance}
The input archive contains unmodified FASTA/GFF3/BLAST files, a SHA256 manifest
and the shared example script. The plot manifest records the code hash, PNG device,
width, height and dpi. The accession names identify source records; a drawing is
not a replacement for the original annotations or alignment evidence.

## Select fields before passing data {#minimal-fields}
The `*_source` objects keep original fields and provenance; `seq_real` has two
columns, `feature_real` only geometry/type/text, and `link_real` adds identity and
reported length only because this workflow filters and colours by them. Standard
12-column BLAST conversions are used, not a requirement for the historical 17 fields.
The bundled gene fixture preserves its original strand in `source_strand`; the
short examples restore that biological strand before selecting four gene columns.
See [alignment input](alignments.html) for exact format and unit boundaries.

## Check your result {#check-result}
Check that IDs agree across the tables, inspect filtering reports and verify that guides encode identity and strand. Use the overview to locate regions, then focus a neighbourhood to study a CDS.

<div class="next-step">Continue: <a href="homology.html">Explicit homologs and coding-coordinate links →</a></div>

## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

```r
ggsave("phage-comparison.png", p_real_chord, width = 10, height = 8, units = "in", dpi = 144)
```
