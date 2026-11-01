<p class="lesson-intro">Start with a plot you can interpret, then choose the next step for your data and question.</p>

## First time: follow this sequence {#start}
1. [Draw your first plot](installation.html): install, select minimal inputs, draw, interpret and save.
2. [Choose a plot](modes.html): distinguish overall links, single-sequence annotation and local position comparisons.
3. [Understand data and coordinates](data.html): check complete lengths, one-based closed intervals and strands.

## Finish a real task {#workflows}
| Goal | Preparation | What you will produce |
| --- | --- | --- |
| [Compare four sequences](real-data.html) | FASTA, GFF3 and standard BLAST tables from the input bundle | A filtered chord comparison and full-length overview |
| [Build a complete plasmid map](plasmids.html) | Bundled pUC19 DNA | A layered circular map with primers and selected enzyme sites |
| [Compare corresponding CDS neighbourhoods](homology.html) | Two bundled plasmids and their coding DNA | Supported membership, local windows and a source-coordinate plot |
| [Compute and add GC tracks](tracks.html) | Bundled pUC19 DNA | Bars, a line and a category band with units, windows and value domains |

Each page provides its own preparation code. For raw-file workflows, [download the input bundle](../assets/tutorial-inputs.zip) and retain its `examples/` paths. Bundled-data workflows need no download.

## Have a plot? Solve one problem {#tasks}
| Problem | Read |
| --- | --- |
| Different column names or coordinate conventions | [Import and conversion](formats.html), [nucleotide and protein alignments](alignments.html) |
| Overlapping annotations or crowded text | [Feature stacking](annotations.html), [label layouts](labels.html) |
| Local windows, segments or display orientation | [Focused views](views.html), [per-sequence parameters](parameters.html) |
| Too many links; filtering or summarization needed | [Prepare links](ribbons.html) |
| Unclear colours, guides or visual hierarchy | [Style mappings](styling.html), [common adjustments](adjustments.html) |
| Protein intervals need genomic coordinates | [CDS projection](projection.html) |
| Export or interpretation of diagnostics | [Output and preview](output.html), [FAQ](faq.html) |

## Look up concepts and exact arguments {#reference}
Use the [function reference](reference.html) for precise contracts, the [glossary](glossary.html) for concepts and the [gallery](gallery.html) to find tutorials by outcome. For older API code, consult [migration](migration.html) and the [changelog](news.html).

## Capability index {#coverage}
{{coverage}}

This is a lookup index, not a required learning sequence. Consult function reference pages for exact returns and restrictions.
