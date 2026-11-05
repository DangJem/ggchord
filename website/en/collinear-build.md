Use the same four sequences and add one alignment and one CDS. Run the blocks in order.

## Recreate the four backbones {#bundled}

<span id="prepare"></span><span id="switch"></span><span id="backbones"></span>

{{code:frontlinearbase}}
{{figure:frontlinearbase}}

## Add an aligned interval pair {#links}

{{code:frontlinks,frontlinearlink}}

{{figure:frontlinearlink}}

The ribbon spans 26,298–51,139 bp on the first row and 7,121–31,959 bp on the second. Its horizontal shift reflects the position of each interval on its own sequence.

## Inspect the CDS and add an arrow {#inspect-genes}

{{code:frontgenes,frontselectedgene}}

{{output:frontselectedgene}}

This plus-strand integrase CDS spans 43,765–46,140 bp on the third row. Its original annotation is in `anno`.

### Add it to the plot {#add-genes}

{{code:frontlinearcds}}

{{figure:frontlinearcds}}

<span id="labels"></span>The red arrow appears in the right half of the third, green backbone. Add `geom_gene_label(data = chord_genes_example)` to show its name. Next, [adjust row spacing and add a scale bar](collinear-tuning.html).
