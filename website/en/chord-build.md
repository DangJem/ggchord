Add one alignment and one CDS to the same chord plot. Run the blocks in page order.

## Start with four backbones {#bundled}

<span id="backbones"></span>Use the four sequences from the [previous section](chord-start.html). The code below recreates the four backbones .

{{code:frontbackbone}}
{{figure:frontbackbone}}

## Add an alignment ribbon {#links}

Take the first row of the bundled alignment table. It connects two of the four sequences.

{{code:frontlinks,frontribbon}}

{{figure:frontribbon}}

Follow the ribbon to both ends. It spans 26,298–51,139 bp on `MT108731.1` and 7,121–31,959 bp on `MT118296.1`. Each end covers an interval, so this layer uses `geom_link_ribbon()`.

The default colour reflects the alignment identity for this row, about 98.6%. The crossed shape comes from how the interval ends are joined. Use the alignment coordinates to assess inversion.

## Add a CDS {#add-genes}

<span id="prepare"></span><span id="inspect-genes"></span><span id="genes"></span>Select the original annotation starting at 43,765 bp on `OQ646790.1` and recover its biological strand.

{{code:frontgenes,frontfeature}}

{{figure:frontfeature}}

Find the red arrow on green `OQ646790.1`. It marks the plus-strand integrase CDS from 43,765 to 46,140 bp. The ribbon and this CDS use independent input tables.

<span id="labels"></span>To display the CDS name, add `geom_gene_label(data = chord_genes_example)`. It reads the label from `anno`. Next, use the same backbones to compare [rotation and gaps](chord-tuning.html).
