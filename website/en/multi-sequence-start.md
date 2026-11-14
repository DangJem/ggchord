Build backbones, a link, a feature and labels in one section using bundled source records. Run in order; this section uses ggchord and optional dplyr.

## Draw the backbones {#draw}

Run this section in a fresh R session after following the [installation instructions](installation.html).

{{code:courseblockac}}

{{code:frontbackbone}}

{{figure:frontbackbone}}

## Read the plot {#read}

Each arc represents a row of the sequence table. Find `MT108731.1` by its label and compare it with `MT118296.1`. These two sequences are 64,323 bp and 32,090 bp, so the first arc is about twice as long as the second.

`geom_seq()` draws the backbones, `geom_seq_label()` adds their IDs, and `coord_chord()` arranges the backbones as arcs.

## What the label layer controls {#try}

Remove the `geom_seq_label()` line and run the code again. The arcs stay in place and the sequence names disappear.

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

To display the CDS name, add `geom_gene_label(data = chord_genes_example)`. It reads the label from `anno`. 

## Display the feature name {#labels}

{{code:frontgenelabels}}

{{figure:frontgenelabels}}

## Switch to a linear layout {#collinear}

{{code:courseblockad}}

![Collinear](../assets/frontlabelledlinear.png)

This retains the ribbon, CDS and feature-name layers; Collinear selects adjacent rows by default. Compare all four views in the [next section](multi-sequence-views.html).

## Use your own sequences {#adapt}

Replace `seq_data_example` with a table containing `accver` and `length`. See [Basic input tables](data.html#prepare) for the column definitions.

For local settings, see [Flexible parameter formats](parameters.html).
