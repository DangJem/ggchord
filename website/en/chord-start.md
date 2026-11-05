Start with four bundled sequences. Read their positions and lengths before adding links and features.

## Draw the backbones {#draw}

Run this lesson in a fresh R session after installing the version described on the [installation page](installation.html).

```r
library(ggchord)
```

{{code:frontbackbone}}

{{figure:frontbackbone}}

## Read the plot {#read}

Each arc represents a row of the sequence table. Find `MT108731.1` by its label and compare it with `MT118296.1`. These two sequences are 64,323 bp and 32,090 bp, so the first arc is about twice as long as the second.

`geom_seq()` draws the backbones, `geom_seq_label()` adds their IDs, and `coord_chord()` arranges the backbones as arcs.

## Try removing the labels {#try}

Remove the `geom_seq_label()` line and run the code again. The arcs stay in place and the sequence names disappear.

## Use your own sequences {#adapt}

Replace `seq_data_example` with a table containing `accver` and `length`. See [Basic input tables](chord-data.html#tables) for the column definitions. Next, [add links and features](chord-build.html).
