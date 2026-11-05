Start with four bundled sequences in horizontal rows. Use the same input as in the chord chapter and look at the backbones first.

## Draw four backbones {#draw}

Run this lesson in a fresh R session after installing the version described on the [installation page](installation.html).

```r
library(ggchord)
```

{{code:frontlinearbase}}

{{figure:frontlinearbase}}

## Read position and length on a common scale {#read}

Each row represents one sequence. Compare `MT108731.1` at about 64 kb with `MT118296.1` at about 32 kb. A shared bp scale lets you compare horizontal lengths directly.

## Try removing the sequence names {#try}

Remove `geom_seq_label()`. The four backbones and their relative lengths remain visible.

## Use your own sequences {#adapt}

Supply a sequence table with `accver` and full `length`. Next, [add an alignment and a CDS](collinear-build.html).
