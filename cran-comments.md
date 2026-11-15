## Resubmission

This is a revised submission of ggchord 0.28.1. The version number is
retained; the help examples and DESCRIPTION have been revised since the
previous submission. Existing website input records are also bundled for
reproducible plotting examples, retaining original values and plotting defaults.

### Version jump

The last CRAN version is 0.16.6. The submitted version includes cumulative
work developed since that CRAN release, including additional sequence views,
input adapters, cross-panel compositions and sequence-detail displays.
Intermediate development versions were not all submitted to CRAN. The version
jump from 0.16.6 to 0.28.1 is intentional. NEWS.md describes the cumulative changes and migration information.

### DESCRIPTION terminology

The flagged words are established terms in sequence analysis:

* HSP means high-scoring segment pair, an aligned segment pair reported by
  sequence alignment tools. The abbreviation is now expanded in DESCRIPTION.
* dotplots are plots of pairwise sequence alignment coordinates.
* synteny describes conserved relationships/order of genomic regions.

These words are intentional, not spelling mistakes.

### Example timings

The plotting examples now use small, explicit subsets of real packaged data.
The Ribbon help topic draws one connection between two source sequences in
one plot. The Adaptive example uses only two backbones and sequence names.
Complete research examples remain on the website.

The VCF example previously parsed all 10,376 records in the optional decoder's
example file. It now copies the original headers and first three unmodified
records to a temporary VCF, reads it and prints the resulting three-row table.
No extra VCF data file is distributed. The BigWig example performs one regional
read of the existing packaged signal. Both are ordinary examples with optional
dependency checks; BigWig retains its Windows platform guard.

In separate fresh local R processes (R 4.6.0, macOS arm64), VariantAnnotation
1.58.0 first loaded in 2.914 s; reading the complete VCF took 4.723 s and a
six-row preview took 0.001 s. With the three-record input, first loading took
2.824 s, preparation 0.020 s, reading 0.109 s and output 0.001 s.
rtracklayer 1.72.0 first loaded in 2.183 s; full BigWig reading took 0.035 s,
a regional read 0.033 s and preview output was below the timer's resolution.
Namespace and source-package loading are separated from reader timings.
Filesystem caches were not flushed. These are local observations, not new
Windows or Debian incoming results.

### Validation status

Direct execution of the revised help topics succeeded. Whole-topic elapsed
times locally were 1.701 s for Ribbon, 0.955 s for Adaptive, 3.097 s for VCF
and 2.254 s for BigWig. The two reader topics were each executed in a separate
R process with their decoder namespace initially unloaded.

This maintenance pass does not establish a new full R CMD check result.
Updated Windows and Debian incoming results remain to be obtained before
claiming that the timing NOTE is resolved on those machines. R CMD check
--as-cran normally also executes donttest examples unless configured otherwise;
remaining extended-layout donttest examples are not claimed to be skipped by
CRAN or fixed by this reader/example revision.
