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

On 2026-11-15, the current source was built with R CMD build and the resulting
archive checked with R CMD check --as-cran on R 4.6.0, macOS 27.2 arm64.
HTML Tidy 5.8.0 was selected from the existing local installation. No check
settings were used to disable incoming checks, examples, donttest examples,
tests or vignette rebuilding.

Result: 0 ERRORs, 0 WARNINGs, 1 NOTE. The only NOTE is the intentional minor
version jump from 0.16.6 to 0.28.1, explained above. Both ordinary examples
and the --run-donttest pass completed without a timing NOTE. All 34 testthat
assertions passed, both vignette outputs rebuilt, and PDF/HTML manual checks
passed. The final recorded full-example pass took 0.248 s for Ribbon, 0.834 s
for Adaptive, 0.514 s for VCF and 2.277 s for BigWig; its longest help topic
was geom_feature_label_repel at 4.874 s. These are sequential check timings;
decoder dependencies may already be loaded and they are not cold-load benchmarks.

Two vignette URLs pointing at an undeployed website chapter were found and
fixed before the final build. Installed tutorials now refer to the bundled
dataset help for provenance; the complete website chapter remains in its source.

Existing example output still includes PDF-device substitution of an en dash
and one physical-geometry refinement warning in the extended Fusion example.
They did not produce an R CMD check WARNING or NOTE, and this maintenance
revision does not claim to fix them.

Updated Windows and Debian incoming results remain to be obtained before
claiming that the timing NOTE is resolved on those machines. No package has
been resubmitted as part of this local validation.

### Source archive size

The rebuilt archive is 5,066,598 bytes (4.832 MiB), compared with 4,862,281
bytes for the original 0.28.1 archive: an increase of 204,317 bytes (4.20%).
The main addition is six compressed, ready-to-plot real-data files totalling
201,736 bytes; plotting no longer requires preparing those inputs from the
repository. Developer website sources and data-generation tools are excluded.
The installed size reported by check is 8.7 Mb (informational, not a NOTE).
