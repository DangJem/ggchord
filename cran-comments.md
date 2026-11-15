## Resubmission

This is a revised submission of ggchord 0.28.1. The version number is
retained; the help examples and DESCRIPTION have been revised since the
previous submission.

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

The previously reported plotting examples have been reduced to small,
explicit subsets of the packaged real data. Feature examples restrict
sequence discovery to two named DNA features, and the primer example uses
reference binding sites instead of searching the full primer catalogue.
Core plotting examples remain executable during ordinary checks. Additional
complete label/layout demonstrations are retained in donttest blocks.

The VCF and BigWig examples require optional Bioconductor decoders whose
initial namespace loading can be expensive. They remain documented, guarded
by dependency/platform checks, and are now in donttest blocks. This does not
claim an improvement in the decoder implementations. We understand that
CRAN may also execute donttest examples.

### Validation status

Targeted local execution and timing of the affected examples is recorded in
design-roadmap/versions/v0.28.1.md in the source repository. This maintenance
pass does not establish a new full R CMD check result. The revised ordinary plotting examples completed locally in 0.48-2.27
seconds per help topic. Optional decoder and additional layout examples
were also executed manually. The extended layouts still emit an existing
physical-geometry refinement warning; this was absent from the ordinary
examples and has not been fixed in this documentation revision. Updated
Windows and Debian incoming results remain to be obtained before claiming
that the reported timing NOTE is resolved on those machines.
