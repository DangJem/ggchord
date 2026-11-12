# Source-derived continuous signals

Canonical tiny inputs are packaged once in `inst/extdata/pUC19-gc-100.bedGraph`
and `inst/extdata/pUC19-gc-100.bw`, with `signal-sources.json` provenance.
Regenerate them using `data-raw/prepare-signal-inputs.R` and optional development
dependencies declared by that script. They contain actual GC percentages from
unchanged `examples/plasmid/pUC19.fna`, not experimental coverage. The full
string must match the bundled fixture before explicitly mapping the original
header `pUC19.dna` to fixture ID `pUC19`; readers never perform this renaming.

Format references: [UCSC bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html),
[UCSC BigWig](https://genome.ucsc.edu/goldenPath/help/bigWig.html),
[rtracklayer manual](https://bioconductor.org/packages/release/bioc/manuals/rtracklayer/man/rtracklayer.pdf).
No third-party implementation is copied. BigWig decoding uses optional
rtracklayer and is unsupported on Windows; the text input remains available.
