{{code:courseblockaw}}

{{prepare:setup}}

Extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains examples/. Each section names its optional dependencies.

## Continuous signals: bedGraph and BigWig {#signals}

This installed-package example contains actual GC percentages computed from the unchanged bundled pUC19 sequence in 100 bp windows (last window 86 bp). It is sequence composition, not measured read coverage. `data-raw/prepare-signal-inputs.R` records the calculation and `inst/extdata/signal-sources.json` records hashes and the explicit original-header mapping. Readers preserve IDs.

{{code:courseblockax}}

bedGraph uses four tab-separated columns, zero-based half-open coordinates and finite signed values. Import converts coordinates once; original endpoints, values, physical lines and opaque display headers remain traceable. One file represents one signal; split multiple track definitions. Literal sequence IDs `track` and `browser` are valid. Headers never set plot styles. Overlaps remain rows, without an implicit last-wins rule.

Rectangles use `start - 1` and `end` as cell boundaries; no line is interpolated across absent ranges. BigWig requires optional **rtracklayer**, whose decoder does not support Windows. Use an existing bedGraph file there. Binary GRanges coordinates are already one-based inclusive. Imported runs retain binary floating-point precision; run counts are not original measurement counts. Full lengths come from the file header. Query windows form an unstranded union, results are clipped to it, and empty selection returns zero rows. Missing intervals stay absent, distinct from a stored zero; no coverage, identity, support, units or pre-clipping observation IDs are invented. Each Track's mapping, value range, units and Guide remain explicit.

## Aligned reads: an optional coverage bridge {#aligned-reads}

Use optional **GenomicAlignments**, **Rsamtools** and **GenomicRanges** for an existing BAM. Keep filtering and counting choices explicit; ggchord draws the resulting signal. This example uses the dependency's unchanged Drosophila RNA-seq subset. Its source BAM header retains the original TopHat/reference and SRR provenance; the fixed source URL and checksum are in `tools/long-term-audit/v0.27.0-read-coverage-results.json`. No read alignment or variant calling occurs.

{{code:courseblockay}}

Each retained **alignment record** contributes one over its aligned reference blocks. Skipped regions (`N`) and, with this explicit choice, deletions (`D`) contribute zero; supplementary alignments are retained. No base-quality threshold or overlapping-mate/molecule deduplication is applied. These values are not automatically fragment counts, unique molecules, expression, Identity or variant support. Zero runs are stored zeros for this declared computation, unlike missing experimental intervals. For large inputs use the dependency's streaming `BamFile` coverage method with the same explicit filters.

Pass these one-based-inclusive ranges to `geom_track(geom = "rect")` using `xmin = start - 1`, `xmax = end`, `ymin = 0`, `ymax = value`, as in [continuous signals](#signals). Full reference lengths come from `scanBamHeader(bam_file)[[1]]$targets`, not the last read end. Specify source windows explicitly before displaying a small region.

For CRAM, decode using a mature installed HTSlib/pysam or samtools environment and the **actual matching reference** when required, retaining the original header, filters and conversion log. Never invent a reference, auto-download unrequested resources or claim all CRAM files are independently decodable. The directed prototype uses `no_ref=1` solely to re-encode the same real BAM records without external reference access, then verifies all record fields/tags, reference dictionary and coverage after decoding. It is a format round trip, not an original CRAM observation or universal reference-based CRAM certification. The upstream header's `SO:sorted` and absent temporary index generate retained decoder messages; sequential reading succeeds without rewriting them.
