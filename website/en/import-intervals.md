{{code:courseblockao}}

{{prepare:setup}}

Extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains examples/. Each section names its optional dependencies.

## BED and PAF: convert coordinates once {#zero-based}
BED starts at zero and excludes the right endpoint. Read the pUC19 origin interval and compare the coordinates before and after import.

{{code:formats}}
{{output:formats}}

PAF gives the full sequence lengths and alignment spans. Select alignments with mapq ≥ 30, then collect both sequences and their lengths for plotting.

{{code:formatpaf}}
{{output:formatpaf}}

{{code:drawpaf}}
{{figure:paf}}

A reverse alignment has descending subject coordinates. Use the converted coordinates directly for plotting.
## Optional interval objects {#ranges}

Convert one interval and then an interval pair. Full lengths come from `seqlengths`. Install GenomicRanges for the first example and InteractionSet for the second.

{{code:frontimportgr}}

{{output:frontimportgr}}

{{code:frontimportgi}}

{{output:frontimportgi}}
## Paired BED intervals {#bedpe}

`read_bedpe()` reads six required tab-separated fields, optional name/score/two strands and additional metadata. The first anchor stays query and the second stays subject. Scores remain strings; neither Identity nor read support is inferred. The format follows the [BEDPE description](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).

This example re-encodes the intervals of bundled observed alignments for a coordinate round trip; it is not paired-read evidence. Six columns contain no strand evidence, so the imported anchor strands and relative direction remain unknown. Supply real strand columns when known; do not infer them from ascending coordinates.

{{code:courseblockap}}

Negative-strand anchors use decreasing canonical endpoints, consistent with paired-interval adapters; the original start/end and strand fields are also retained. Names are not assumed unique. Extra columns remain in `bedpe_extra`, named `column_11` onwards. Unknown chromosomes (`.`) and boundaries (`-1`) fail by default. Explicit `unmapped="keep"` retains incomplete rows for auditing; filter `q_mapped & s_mapped` before adaptation or plotting. Invalid known coordinates and zero-width intervals always fail. For interval ribbons without Identity, set a constant fill explicitly rather than manufacturing pident.
