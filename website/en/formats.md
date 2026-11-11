<div class="task-intro"><p>Convert tables, BED, PAF or GenBank into ggchord inputs while retaining record identity and converting coordinate conventions only once.</p></div>

## Choose by the files you have {#choose-format}

| Files | Start here |
| --- | --- |
| FASTA sequences, GFF3 annotations and BLAST alignments | [Draw a comparison from source files](real-data.html#read) |
| BED intervals or PAF (Pairwise mApping Format) alignments | [Check zero-based coordinates](#zero-based) |
| GenBank sequence with annotations | [Preserve record IDs and segments](#genbank) |
| Existing R tables | [Map required fields](#adapters) |

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed. For file-based examples, first extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains `examples/`.

{{prepare:setup}}

## Adapt an existing table {#adapters}
The preparation above loads the input tables. These examples rename columns of
bundled source-backed sequences, links and genes; they do not invent alignments.

{{code:adapters}}

{{console:adapters}}

Adapters return ordinary frames with `ggchord_conversion` provenance. IDs remain
unchanged. Supply reliable full lengths rather than using maximum feature ends.
See [alignment input](alignments.html) for minimal fields and coordinate units.

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

## GenBank: retain record identity {#genbank}
{{code:genbank}}
{{console:genbank}}

This is the original NCBI L09137.2 record, named pUC19c; do not rename it to the
bundled pUC19. Results include `seq_data`, `features$annotation` and `report`,
retaining segmented locations and record/feature provenance. The selected Python
needs Biopython already installed; no automatic installation occurs.
`unsupported = "report"` records unsupported locations rather than guessing bounds.

The GenBank example was rerun with Python Biopython 1.85 enabled: it returned L09137.2, length 2686 bp and one annotation row. For a separate documentation environment, install Biopython with that Python’s `python3 -m pip install biopython`, then put its executable on PATH before starting R. Do not treat an optional-dependency skip or NULL as successful import.

## Optional interval objects {#ranges}

Convert one interval and then an interval pair. Full lengths come from `seqlengths`. Install GenomicRanges for the first example and InteractionSet for the second.

{{code:frontimportgr}}

{{output:frontimportgr}}

{{code:frontimportgi}}

{{output:frontimportgi}}

## Paired BED intervals {#bedpe}

`read_bedpe()` reads six required tab-separated fields, optional name/score/two strands and additional metadata. The first anchor stays query and the second stays subject. Scores remain strings; neither Identity nor read support is inferred. The format follows the [BEDPE description](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format).

This example re-encodes the intervals of bundled observed alignments for a coordinate round trip; it is not paired-read evidence. Six columns contain no strand evidence, so the imported anchor strands and relative direction remain unknown. Supply real strand columns when known; do not infer them from ascending coordinates.

```r
data(ribbon_data_example)
paired_rows <- ribbon_data_example |>
  transmute(
    chrom1 = qaccver, start1 = pmin(qstart, qend) - 1, end1 = pmax(qstart, qend),
    chrom2 = saccver, start2 = pmin(sstart, send) - 1, end2 = pmax(sstart, send)
  )
paired_file <- tempfile(fileext = ".bedpe")
write.table(paired_rows, paired_file,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
paired_links <- read_bedpe(paired_file)
unlink(paired_file)
paired_links |>
  select(qaccver, saccver, qstart, qend, sstart, send, q_mapped, s_mapped)
```

Negative-strand anchors use decreasing canonical endpoints, consistent with paired-interval adapters; the original start/end and strand fields are also retained. Names are not assumed unique. Extra columns remain in `bedpe_extra`, named `column_11` onwards. Unknown chromosomes (`.`) and boundaries (`-1`) fail by default. Explicit `unmapped="keep"` retains incomplete rows for auditing; filter `q_mapped & s_mapped` before adaptation or plotting. Invalid known coordinates and zero-width intervals always fail. For interval ribbons without Identity, set a constant fill explicitly rather than manufacturing pident.

## Check your result {#check-result}
Inspect returned fields and conversion reports before selecting plotting fields. Preserve endpoint direction for reverse PAF hits; do not rename a GenBank record to pretend it is a different bundled plasmid.

<div class="next-step">Continue: <a href="alignments.html">Import nucleotide and protein alignments →</a></div>


## Retain actual FASTA characters {#characters}

Use `keep_sequence = TRUE` to retain characters with their physical header line and optional file provenance. Only whitespace is removed; case is retained. Declare bp or aa explicitly instead of guessing the alphabet. Focusing crops the string as well as the interval; local position 1 below is original position 1001. A decreasing locus selects the same interval and does not reverse-complement it. Lengths alone cannot supply characters; HSP endpoints cannot supply mismatch or gap columns. Gapped FASTA counts gap characters and needs a separate alignment-to-source map before genomic use.

```r
character_source <- read_fasta_lengths(
  "examples/plasmid/pUC19.fna",
  sequence_unit = "bp", keep_sequence = TRUE, source_file = TRUE
)
character_locus <- focus_ggchord_data(
  character_source,
  loci = data.frame(accver = character_source$accver, start = 1001, end = 1060)
)
character_locus$seq_data |>
  select(accver, length, sequence, .source_origin)
```
