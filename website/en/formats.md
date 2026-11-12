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


## Existing alignments: keep columns and source positions separate {#aligned-characters}

These examples require the source checkout's `examples/alignment/` directory. It contains unchanged Biopython DNA/protein alignment inputs, the upstream permission notice, hashes and a reproducible FASTA re-encoding rule. No aligner is run. Full record IDs and original characters remain unchanged; choose rows and columns explicitly.

`read_aligned_fasta()` reads one existing alignment. Declare DNA, RNA or AA; equal lengths alone do not establish homology. `alignment_column` includes `-` and `.` gaps, whose `source_position` is NA. Default positions count the ungapped input fragment from one; they are not complete-genome coordinates. Lowercase letters consume positions, so this reader must not be used as an A3M insertion parser.

```r
aligned_dna <- read_aligned_fasta(
  "examples/alignment/opuntia-aligned.fasta", "DNA", source_file = TRUE
)
selected_rows <- unique(aligned_dna$accver)[c(1, 3, 7)]
local_dna <- aligned_dna |>
  filter(accver %in% selected_rows, between(alignment_column, 321, 380))

map_alignment_position(local_dna, selected_rows[1], c(342, 343, 344))
map_alignment_position(local_dna, selected_rows[1], 343, inverse = TRUE)

p_characters <- ggplot2::ggplot(
  local_dna, aes(alignment_column, accver, label = symbol)
) +
  geom_residue(aes(fill = symbol), size = 3, overflow = "error") +
  ggplot2::scale_fill_manual(
    values = c(
      A = "#B8D8BA", C = "#A6CBE3", G = "#F2D28C", T = "#EAA8A8",
      N = "#D8D8D8", "-" = "#F5F5F5", "." = "#F5F5F5"
    ),
    guide = "none"
  ) +
  ggplot2::labs(x = "Alignment column", y = NULL) +
  ggplot2::theme_minimal()
p_characters
```

The first selected row has a gap at alignment column 343. Column lookup therefore returns `matched = TRUE`, `is_gap = TRUE` and a missing source position. Inverse lookup of source position 343 finds its actual alignment column, without interpolation. Positions removed by filtering return `matched = FALSE`; filtering never renumbers retained columns. `map_alignment_position()` selects an alignment-row ID even when a separately declared source ID differs.

`geom_residue()` is a standard Cartesian ggplot layer. Styles and user scales follow ggplot2. It measures glyphs against actual cell width and height at draw time: `overflow = "error"` requires readable letters; `"tile"` retains colour cells where text cannot fit, without shrinking fonts. Focus fewer columns or increase physical width to show letters. For example, save `p_characters` with ordinary `ggsave(width = 12, height = 3.4, units = "in")`; simply raising DPI adds no physical room. Discrete rows, identity/reverse position scales and Cartesian coordinates are supported; this layer does not attach characters to ggchord's backbone Coords or Fusion panels.

### Protein excerpt {#aligned-protein}

The same data contract accepts actual protein characters. This excerpt contains upstream gaps and residue differences; fill distinguishes literal symbols, not inferred functional domains.

```r
aligned_protein <- read_aligned_fasta(
  "examples/alignment/protein-aligned.fasta", "AA"
)
protein_rows <- unique(aligned_protein$accver)[2:4]
local_protein <- aligned_protein |>
  filter(accver %in% protein_rows, between(alignment_column, 221, 260))

ggplot2::ggplot(local_protein, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 3, overflow = "tile") +
  ggplot2::scale_fill_hue(l = 85, c = 35, guide = "none") +
  ggplot2::labs(x = "Alignment column", y = NULL) +
  ggplot2::theme_minimal()
```

### Explicit source coordinates {#aligned-source}

Provide one `source_data` row per aligned record with `accver`, `start` and `strand`. `start` is the first nongap source coordinate. Negative-strand positions decrease; strings must already be in alignment orientation. An optional `source_accver` identifies the original molecule. If you provide a full forward `sequence`, every nongap character is checked, using nucleotide complementation on the negative strand and order reversal only for AA. Missing source strings remain unverified; no origins or strands are guessed.

The following example checks consistency against ungapped fragments derived from the same alignment. It is deliberately **not independent genomic validation**. For a real genomic projection, supply independently known source sequences and interval metadata instead.

```r
original_strings <- read_fasta_lengths(
  "examples/alignment/opuntia-aligned.fasta", keep_sequence = TRUE
)
fragment_sources <- original_strings |>
  transmute(
    accver, source_accver = accver, start = 1, strand = "+",
    sequence = gsub("[-.]", "", sequence)
  )
checked_cells <- as_ggchord_alignment(
  original_strings, "DNA", source_data = fragment_sources
)
```

No mismatch/gap detail can be recovered from six HSP endpoints alone. This workflow accepts existing columns, not CIGAR/A3M/Clustal/Stockholm parsing, CDS translation/phase inference, consensus/Logo analysis or alignment computation. It leaves source identity and positions available for later explicit macro/detail mapping.
