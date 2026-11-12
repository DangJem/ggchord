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

No mismatch/gap detail can be recovered from six HSP endpoints alone. This workflow accepts existing columns, not CIGAR/A3M/Clustal/Stockholm parsing, phase inference, consensus/Logo analysis or alignment computation. Explicit CDS translation is described below; exact macro/detail registration is available in the Fusion lesson.

## Source-backed CDS translation {#cds-translation}

`build_ggchord_cds()` defines observed coding coordinates; `translate_ggchord_cds()` additionally needs actual complete forward DNA/RNA strings and optional **Biostrings**. This real EMALE negative-strand CDS has source phase 0 and an ATG start. Table 11 is an explicit choice checked against the unchanged upstream protein FASTA. Input provenance is in `examples/format-inputs/README.md`.

```r
library(dplyr)
# Requires the source checkout, dplyr and optional Biostrings.
genome <- read_fasta_lengths(
  "examples/format-inputs/emales-BVI_023A.fna", keep_sequence = TRUE
)
coding <- read_gff3("examples/format-inputs/emales-BVI_023A.gff3") |>
  filter(feature_id == "BVI_023A_0009")
cds <- build_ggchord_cds(coding, genetic_code = 11)
observed <- translate_ggchord_cds(cds, genome, initiation = "start")
expected <- read_fasta_lengths(
  "examples/format-inputs/protein-pair.faa", keep_sequence = TRUE
) |>
  filter(accver == "BVI_023A_0009")
stopifnot(sub("[*]$", "", observed$models$protein_sequence) == expected$sequence)
observed$residues |>
  select(alignment_column, symbol, codon, source_positions) |>
  head()
```

`models` reports coding/protein strings, genetic code, initiation/ambiguity policies, partial ends and terminal remainder. `nucleotides` retains each genomic base; `residues$source_positions` contains the three genomic positions of each codon in biological order, even across exon boundaries. Residue positions use the CDS feature ID's observed-protein namespace; they are not an inferred external protein accession. Terminal and internal stops remain `*`. The example removes only the final stop for comparison with the upstream FASTA convention.

The default `initiation = "literal"` does not recode alternative starts; `"start"` requires a complete 5-prime end. `ambiguous = "error"` rejects fuzzy codons; explicit `"solve"` delegates to Biostrings, using X if unresolved. Unknown codes, missing/short sources and declared translation exceptions fail. `partial = "trim"` belongs to the CDS map: only the initial phase and incomplete terminal codon are excluded; later segment phases validate continuity without deleting junction bases. No missing sequence is reconstructed. RNA must be declared explicitly. Ordered `segments` support the existing spliced/origin-join contract.

A complete runnable nucleotide/codon/macro example is `examples/m2-detail.R`; its negative-strand labels are centred on their three coding bases. Real spliced plus/minus controls, source hashes and original expected translations are retained in `examples/cds/`. [Register the character plot](fusion.html#character-detail) to connect its exact source boundaries.

## Continuous signals: bedGraph and BigWig {#signals}

This installed-package example contains actual GC percentages computed from the unchanged bundled pUC19 sequence in 100 bp windows (last window 86 bp). It is sequence composition, not measured read coverage. `data-raw/prepare-signal-inputs.R` records the calculation and `inst/extdata/signal-sources.json` records hashes and the explicit original-header mapping. Readers preserve IDs.

```r
signal_file <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
signal <- read_bedgraph(signal_file, source_file = TRUE)
data(plasmid_example_pUC19)
gc_band <- position_track("GC", side = "inside", width = unit(7, "mm"))
p_signal <- ggchord(plasmid_example_pUC19) +
  geom_seq(show.legend = FALSE) +
  geom_track(
    aes(xmin = start - 1, xmax = end, ymin = 0, ymax = value, fill = value),
    data = signal, geom = "rect", track = gc_band
  ) +
  scale_track_value_continuous("GC", limits = c(0, 100), unit = "GC (%)") +
  scale_track_fill_gradientn(
    "GC", colours = c("#F1E51D", "#31849B"), limits = c(0, 100), name = "GC (%)"
  ) +
  coord_circular()
p_signal

if (.Platform$OS.type != "windows" && requireNamespace("rtracklayer", quietly = TRUE)) {
  indexed_file <- system.file("extdata", "pUC19-gc-100.bw", package = "ggchord")
  selected_signal <- read_bigwig(
    indexed_file,
    which = data.frame(accver = "pUC19", start = c(50, 100), end = c(150, 250)),
    source_file = TRUE
  )
  selected_signal |> select(accver, start, end, value)
}
```

bedGraph uses four tab-separated columns, zero-based half-open coordinates and finite signed values. Import converts coordinates once; original endpoints, values, physical lines and opaque display headers remain traceable. One file represents one signal; split multiple track definitions. Literal sequence IDs `track` and `browser` are valid. Headers never set plot styles. Overlaps remain rows, without an implicit last-wins rule.

Rectangles use `start - 1` and `end` as cell boundaries; no line is interpolated across absent ranges. BigWig requires optional **rtracklayer**, whose decoder does not support Windows. Use an existing bedGraph file there. Binary GRanges coordinates are already one-based inclusive. Imported runs retain binary floating-point precision; run counts are not original measurement counts. Full lengths come from the file header. Query windows form an unstranded union, results are clipped to it, and empty selection returns zero rows. Missing intervals stay absent, distinct from a stored zero; no coverage, identity, support, units or pre-clipping observation IDs are invented. Each Track's mapping, value range, units and Guide remain explicit.

## Variant records and structural breakpoints {#variants}

Install optional **VariantAnnotation** for VCF and additionally **StructuralVariantAnnotation** for explicit structural decoding. These examples read the dependencies' unchanged source-backed examples: the chr22 1000 Genomes subset and a real caller VCF distributed with StructuralVariantAnnotation. Data stays in those packages; no third-party parser implementation or input is copied into ggchord. `check-variant-import.R` records exact package versions and source hashes.

```r
if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
  variant_file <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
  imported_variants <- read_vcf(variant_file, source_file = TRUE)
  imported_variants$features$variant |>
    select(accver, start, end, feature_id, REF, ALT) |>
    head()
}

if (requireNamespace("StructuralVariantAnnotation", quietly = TRUE)) {
  structural_file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
  structural <- read_vcf(structural_file, structural = TRUE)
  structural$links |>
    select(qaccver, saccver, qstart, qend, sstart, send,
      q_breakend_orientation, s_breakend_orientation)
  structural$report$records[[1]]$warnings
}
```

Each placed record remains one feature row with an ALT list, including multiallelic calls. Features show the reference allele span, including the padding base, rather than an inferred symbolic SV footprint. Complete INFO, sample/genotype arrays and headers remain in `report$records[[1]]$parsed`. Inspect actual caller fields there; neither missing support nor Identity becomes zero or a fabricated count. Unknown full lengths remain NA; provide reliable reference lengths before plotting. Telomeric/out-of-bounds spans stay in `unplaced_features` and the source, without clipping.

Structural decoding is opt-in, keeps confidence intervals, disables absent-mate inference, and creates one link per reciprocal pair with both original record indices. Warnings identify records not paired by the decoder, while all records remain in the report. `.plottable = FALSE` marks breakends outside known reference bounds; affected pairs remain in `unplaced_links`. Exact single-base breakpoints suit `geom_link_line()` after explicitly selecting their positions; uncertainty intervals suit interval geometry. `q_breakend_orientation` and `s_breakend_orientation` describe the derivative sequence, not relative homology strand. Never manufacture alignment Identity to colour these links.

For modern BCF (`.bcf` or `.bcf.gz`), supply `python` naming an existing interpreter with optional **pysam**: `read_vcf("calls.bcf", python = "/path/to/environment/bin/python")`. Only BCF invokes it, transcodes to temporary VCF via HTSlib, and records decoder/version/messages. Virtual-environment interpreter paths retain their environment. Original record order, values and genotypes are preserved; transcoded text is not an original physical BCF line. No installation occurs. This reader loads selected local files into memory; select or subset large datasets upstream.
