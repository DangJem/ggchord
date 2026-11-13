<div class="task-intro"><p>Read the default BLAST columns in observed alignments and select the six endpoint fields needed for plotting.</p></div>

## Why begin with BLAST input? {#prepare}

ggchord was originally designed to visualize BLASTn alignments as ribbons between nucleotide sequences. BLAST tables are its default recommended alignment input, and the link fields retain query/subject naming as the plotting tools have expanded.

Plotting uses full sequence lengths and precomputed alignments. FASTA supplies backbone lengths and BLAST supplies corresponding intervals. Begin with the default BLAST table, then choose nucleotide or protein coordinates.

`-outfmt 6` writes tab-separated rows without a header. `-outfmt 7` adds comment lines beginning with `#`. Omitting `-outfmt` produces human-readable pairwise text. These two commands illustrate nucleotide and protein output.

```sh
blastn -query query.fna -subject subject.fna -outfmt 6 -out hits-dna.tsv
blastp -query proteins.faa -subject reference.faa -outfmt 7 -out hits-protein.tsv
```

## Use the default tabular columns {#standard}

Both tabular formats use these 12 columns, in this order, unless you provide a field list.

| Order | Fields | Meaning |
| --- | --- | --- |
| 1–2 | `qseqid`, `sseqid` | Query and subject sequence IDs. The reader names these roles `qaccver` and `saccver` and preserves the values |
| 3–6 | `pident`, `length`, `mismatch`, `gapopen` | Percent identity, alignment columns, mismatches and gap openings |
| 7–10 | `qstart`, `qend`, `sstart`, `send` | One-based inclusive endpoints on each sequence. Descending coordinates preserve alignment direction |
| 11–12 | `evalue`, `bitscore` | Expected count and score in bits. Lower E-values and higher bit scores generally indicate stronger hits |

Extract the [input bundle](../assets/tutorial-inputs.zip) and work from the directory containing `examples/`. The example file preserves the first 12 fields of observed HSPs.

{{code:standardblast}}
{{output:standardblast}}

The reader puts the six plotting fields first in the returned table.

## Select six fields for plotting {#minimal}

{{code:blastendpoints}}
{{output:blastendpoints}}

A ribbon uses the two sequence IDs and their start and end positions. Keep `pident` to colour or filter by identity and `length` to filter by alignment length. For a custom headerless file, pass the actual field order as `col_names` to `read_blast()`. For a table already in R, see [link conversion](import-sequences.html#adapters).

## Choose the coordinate unit first {#protein}

| Objects being compared | Native unit | Alignment input | Next step |
| --- | --- | --- | --- |
| Nucleotide sequences | bp | BLASTn | Supply full FASTA lengths and a link table |
| Protein sequences | aa | BLASTp | Declare sequence_unit = "aa" with full protein lengths |
| Protein hits on a nucleotide genome | bp | BLASTp + CDS models | Project protein intervals to nucleotide coordinates |

One plot uses one coordinate unit. Full sequences supply protein lengths, while HSP endpoints describe only aligned intervals. BLASTx and tBLASTn have endpoints in different units and require a chosen common display domain. Nucleotide strand and protein N→C order have distinct biological meanings.

[Draw a native BLASTp protein comparison](protein-alignments.html#protein), or continue to [the four plotting modules](modes.html).

<span id="projection"></span>To locate protein hits on a nucleotide genome, read [CDS coordinate projection](projection.html).

## Sources and limits {#source}
Original EMALE records come from the public gggenomes EMALE dataset.
Only data records are reused; all figures are drawn with ggchord, without copying
upstream figures or plotting code. The input bundle's `sources.json` records URLs,
archive members, line ranges, hashes and selection. Parsing neither reruns
alignments nor establishes biological function.

## Check your result {#check-result}
Choose nucleotide or protein coordinates, then check table fields, full lengths and units.

<div class="next-step">Continue: <a href="projection.html">CDS coordinate projection →</a></div>
