Start with sequences and analysis results, then read them with ggchord. Missing comparisons must not be invented to populate a figure.

## Start with FASTA and GFF3 {#files}
Extract the [tutorial bundle](../assets/tutorial-inputs.zip) so the working directory contains `examples/`. FASTA supplies full lengths; GFF3 supplies intervals and strands. Retain source records separately from display tables.

{{prepare:sequencepalette,docstyle}}
{{code:showcasedata}}

`show_alignments` retains 32 observed HSPs and `show_links` selects 12. `show_genes` retains 266 CDS; `show_features` selects CDS of at least 500 bp from four groups across annotation order per sequence, additionally retaining source-backed negative-strand endpoints. Selection neither reruns alignment nor establishes functional homology.

For coordinate inspection, reshape the four endpoint columns with tidyr while retaining the original link table:

{{code:inspectends}}

## Produce your own BLAST results {#blast}
These producer commands are for analysts with BLAST installed and their own input files; the website builder does not execute them.

```sh
blastn -query query.fna -subject subject.fna -outfmt 6 -out result.tsv
```

```r
library(ggchord)
blast_source <- read_blast(file = "result.tsv")
```

Without a custom field list, outfmt 6 uses the standard 12-column table; 7 adds comments. Omitting outfmt entirely requests another format. Query/subject IDs map to qaccver/saccver. Consult the [NCBI field documentation](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/table/appendices.T.options_common_to_all_blast/) and [nucleotide/protein guide](alignments.html). Historical tutorial HSPs are converted records, not newly run BLAST results.

## When another format already exists {#alternatives}
| Source | Entry point | Check |
| --- | --- | --- |
| Existing table | as_ggchord_seq/feature/link | Explicit role mapping; retain the original |
| BED, PAF | read_bed, read_paf | Convert zero-based half-open coordinates only once |
| GenBank | read_genbank | Biopython required; preserve record and segmented-location provenance |
| GRanges, Seqinfo | as_ggchord_feature, as_ggchord_seq | Full lengths from Seqinfo; unknown strand stays unknown |
| GInteractions | as_ggchord_link, as_ggchord_seq | First anchor is query; full length context is still required |

These adapters address the request in [issue #2](https://github.com/DangJem/ggchord/issues/2). Current code provides the entries; the historical discussion does not certify every current combination. Missing optional packages are explicitly skipped locally. See [format conversion](formats.html#ranges) for runnable optional-object examples and limitations.

Continue to [1.3 Build from backbones](chord-build.html).
