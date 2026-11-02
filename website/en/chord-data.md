Learn the tables before selecting layers. Chord comparison was ggchord's original core view. Its early BLAST workflows explain the compact field names retained today. A multi-sequence chord view displays existing interval comparisons among sequences; drawing does not run a multiple-sequence alignment algorithm.

## Minimal shapes of three tables {#tables}
| Table | Required columns | One row represents |
| --- | --- | --- |
| seq_data | accver, length | One sequence and its full length in bp |
| ribbon_data | qaccver, saccver, qstart, qend, sstart, send | One observed pair of aligned intervals |
| gene_data | accver, start, end, strand | One directed gene interval |

Backbones alone need neither links nor genes. Use `geom_feature()` for annotations with unknown strand rather than assigning a direction.

## Names and optional columns {#names}
`accver` originated from accession.version, but any stable contig identifier is allowed; versions are not added or removed. `q` means query and `s` means subject. Both IDs must occur in the sequence table.

`pident` is percentage identity, needed only for identity mapping/filtering. Alignment `length` differs from full sequence `length`. Keep `evalue`, `bitscore`, `qlen` and `slen` when required by an analysis. Gene `anno` is optional annotation text; supply it for labels. Additional category columns can feed `aes()` without becoming mandatory inputs.

## Inspect real tables {#inspect}
```r
library(ggchord)
library(dplyr)
seq_data_example |> select(accver, length) |> slice_head(n = 3)
ribbon_data_example |> select(qaccver, saccver, qstart, qend, sstart, send) |> slice_head(n = 3)
gene_data_example |> transmute(accver, start, end, strand = source_strand, anno) |> slice_head(n = 3)
```

These are bundled source-backed records. This historical teaching fixture stores biological strand separately in `source_strand`; restore it explicitly. Newly read GFF3 annotations use their ordinary `strand` column.

## Coordinates and directions {#coordinates}
Use one-based inclusive bp coordinates. Genes have ordered bounds with a separate strand. Endpoint order can represent reverse alignments; do not sort every link. Maximum feature endpoints do not establish full sequence lengths. See the [coordinate contract](data.html#coordinates).

Continue to [1.2 From analysis results to input tables](chord-import.html).
