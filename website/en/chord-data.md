`ggchord()` uses the sequence table to establish the backbones. Add links and features as needed. The rows below come from the bundled data. Run the code to reproduce them.

## Sequence table `seq_data` {#tables}

{{code:frontseqtable}}

{{output:frontseqtable}}

| Required column | Meaning |
| --- | --- |
| `accver` | Unique sequence ID. You can use your own contig names. Other tables match on this value. |
| `length` | Full sequence length in bp. |

**Optional columns** can support display or grouping, such as a separate short display name.

## Link table `ribbon_data` {#inspect}

{{code:frontlinktable}}

{{output:frontlinktable}}

| Required column | Meaning |
| --- | --- |
| `qaccver`, `saccver` | IDs of the query and subject sequences. Both must occur in the sequence table. |
| `qstart`, `qend` | The two interval endpoints on the query sequence. |
| `sstart`, `send` | The two interval endpoints on the subject sequence. |

**Optional column** `pident` is percent identity and can be used for filtering or colour. The first row is about 98.6%. A link table can also have a `length` column for the alignment length.

## Genes and generic features {#names}

{{code:frontgenetable}}

{{output:frontgenetable}}

| Column | Purpose |
| --- | --- |
| `accver` | Required. ID of the sequence containing the interval. |
| `start`, `end` | Required. Interval start and end, with `start <= end`. |
| `strand` | Gene arrows require `+` or `-`. Use a generic feature layer for annotations without known direction. |
| `anno` | Optional. Can supply label text. |
| `type` | Optional. Category of a generic feature, such as CDS, promoter or replication origin. |

The bundled example stores the original annotated strand in `source_strand`. For your own imported GFF3 records, use their `strand` directly.

## Reading coordinates {#coordinates}

Input coordinates are 1-based and inclusive. For example, 10 through 12 covers three bases. Gene intervals carry direction in `strand`. Link intervals can encode direction in ascending or descending endpoint order, so preserve the order at each end after import.

Continue to [Plotting concepts](modes.html), or see [Input formats and conversion](formats.html) when reading files.

The BLAST-style fields reflect ggchord’s origins in alignment visualization. Query and subject use q/s prefixes, and custom IDs are accepted. See [alignment inputs](alignments.html) for actual BLAST rows.

## Protein sequence units {#protein-unit}

Nucleotide tables may omit the unit column and use bp. A complete protein length table declares `sequence_unit = "aa"`. A plot uses one unit. See [protein alignments](alignments.html#protein) for full inputs and figures.
