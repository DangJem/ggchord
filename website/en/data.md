<div class="task-intro"><p>Turn file observations into plotting roles. Check minimal fields and coordinates, then use validation reports to decide on cleaning; do not guess missing lengths or strands.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed. For file-based examples, first extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains `examples/`.

{{prepare:setup}}

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

## Reading coordinates {#reading-coordinates}

Input coordinates are 1-based and inclusive. For example, 10 through 12 covers three bases. Gene intervals carry direction in `strand`. Link intervals can encode direction in ascending or descending endpoint order, so preserve the order at each end after import.

Continue to [Plotting concepts](modes.html), or see [Input formats and conversion](formats.html) when reading files.

The BLAST-style fields reflect ggchord’s origins in alignment visualization. Query and subject use q/s prefixes, and custom IDs are accepted. See [alignment inputs](alignments.html) for actual BLAST rows.

## Protein sequence units {#protein-unit}

Nucleotide tables may omit the unit column and use bp. A complete protein length table declares `sequence_unit = "aa"`. A plot uses one unit. See [protein alignments](alignments.html#protein) for full inputs and figures.

## Coordinate contract {#coordinates}
Genomic coordinates are **one-based, inclusive bp**. BED/PAF readers convert their
source conventions; do not convert twice. Gene endpoints are ordered with
strand stored separately; link endpoint order may encode direction. A one-bp
link without direction evidence stays unknown. Generic features retain unknown
strand. The legacy `read_gff3()` default is an exception: `unstranded = "plus"`
converts it to `+` while retaining `source_strand`. Use explicit
`unstranded = "keep"` and `geom_feature()` to retain unknown strand.

Protein backbones use native aa coordinates; see [protein alignments](protein-alignments.html). Use [CDS projection](projection.html) to place protein intervals on nucleotide backbones.
Native Cartesian layout units, genomic bp, angles in degrees, physical mm and
track value units are different domains. See the full
[data contract](reference/ggchord_data_contract.html).

## Import original records {#import}
Run this repository-only example from its root. Original files are not rewritten.

{{code:imports}}

{{console:imports}}

`read_genbank()` retains record/feature provenance and segmented locations.
`read_bed()` and `read_paf()` retain their conversion metadata. `as_ggchord_seq()`,
`as_ggchord_feature()` and `as_ggchord_link()` adapt supported S4 objects; provide
full lengths and IDs explicitly when the object lacks them. `file` and `files` both accept path vectors or globs; repeated resolved paths
are read once in first-occurrence order. Choose one clear entry for a workflow. Missing sequence or
homology data must be supplied by the analyst, not invented by the plotting layer.

## Validate, clean, prepare {#validation}
After this page’s preparation code:

{{code:validation}}

{{console:validation}}

Validation reports issues; cleaning makes explicit repairs with a report.
Preparation combines workflow steps without inferring biological membership.
Inspect `validation`, `cleaned$report` and `prepared$report` before plotting.
`source_row` in geometry refers to the resolved layer input, while `.source_line`
means a physical file line; GenBank has separate record/feature provenance.
Aggregates may have `source_rows` and `.length_estimated` flags. Preserve these
when filtering or merging links.

## Inspect a problematic input and its repair report {#repairs}
{{code:cleaning}}
{{results:cleaning}}

The demonstration intentionally includes an unknown ID and
identity outside 0–100 in copied source rows. These explicit repairs are not a recommendation to clip
scientific observations silently. Inspect `problem_report` and `repaired$report`,
retain the original table, and choose policies appropriate to actual evidence.
A legitimate reverse alignment can keep endpoint order instead of being sorted.
`prepare_ggchord_data()` therefore retains its workflow-specific default policy;
standalone cleaning and preparation should not be assumed identical.

## Declare genome and contig keys {#hierarchy}
{{code:hierarchy}}

{{console:hierarchy}}

`hierarchy$contigs` summarizes declared membership without concatenating contigs.
No membership is inferred from an accession prefix. Provide globally unique IDs
when separate samples reuse local contig names. Use [file formats and adapters](formats.html)
for format/optional-object details, and the [real-record manual](real-data.html)
to carry reports into a complete figure.

## Minimal input and protein alignments {#minimal}
{{results:minimal}}

The second positional constructor argument is `ribbon_data`, not `gene_data`.
Use explicit argument names as in the shared examples. Optional identity, length
and coverage fields are needed only by operations that consume them. See
[standard BLAST/protein input](alignments.html) and [parameter formats](parameters.html).

## Result checks {#check-result}
Input IDs must match complete sequence records. Keep original records separate from repaired copies. Decreasing endpoints can be legitimate reverse hits; warnings do not all mean rows must be removed.

<div class="next-step">Continue: <a href="formats.html">File formats and explicit adapters →</a></div>

## Source tables and display copies {#source}

This chapter uses bundled sequence lengths, observed HSPs and source-strand
annotations. `feature_data` is a common annotation source. Generic feature layers
inherit the whole table; default gene layers require explicit gene/CDS classes.
An explicit layer table or data function takes priority. Use disjoint subsets
when combining the two kinds of geometry.

{{code:foundationdata}}

Prepare, clean, focus and segment helpers also accept `feature_data` and return
that member. Do not combine it with `gene_data`. Unknown strands remain unknown;
the legacy gene input retains its directional semantics and positional arguments.
