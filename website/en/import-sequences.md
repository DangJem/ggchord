{{code:courseblockan}}


Examples use bundled tables or small installed files. Each block identifies the packages it uses.

<a id="adapters"></a>

## Adapt an existing table {#adapters}
The preparation above loads the input tables. These examples rename columns of
bundled source-backed sequences, links and genes; they do not invent alignments.

{{code:adapters}}

{{console:adapters}}

Adapters return ordinary frames with `ggchord_conversion` provenance. IDs remain
unchanged. Supply reliable full lengths rather than using maximum feature ends.
See [alignment input](alignments.md) for minimal fields and coordinate units.

<a id="genbank"></a>

## GenBank: retain record identity {#genbank}
{{code:genbank}}

This is the original NCBI L09137.2 record, named pUC19c; do not rename it to the
bundled pUC19. Results include `seq_data`, `features$annotation` and `report`,
retaining segmented locations and record/feature provenance. The selected Python
needs Biopython already installed; no automatic installation occurs.
`unsupported = "report"` records unsupported locations rather than guessing bounds.

This section requires a Python interpreter with Biopython. A successful read returns one sequence record, L09137.2, with a full length of 2,686 bp.

An explicit `python` path keeps the selected virtual environment, including interpreter symlinks; no environment is installed or modified.



<a id="files"></a>

## FASTA lengths and GFF3 annotations {#files}

{{code:frontimportseq}}
{{output:frontimportseq}}



FASTA supplies full lengths. GFF3 supplies intervals, types and source strands. Select CDS for plotting as needed.

{{code:frontimportgff}}

{{output:frontimportgff}}

<a id="blast"></a>

## Read a BLAST table {#blast}

{{code:frontimportblast}}

{{output:frontimportblast}}

Each row links two intervals. For a plot, check that both IDs occur in its backbone table. The pUC19 FASTA demonstrates length import. The BLAST example comes from a separate dataset. See [alignment inputs](alignments.md).

<a id="alternatives"></a>

## Convert Bioconductor objects directly {#alternatives}

With GenomicRanges installed, convert one interval into a feature table. Known full lengths are supplied in `seqlengths` and can be extracted with `as_ggchord_seq(gr)`.

{{code:frontimportgr}}

{{output:frontimportgr}}

The first GInteractions anchor is the query, and the second is the subject. This minimal example joins two intervals. It also requires InteractionSet.

{{code:frontimportgi}}

{{output:frontimportgi}}

See [import and conversion](import-intervals.md#ranges) for other inputs. Use the converted tables to [add layers](multi-sequence-start.md#draw).
