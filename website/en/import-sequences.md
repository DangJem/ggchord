{{code:courseblockan}}

{{prepare:setup}}

Extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains examples/. Each section names its optional dependencies.

## Adapt an existing table {#adapters}
The preparation above loads the input tables. These examples rename columns of
bundled source-backed sequences, links and genes; they do not invent alignments.

{{code:adapters}}

{{console:adapters}}

Adapters return ordinary frames with `ggchord_conversion` provenance. IDs remain
unchanged. Supply reliable full lengths rather than using maximum feature ends.
See [alignment input](alignments.html) for minimal fields and coordinate units.

## GenBank: retain record identity {#genbank}
{{code:genbank}}
{{console:genbank}}

This is the original NCBI L09137.2 record, named pUC19c; do not rename it to the
bundled pUC19. Results include `seq_data`, `features$annotation` and `report`,
retaining segmented locations and record/feature provenance. The selected Python
needs Biopython already installed; no automatic installation occurs.
`unsupported = "report"` records unsupported locations rather than guessing bounds.

The GenBank example was rerun with Python Biopython 1.85 enabled: it returned L09137.2, length 2686 bp and one annotation row. For a separate documentation environment, install Biopython with that Python’s `python3 -m pip install biopython`, then put its executable on PATH before starting R. Do not treat an optional-dependency skip or NULL as successful import.

An explicit `python` path keeps the selected virtual environment, including interpreter symlinks; no environment is installed or modified.

Extract the [tutorial bundle](../assets/tutorial-inputs.zip) so the working directory contains `examples/`. Read one format at a time, inspect a few rows, then match sequence IDs.

## FASTA lengths and GFF3 annotations {#files}

{{code:frontimportseq}}

{{output:frontimportseq}}

FASTA supplies full lengths. GFF3 supplies intervals, types and source strands. Select CDS for plotting as needed.

{{code:frontimportgff}}

{{output:frontimportgff}}

## Read a BLAST table {#blast}

{{code:frontimportblast}}

{{output:frontimportblast}}

Each row contains two intervals. Keep their endpoint order and check both sequence IDs against the FASTA records. See [alignment inputs](alignments.html) for standard BLAST output.

## Convert Bioconductor objects directly {#alternatives}

With GenomicRanges installed, convert one interval into a feature table. Known full lengths are supplied in `seqlengths` and can be extracted with `as_ggchord_seq(gr)`.

{{code:frontimportgr}}

{{output:frontimportgr}}

The first GInteractions anchor is the query, and the second is the subject. This minimal example joins two intervals. It also requires InteractionSet.

{{code:frontimportgi}}

{{output:frontimportgi}}

See [import and conversion](import-intervals.html#ranges) for other inputs. Use the converted tables to [add layers](multi-sequence-start.html#draw).
