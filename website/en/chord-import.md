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

See [import and conversion](formats.html#ranges) for other inputs. Use the converted tables to [add layers](chord-build.html).
