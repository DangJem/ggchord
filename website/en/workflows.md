## Choose by input and goal {#tools}

Identify the original input and the intended plot. Every path preserves source IDs, complete lengths and coordinate units. Readers do not run alignments or generate missing relationships or annotations.

| Input available | Shortest processing path | Result and plotting entry |
| --- | --- | --- |
| Existing R tables | Organize [minimal fields](data.html#tables), check endpoint/annotation IDs and inspect [validation reports](data.html#validation) | Length, HSP and annotation tables → [shared data and layers](multi-sequence-start.html) |
| FASTA, GFF3 and BLAST/PAF | [Import sequences and annotations](import-sequences.html), then [alignments](alignments.html); check lengths, units and IDs, retaining reverse endpoints | Source tables → [relationship plot](multi-sequence-start.html) or [HSP Dotplot](dotplot-start.html#dotplot) |
| One GenBank record | Use `read_genbank()` in [sequence import](import-sequences.html); inspect lengths, annotation hierarchy and original segments | One sequence and annotations → [circular backbone and annotations](circular-start.html) |
| Protein BLASTp and complete lengths | [Declare aa units and read alignments](protein-alignments.html); verify both protein IDs; supply actual CDS/phase before projecting to DNA | Protein positions → [protein comparison](protein-alignments.html); DNA positions → [projection](projection.html) |
| BED/BEDPE intervals | [Read intervals or paired intervals](import-intervals.html) and check the reader's one-based conversion without converting twice | Annotation intervals or paired links → [features](annotations.html) and [links](multi-sequence-links.html) |
| bedGraph/BigWig signals | [Import signals](import-signals.html), retaining units, missingness and source intervals without filling unknown regions with zero | Signals → [tracks](tracks.html); counts/bins → [track statistics](tracks-statistics.html) |
| Existing character alignment | [Import aligned characters](aligned-sequences.html), checking gaps, source positions and alphabet | Observed alignment columns → [magnification and character details](fusion-details.html) |

## Checks before plotting {#checks}

1. Sequence `accver` values are unique and lengths come from complete sources, never the last annotation endpoint.
2. Link `qaccver`/`saccver` and annotation `accver` match the sequence table. Explicitly distinguish repeated contig names across samples; readers do not guess membership.
3. DNA uses bp and proteins use aa; BED/PAF readers already convert coordinates. Ascending/descending link endpoints can encode direction and must not be uniformly sorted.
4. Draw a minimal plot before adding layers, [local settings](parameters.html) and [target-size output](output.html). Use [symptom-based diagnostics](diagnostics.html#symptoms) for problems.

## Complete source-backed example {#examples}

[Sequence comparison from FASTA/GFF3/BLAST](real-data.html) retains paths, display filters, inputs and results. Obtain originals from the [input archive](formats.html), extract it and run downloaded scripts from a directory containing `examples/`.
