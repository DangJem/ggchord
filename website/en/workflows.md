<a id="tools"></a>

## Choose by input and goal {#tools}

Identify the original input and the intended plot. Every path preserves source IDs, complete lengths and coordinate units. Readers do not run alignments or generate missing relationships or annotations.

| Input available | Shortest processing path | Result and plotting entry |
| --- | --- | --- |
| Existing R tables | Organize [minimal fields](data.md#tables), check endpoint/annotation IDs and inspect [validation reports](data.md#validation) | Length, HSP and annotation tables → [shared data and layers](multi-sequence-start.md) |
| FASTA, GFF3 and BLAST/PAF | [Import sequences and annotations](import-sequences.md), then [alignments](alignments.md); check lengths, units and IDs, retaining reverse endpoints | Source tables → [relationship plot](multi-sequence-start.md) or [HSP Dotplot](dotplot-start.md#dotplot) |
| One GenBank record | Use `read_genbank()` in [sequence import](import-sequences.md); inspect lengths, annotation hierarchy and original segments | One sequence and annotations → [circular backbone and annotations](circular-start.md) |
| Protein BLASTp and complete lengths | [Declare aa units and read alignments](protein-alignments.md); verify both protein IDs; supply actual CDS/phase before projecting to DNA | Protein positions → [protein comparison](protein-alignments.md); DNA positions → [projection](projection.md) |
| BED/BEDPE intervals | [Read intervals or paired intervals](import-intervals.md) and check the reader's one-based conversion without converting twice | Annotation intervals or paired links → [features](annotations.md) and [links](multi-sequence-links.md) |
| bedGraph/BigWig signals | [Import signals](import-signals.md), retaining units, missingness and source intervals without filling unknown regions with zero | Signals → [tracks](tracks.md); counts/bins → [track statistics](tracks-statistics.md) |
| Existing character alignment | [Import aligned characters](aligned-sequences.md), checking gaps, source positions and alphabet | Observed alignment columns → [magnification and character details](fusion-details.md) |

<a id="checks"></a>

## Checks before plotting {#checks}

1. Sequence `accver` values are unique and lengths come from complete sources, never the last annotation endpoint.
2. Link `qaccver`/`saccver` and annotation `accver` match the sequence table. Explicitly distinguish repeated contig names across samples; readers do not guess membership.
3. DNA uses bp and proteins use aa; BED/PAF readers already convert coordinates. Ascending/descending link endpoints can encode direction and must not be uniformly sorted.
4. Draw a minimal plot before adding layers, [local settings](parameters.md) and [target-size output](output.md). Use [symptom-based diagnostics](diagnostics.md#symptoms) for problems.

<a id="examples"></a>

## Complete source-backed example {#examples}

[Sequence comparison from FASTA/GFF3/BLAST](real-data.md) retains paths, display filters, inputs and results. Obtain originals from the [input archive](formats.md), extract it and run downloaded scripts from a directory containing `examples/`.

## Choose the operation by its output {#choose-operation}

| Desired result | Entry point and interpretation |
| --- | --- |
| Local tables, several windows, or an overview with details | `focus_ggchord_data()` clips to a local interval. `segment_ggchord_data()` prepares windows, which can represent one source row more than once when they overlap. `magnify_ggchord()` returns an overview and detail composition. A clipped HSP length estimates the displayed span. See [local views](views.md). |
| Move a feature or its text | `position_feature_stack()` allocates lanes to overlapping features. See [label layout](labels.md) for placing an individual label and for automatic placement of a group. |
| Select or summarize links | Filtering keeps qualifying rows, deduplication selects representatives, interval merging produces new boundaries, and bundling records bin members. In plots, `after_stat(density)` is normalized by the largest bin weight within the sequence pair and expresses relative density. See [ribbon statistics](ribbons.md#statistics). |
