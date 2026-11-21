<a id="inputs"></a>

## Find the table you need {#inputs}

| What you have | What to produce | Where to start |
| --- | --- | --- |
| Complete sequence lengths, feature intervals or two-ended records | Tables with matching accessions and known units | [Reference and layer tables](data.md) |
| FASTA, GenBank or GFF3 | A length table or annotated intervals | [Read sequence and annotation files](import-sequences.md) |
| BLAST, PAF or feature-group membership | Alignment endpoints or explicitly identified members | [Read relationships](import-relations.md) |
| bedGraph or BigWig | Interval coordinates and a measured value | [Read a signal](import-signals.md) |
| An amino acid interval and its parent CDS | The corresponding genomic bases and strand | [Project residues onto a CDS](projection.md) |

<span id="reader-results"></span>

<a id="results"></a>

## Follow a record into a figure {#results}

Register the complete reference with `ggchord()`, and pass observations to each layer explicitly. Use `ggchord_data()` to store related reference, feature and link tables together. Retain identifiers and original endpoints when filtering. Use `validate_ggchord_data()` to check input roles, and `data_report()` to inspect evidence attached by a reader or data operation. [Layout and source queries](diagnostics.md) connect a plotted layer to its input. [Preview and export](output.md) covers final figure dimensions.
