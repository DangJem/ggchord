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

Start a plot from its complete reference table. Supply observations to their individual layers, or assemble a prepared object with `prepare_ggchord_data()`. Retain the accession and original endpoints through any display filtering. To check a plotted layer or a statistical window against its source observations, use [layout and source queries](diagnostics.md). For the dimensions of the final image, see [preview and export](output.md).
