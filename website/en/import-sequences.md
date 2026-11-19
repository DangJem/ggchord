<a id="gff3"></a>

## Read annotated intervals from GFF3 {#gff3}

The installed GFF3 contains CDS annotations for MT108731.1. Keep positions, strand and the original line number together so that a plotted feature can be traced to its record. `unstranded = "keep"` preserves an unknown strand where one occurs.

{{example:installed-gff}}

<a id="reference"></a>

## Add a complete reference {#reference}

Use `read_fasta_lengths()` for a table of FASTA IDs and lengths. GFF3 annotations alone may end well before the end of a sequence. For this GFF3, select MT108731.1 from `seq_data_example` to obtain its complete reference length, then pass the GFF3 table to `geom_feature(data = features)`. When reading your own files, check the exact ID spelling in both tables.

A GenBank record contains both the sequence and its annotations. `read_genbank()` returns the reference in `seq_data` and the feature rows in `features$annotation`. Retain joined segments and their source strand when using a coding feature for [coordinate projection](projection.md).
