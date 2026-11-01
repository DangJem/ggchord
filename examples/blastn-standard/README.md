# Standard tabular BLAST input

These 12-column files preserve exactly the first 12 fields and all HSP rows of
`examples/blastn/*.o7` (BLASTN 2.16.0+). They are column-format conversions,
not newly computed alignments. Rebuild with `Rscript data-raw/prepare-standard-blast.R`;
`provenance.csv` records source/output hashes and row counts.

A new run can use `blastn -query query.fna -subject subject.fna -outfmt 6 -out result.tsv`.
`-outfmt 7` uses the same default columns with comments. Omitting `-outfmt`
completely produces pairwise text; that is a different format.

Columns: query ID, subject ID, pident, length, mismatch, gapopen, qstart, qend,
sstart, send, evalue, bitscore. ggchord preserves identifier values and calls
their roles `qaccver` / `saccver`; an identifier need not be an accession.
For drawing, only the two identifiers and four endpoints are required.
