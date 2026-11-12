# Source-backed alignments

`opuntia.aln` and `protein.aln` are unchanged Biopython example/test inputs.
The first contains Opuntia nucleotide records (AF191659.1 and related accessions);
the second contains protein records including CHSY_EQU / AAK54648.1 and related
chalcone-synthase records. Identifiers and alignment hypotheses are upstream facts,
not newly validated homology, complete genome lengths or inferred coordinates.

The `*-aligned.fasta` files concatenate each record's original CLUSTAL blocks,
preserving all characters/gaps/order. No aligner is run and no character is edited.
See `sources.json` for exact URLs, hashes, retrieval date and transformations.
Reproduce with `python3 data-raw/prepare-alignment-inputs.py`; input hash changes fail.

Biopython files are covered by the accompanying `BIOPYTHON-LICENSE.rst`;
retain that permission notice. No Biopython parser or plotting code is copied.
Tutorials choose explicit alignment columns and records without renumbering them.
Source positions default to the ungapped input fragment, starting at one, unless
users provide independently known `source_data`. They are not complete-genome coordinates.
