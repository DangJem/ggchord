# Existing genomic alignment sources

These two unchanged Biopython 1.85 inputs retain their original identifiers,
coordinates, characters and scores. Exact URLs and SHA256 are in `sources.json`.
The permission notice is retained in `BIOPYTHON-LICENSE.rst`;
no upstream parser implementation is copied. Verify or explicitly retrieve the
pinned originals with `python3 data-raw/prepare-genome-alignment-inputs.py`.

- `ucsc_mm9_chr10.maf` contains existing multiway genomic alignment blocks,
  mouse chromosome 10 as the first component and other declared source IDs,
  source lengths, gaps, positive/negative strands, qualities/context/empty
  component records. Full lengths do not independently verify characters.
- `dna_rna.chain` contains source-backed transcript/genome chains as well as
  explicitly `_modified` format controls from upstream. Scientific examples
  select source block 1 (NR_046654.1 versus chr3), without claiming the modified
  controls are new biological observations. The original file is never edited.
  Chain encodes alignment blocks/gaps, not nucleotide mismatch characters.

No aligner, assembly, lift-over or homology detection is run by the readers.
