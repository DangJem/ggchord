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

- `nucmer.delta` and `promer.delta` are unchanged upstream MUMmer examples at
  tag v4.0.0rc1. They contain existing alignment coordinates, error counts and
  signed gap scripts, not source FASTA strings. The preamble's historical
  FASTA paths are provenance only; the reader never opens or retrieves them.
  The nucleotide example's original organism name is B_anthracis; the protein
  comparison uses Drosophila. This is file-format/coordinate validation, not
  experimental work, sequence design or a newly run comparison.
  `DELTA-TOOLS-LICENSE.md` retains the exact upstream Artistic 2.0 notice.
  No upstream implementation or executable is distributed in ggchord.
  `read_delta()` requires an existing user-selected show-coords executable.
  PROMER identity is an AA alignment statistic, while source endpoints and
  complete source lengths remain bp; these units must not be conflated.

After the 2026-11-12 sample constraint, `nucmer.delta` is **historical only**:
its original bytes, organism provenance, original source URL and past results
remain intact, but active checks and regeneration no longer load/download it.
The active real example is the unchanged Drosophila `promer.delta`; nucleotide
mode uses clearly marked synthetic schema/direction controls, not renamed
biological data. A new benign real NUCMER sample has not been validated here.
Historical `v0.27.0-delta-results.json` is retained; new runs write the separate
`v0.27.0-delta-active-results.json` without overwriting that historical evidence.
