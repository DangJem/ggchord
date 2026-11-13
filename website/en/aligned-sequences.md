{{code:courseblockaq}}

{{prepare:setup}}

Extract the [input bundle](../assets/tutorial-inputs.zip) so the working directory contains examples/. Each section names its optional dependencies.

## Retain actual FASTA characters {#characters}

Use `keep_sequence = TRUE` to retain characters with their physical header line and optional file provenance. Only whitespace is removed; case is retained. Declare bp or aa explicitly instead of guessing the alphabet. Focusing crops the string as well as the interval; local position 1 below is original position 1001. A decreasing locus selects the same interval and does not reverse-complement it. Lengths alone cannot supply characters; HSP endpoints cannot supply mismatch or gap columns. Gapped FASTA counts gap characters and needs a separate alignment-to-source map before genomic use.

{{code:courseblockar}}
## Existing alignments: keep columns and source positions separate {#aligned-characters}

These examples require the source checkout's `examples/alignment/` directory. It contains unchanged Biopython DNA/protein alignment inputs, the upstream permission notice, hashes and a reproducible FASTA re-encoding rule. No aligner is run. Full record IDs and original characters remain unchanged; choose rows and columns explicitly.

`read_aligned_fasta()` reads one existing alignment. Declare DNA, RNA or AA; equal lengths alone do not establish homology. `alignment_column` includes `-` and `.` gaps, whose `source_position` is NA. Default positions count the ungapped input fragment from one; they are not complete-genome coordinates. Lowercase letters consume positions, so this reader must not be used as an A3M insertion parser.

{{code:courseblockas}}

The first selected row has a gap at alignment column 343. Column lookup therefore returns `matched = TRUE`, `is_gap = TRUE` and a missing source position. Inverse lookup of source position 343 finds its actual alignment column, without interpolation. Positions removed by filtering return `matched = FALSE`; filtering never renumbers retained columns. `map_alignment_position()` selects an alignment-row ID even when a separately declared source ID differs.

`geom_residue()` is a standard Cartesian ggplot layer. Styles and user scales follow ggplot2. It measures glyphs against actual cell width and height at draw time: `overflow = "error"` requires readable letters; `"tile"` retains colour cells where text cannot fit, without shrinking fonts. Focus fewer columns or increase physical width to show letters. For example, save `p_characters` with ordinary `ggsave(width = 12, height = 3.4, units = "in")`; simply raising DPI adds no physical room. Discrete rows, identity/reverse position scales and Cartesian coordinates are supported; this layer does not directly project characters onto backbone Coords. A Cartesian character plot with explicit source mapping can be registered using `register_ggchord_alignment()` for [Fusion character details](fusion-details.html#character-detail).

### Protein excerpt {#aligned-protein}

The same data contract accepts actual protein characters. This excerpt contains upstream gaps and residue differences; fill distinguishes literal symbols, not inferred functional domains.

{{code:courseblockat}}

### Explicit source coordinates {#aligned-source}

Provide one `source_data` row per aligned record with `accver`, `start` and `strand`. `start` is the first nongap source coordinate. Negative-strand positions decrease; strings must already be in alignment orientation. An optional `source_accver` identifies the original molecule. If you provide a full forward `sequence`, every nongap character is checked, using nucleotide complementation on the negative strand and order reversal only for AA. Missing source strings remain unverified; no origins or strands are guessed.

The following example checks consistency against ungapped fragments derived from the same alignment. It is deliberately **not independent genomic validation**. For a real genomic projection, supply independently known source sequences and interval metadata instead.

{{code:courseblockau}}

No mismatch/gap detail can be recovered from six HSP endpoints alone. This workflow accepts existing columns, not CIGAR/A3M/Clustal/Stockholm parsing, phase inference, consensus/Logo analysis or alignment computation. Explicit CDS translation is described below; exact macro/detail registration is available in the Fusion lesson.
## Source-backed CDS translation {#cds-translation}

`build_ggchord_cds()` defines observed coding coordinates; `translate_ggchord_cds()` additionally needs actual complete forward DNA/RNA strings and optional **Biostrings**. This real EMALE negative-strand CDS has source phase 0 and an ATG start. Table 11 is an explicit choice checked against the unchanged upstream protein FASTA. Input provenance is in `examples/format-inputs/README.md`.

{{code:courseblockav}}

`models` reports coding/protein strings, genetic code, initiation/ambiguity policies, partial ends and terminal remainder. `nucleotides` retains each genomic base; `residues$source_positions` contains the three genomic positions of each codon in biological order, even across exon boundaries. Residue positions use the CDS feature ID's observed-protein namespace; they are not an inferred external protein accession. Terminal and internal stops remain `*`. The example removes only the final stop for comparison with the upstream FASTA convention.

The default `initiation = "literal"` does not recode alternative starts; `"start"` requires a complete 5-prime end. `ambiguous = "error"` rejects fuzzy codons; explicit `"solve"` delegates to Biostrings, using X if unresolved. Unknown codes, missing/short sources and declared translation exceptions fail. `partial = "trim"` belongs to the CDS map: only the initial phase and incomplete terminal codon are excluded; later segment phases validate continuity without deleting junction bases. No missing sequence is reconstructed. RNA must be declared explicitly. Ordered `segments` support the existing spliced/origin-join contract.

A complete runnable nucleotide/codon/macro example is `examples/m2-detail.R`; its negative-strand labels are centred on their three coding bases. Real spliced plus/minus controls, source hashes and original expected translations are retained in `examples/cds/`. [Register the character plot](fusion-details.html#character-detail) to connect its exact source boundaries.
