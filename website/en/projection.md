## Build an explicit coding map {#map}
This synthetic split codon spans two coding segments:

{{code:projection}}

`genomic$segments` preserves the actual coding pieces, not an envelope including
introns. `project_ggchord_links()` similarly projects explicit protein links;
it does not compute new homology or translate DNA.

## Phase, partial ends and exceptions {#phase}
CDS segments share `feature_id`, one sequence and a known strand. GFF3 `phase`
counts from each biological 5-prime end to the next codon boundary. Only the
first phase is skipped; later phases check accumulated coding length, retaining
codons split across segments. Explicit `segment_order` can declare an origin join.

Partial ends default to an error. `partial = "trim"` explicitly excludes bases
without complete codons; residue one means the first complete observed codon.
Genetic-code IDs are metadata for coordinate projection, not translation choices.
Declared frameshifts, slippage, translation exceptions and unsupported location
operators are rejected. No unobserved protein length or identity is inferred.

## Inverse mapping {#inverse}
Use `project_ggchord_features(genomic, cds, from = "genome")` to map back.
Non-codon boundaries require an explicit `codon_boundary` policy. Source
coordinates, parents and projection provenance remain attached. Read the full
[function reference](reference/project_ggchord_features.html) before applying
this to partial or multi-segment biological models.

## Project feature links as well as domains {#links}
The domain example above maps an observed residue interval through explicit CDS
segments. [Coding-coordinate links](homology.html#projection) adds a paired model
and demonstrates the separate ungapped assertion. Neither helper translates DNA,
validates genetic code or creates a new alignment. Compare the retained segment
pieces and source metadata before drawing discontinuous genomic anchors.
