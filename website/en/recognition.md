<div class="task-intro"><p>Search known DNA for primers, restriction sites and common features, then choose a display. Computational matches and experimental validation are different.</p></div>

## Prepare this page {#prepare}
This example uses bundled pUC19.

## Search and display are separate {#search}
Use the bundled pUC19 DNA for exact catalogue-primer and selected restriction searches:

Run `library(dplyr)` to select the sequence fields.

{{code:recognition}}
{{figure:recognition}}

Search supports circular-origin binding; set `circular = FALSE` for a linear sequence.
The actual pUC19 hits shown here should be inspected in the returned table. Search hits are computational matches, not experimental validation.
Primer `"three_prime"` mode finds an exact suffix, with no mismatch model;
its Wallace melting-temperature estimate is a composition heuristic.

## Restriction interpretation {#restriction}
`find_restriction_sites()` returns unfiltered matches and cleavage coordinates;
`filter_restriction_sites()` selects a display set. Enzyme site counts refer to
the full supplied sequence and remain stable after filtering. `methylation =
"auto"` uses audited context only for exactly recognized bundled references;
it does not claim unknown samples are unblocked. This example explicitly uses `"none"`; this is a search assumption, not a measurement
of methylation in a physical sample. `geom_restriction_site(leader = "radial")` is the default route.

## General feature recognition {#features}
`find_common_features()` supports curated exact-reference annotations, DNA and
record-specific matching modes. Protein searches use standard-code six-frame
translation; they do not implement splicing, alternative initiation or
origin-crossing protein matches. Source genetic-code metadata does not select
another translation table. For nonstandard codes use DNA matching or external
code-aware analysis. Confidence is a heuristic label, not a probability.
Use [explicit CDS projection](projection.html) for known protein/genome maps.

## Follow an annotated reference workflow {#reference}
This short workflow uses the same source-backed plasmid as the complete map. The [complete pUC19
workflow](plasmids.html) demonstrates catalogue annotations, exact reference
matching, site display filtering, primer labels and feature stacking together.
Changing only a sequence name does not make it a recognized reference; sequence
content determines source-backed recognition. Read the result methods and
methylation provenance instead of assigning reference claims to an unknown sample.

## Check your result {#check-result}
Inspect match positions, methods and provenance. Restriction filtering is a display choice; methylation search conditions are not a measurement of an unknown sample.

<div class="next-step">Continue: <a href="plasmids.html">A complete circular annotation workflow →</a></div>
