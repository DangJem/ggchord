## Search and display are separate {#search}
A tiny synthetic sequence demonstrates exact primer and restriction searches:

{{code:recognition}}

The primer crosses the circular origin. Set `circular = FALSE` for a linear
sequence. Search hits are computational matches, not experimental validation.
Primer `"three_prime"` mode finds an exact suffix, with no mismatch model;
its Wallace melting-temperature estimate is a composition heuristic.

## Restriction interpretation {#restriction}
`find_restriction_sites()` returns unfiltered matches and cleavage coordinates;
`filter_restriction_sites()` selects a display set. Enzyme site counts refer to
the full supplied sequence and remain stable after filtering. `methylation =
"auto"` uses audited context only for exactly recognized bundled references;
it does not claim unknown samples are unblocked. This synthetic example explicitly
uses `"none"`. `geom_restriction_site(leader = "radial")` is the default route.

## General feature recognition {#features}
`find_common_features()` supports curated exact-reference annotations, DNA and
record-specific matching modes. Protein searches use standard-code six-frame
translation; they do not implement splicing, alternative initiation or
origin-crossing protein matches. Source genetic-code metadata does not select
another translation table. For nonstandard codes use DNA matching or external
code-aware analysis. Confidence is a heuristic label, not a probability.
Use [explicit CDS projection](projection.html) for known protein/genome maps.
