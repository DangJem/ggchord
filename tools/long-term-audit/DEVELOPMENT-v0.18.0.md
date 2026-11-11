# v0.18.0 development record

2026-11-11: user approved local fast-forward of v0.17.0 f183d67 into main
and creation of v0.18.0. No push, deployment or release authorization.
Package development version: 0.18.0.9000.

## First implemented slice: HSP Dotplot

`coord_dotplot()` uses standard ribbon endpoints with independently selected,
ordered and reversed axes. `geom_alignment_segment()` is a thin segment geom,
not a midpoint projection of point links. Geometry is prepared independently
of the existing backbone-first pipeline. Native axis defaults reuse ggchord
axis styling; explicit ggplot2 theme/scale overrides survive preparation.
HSP rows stay source-local and per-layer, with four endpoints and metadata.
Exports include source tables, per-layer resolved inputs, axis transforms and
retained/omitted reasons. Explicit unknown direction is omitted with warning;
contradictory normalized direction is rejected. Declared repeated HSP IDs use
existing explicit authorization semantics. No inferred fine alignment trace.

| Combination | Current support |
| --- | --- |
| Alignment segment × Dotplot | Implemented; identity stat/position only |
| Native ggplot2 annotations × Dotplot | Available in concatenated native units |
| Sequence/ribbon/feature/track/custom label geoms × Dotplot | Explicit error; adaptation pending |
| Facets × Dotplot | Explicit error |
| Alignment segment × existing sequence coordinates | Explicit error |
| Existing Chord/Circular/Collinear/Sankey | Original geometry pipeline retained |
| Radial, its tracks/labels/routes | Not implemented in this slice |

## Direct evidence

Run `Rscript tools/long-term-audit/check-dotplot.R` from repository root.

- Real EMALE PAF: 11 HSPs; declared query/target lengths, reverse normalized
  endpoints, four-endpoint round trips and actual gtable passed.
- Real nucleotide BLAST: 7 HSPs; complete FASTA lengths, same checks passed.
- Real protein BLAST: 1 HSP; complete FASTA lengths and explicit aa, same checks passed.
- Explicit synthetic A/B/C: independent axis order/direction, two local layers,
  distinct metadata, self hit, declared repeated HSP authorization, explicit
  user colour scale/native axis theme, selection omissions, unknown/conflicting
  direction, source-coordinate focus round trip and empty layer build passed.
- Existing Chord sequence/ribbon representative build passed.
- Real PAF PNG inspected in /private/tmp, not added to Git. Default axes start
  at 1 and show true sequence endpoints; Identity guide has all six ticks.
- A single public Dotplot behavior case added to existing plot API tests for
  future authorized release acceptance. No testthat run in ordinary development.

This is a runnable first development slice, not completion of v0.18.0 or
formal release acceptance. Dotplot performance, new-view annotation/track
adaptation and spokes Radial remain open. No frozen prototype files, website
or repository PNGs were regenerated. No package build/check, vignette rebuild,
release device matrix, push, tag, deployment or release was performed.
