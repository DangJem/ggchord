## Unreleased RC checkout {#rc-development}
The 0.16.6.9000 source checkout audits the existing candidate contracts and retains
flexible parameter forms, aliases and reexports. It fixes inside primer callouts
being moved outside by the circular exterior allocator. Inside restriction registry rows are also corrected to match their inner geometry. Local finite checks and
benchmarks are distinct from release/platform certification; v0.16.5 remains the
formal install tag. No new migration or deprecation is introduced.

## Current version {#status}
These pages document **ggchord 0.16.5**, including cumulative changes from the 0.16.0–0.16.4 development phases.
Consult [GitHub Releases](https://github.com/DangJem/ggchord/releases) for formal tags, assets and publication dates.

## Changes in this version {#changes}
- Reorganized learning, complete workflows, task guides and reference entry points. The first-plot lesson covers installation, minimal inputs, interpretation, saving and troubleshooting; task pages prepare their own objects.
- Complete Han function reference, including arguments, returns and explanatory prose, with same-page English switching and source-hash checks.
- Redesigned typography, chapter navigation, gallery and data encodings. Choose overall chord, single-sequence circular or local CDS views by question. Persistent navigation, mobile outlines and reserved image dimensions support long-page reading.
- Native pkgdown / Bootstrap 5 package-documentation layout with same-page language switching, search, contents, code disclosure and exact code copying.
- Chord, circular and collinear home-page examples; every documentation figure supports enlargement and actual-size viewing.
- Bundled or attributed real inputs, with explicit row/column selection before function calls.
- GC tracks computed from actual DNA and CDS membership checked by sequence identity.
- Standard 12-column BLAST workflows, real PAF/BLASTp input, CDS projection and format adaptation.
- Sequence/strand parameter formats, supported scope, recommendations and compatibility; native aa plotting remains a separate proposal.
- Fixes for empty primer layers and lost strand/user colours when gene tables omit optional annotation text.
- Full-resolution user artwork with lossless RGBA-preserving encoding. Candidate signatures and default palettes are unchanged; documentation palettes are explicit in code.

## History and migration {#history}
The [complete NEWS source](../assets/NEWS.md) retains past changes.
Historical calls are not current examples; use the [migration guide](migration.html).
Local development checks do not constitute formal release validation.
