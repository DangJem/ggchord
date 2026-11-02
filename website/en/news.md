## Current version {#status}
These pages document **ggchord 0.16.6**, including cumulative changes from the 0.16.0–0.16.4 development phases.
Consult [GitHub Releases](https://github.com/DangJem/ggchord/releases) for formal tags, assets and publication dates.

## RC fixes and progressive tutorials {#rc-development}
- Progressive chord → collinear → plasmid modules, with schemas, import, layer construction and 13 source-backed plasmid cases.
- Refined four-sequence showcases, unified backgrounds and readable arrow legends; explicit guide text sizes are respected.
- Corrected inside primer placement and restriction registry side; tibble optional columns no longer cause spurious warnings.
- Retained candidate signatures, flexible forms, aliases and reexports; no new migration or deprecation. Dense labels can still be hidden or overlap across layers.

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
