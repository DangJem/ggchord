## Shared terminology {#terms}
{{terms}}

## Design and asset provenance {#sources}
The grammar follows [ggplot2](https://ggplot2.tidyverse.org/): mappings, scales,
geometry, positions, coordinates and themes retain distinct responsibilities.
Compact query/subject field names follow BLAST-style tables; ggchord does not run
an alignment tool during drawing. Circular hierarchy and feature readability
are compared with Circos, clinker and DNA Features Viewer, as recorded in the
repository design documents. These are conceptual comparisons, not copied APIs.

The detailed ggchord artwork was supplied by the user on 2026-11-01 for this
website and documentation. Full-size PNG and WebP preserve every decoded RGBA
pixel. Only the small navigation thumbnail is resized. See the [image manifest](../assets/icon.json)
for source hashes and exact sizes. The compact hexagonal study mark below, site
CSS/JS and documentation diagrams are original project assets under the
repository MIT terms. The study mark is supplementary, not a replacement for
the supplied artwork.

![Compact hexagonal study mark](../assets/logo.svg) System fonts are used; no webfont, icon
library or third-party image is downloaded. Plot colours come from the package's
existing defaults or explicit teaching code; those literal colour choices are
listed with the source. The site uses the existing Pandoc and R documentation
converters at build time without bundling their binaries or stylesheets.
