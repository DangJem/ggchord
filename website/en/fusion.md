## Independent views, connected source intervals {#task}

Fusion is the composition layer. This example places a Chord overview, a Collinear view, a Radial view and a local two-sequence detail into explicit panel rectangles. All links refer to observed HSPs and registered source intervals; repeated displays do not become additional observations.

{{manualfigure:flagfusion}}

## Panel composition {#path}

[Panels and source links](fusion-start.html) · [Position, size and guides](fusion-layout.html) · [Cross-panel routing](fusion-routing.html)

## Local magnification {#magnify}

Use explicit regions to connect a global view to independently constructed detail plots. Registered aligned characters can provide exact residue or nucleotide details. [Local magnification and character details](fusion-details.html).

## Nested and circular regions {#nested}

Nested Fusion preserves leaf-panel paths. An origin-crossing region needs explicit wrapping and compatible complete source lengths. [Nested composition and origin-crossing zoom](fusion-advanced.html).

## Supported panel types {#inputs}

Registered genomic views, nested compositions and registered Cartesian character details are supported according to their contracts. Dotplot, facets and unregistered external plots are not inferred or adapted automatically.
## Flagship display selections {#flagship-selection}

Each of the three flagship overview panels retains 32 original HSPs and 266 CDS. Only three actual observations are selected for cross-panel correspondences. The local view selects the shortest observed HSP of at least 2000 bp, with 800 bp of context on either side (clipped at sequence boundaries), showing existing CDS and names. Pale cross-panel bands mark corresponding intervals; within-panel ribbon colours show Identity. Display selections and repeated views are not additional observations.
