<a id="task"></a>

## Independent views, connected source intervals {#task}

Fusion composes explicit panel rectangles containing a Chord overview, Collinear view, Radial view, two-sequence context, local zoom and observed nucleotide detail. Connections use existing HSPs and registered source intervals; displaying a record twice creates no new observation.

{{manualfigure:flagfusion}}

<a id="path"></a>

## Panel composition {#path}

[Panels and source links](fusion-start.md) · [Position, size and guides](fusion-layout.md) · [Cross-panel routing](fusion-routing.md)

<a id="magnify"></a>

## Local magnification {#magnify}

Use explicit regions to connect a global view to independently constructed detail plots. Registered aligned characters can provide exact residue or nucleotide details. [Local magnification and character details](fusion-details.md).

<a id="nested"></a>

## Nested and circular regions {#nested}

Nested Fusion preserves leaf-panel paths. An origin-crossing region needs explicit wrapping and compatible complete source lengths. [Nested composition and origin-crossing zoom](fusion-advanced.md).

<a id="inputs"></a>

## Supported panel types {#inputs}

Registered genomic views, nested compositions and registered Cartesian character details are supported according to their contracts. Dotplot, facets and unregistered external plots are not inferred or adapted automatically.
<a id="flagship-selection"></a>

## Flagship display selections {#flagship-selection}

A, B and C share six selected real HSPs and twelve source CDS; complete 32-HSP/266-CDS inputs remain available. One HSP closest to 4000 bp is selected per directed sequence pair. CDS are selected from three positional bins per sequence, preferring informative names and longer features. A-B, A-C and C-D correspondences use selected original records, with C connected to two panels.

D restores forward full-length sequences. E shows the corresponding HSP and source windows with 800 bp flanking context on either side, selecting two CDS per sequence. The same four CDS selected by source feature identifiers in E are retained in D; D no longer uses an independent selection inside these windows. Paired highlights in D and E are registered by `magnify_ggchord()`; actual magnification is remeasured for the device and panel allocation. Highlights are not additional HSPs and do not change source coordinates or biological strands. The flagship uses `connector = "envelope"`: each sequence in D has its exact source highlight individually outlined, without a padded rectangle enclosing multiple selections. Two straight leaders connect actual source-band edges to E. The E frame encloses sequences, features, ticks and labels as a display boundary only; its heading and description remain outside and close to the sequence content. Full-size and homepage outputs allocate panel space separately, enlarge A, widen B and use more compact homepage legends. The C-D connection avoids E. Cross-panel colours identify correspondences; internal ribbon colours show Identity.

F sits below C and E and shows the first 60 observed coding bases of the second CDS on E’s lower sequence: OQ646790.1, cds-WGH15551.1, forward strand, source positions 11006–11065 bp. Characters are extracted from the complete FASTA and verified against its source string; decorative leaders connect the exact short interval in E to F. Cells meet continuously with dark bold characters and chemistry colours by default. Original HSPs, all four local CDS and full genomic coordinates remain intact.

D and E use the same original source windows by default. E labels original genomic positions: 7122–11192 for MT118296.1 and 10615–14691 for OQ646790.1. Horizontal source-frame extents match these intervals; frame height remains independent. Identity guides in both homepage and full figures preserve numeric interval widths, so 0–50 occupies half the colour bar rather than making unequal intervals appear equally wide.

<a id="shared-guides"></a>

## Shared child legends {#shared-guides}

`ggchord_fusion(center, guides = "collect")` collects native outside child legends at the composition edge, displaying visually identical legends once. Default `guides = "keep"` retains the earlier layout. Original plots are not modified, and explicitly hidden legends remain hidden. Each panel still uses its own scales: shared semantics require matching limits, breaks and guide styles. Legends with the same title but different contents are not forcibly merged.

The flagship places Sequence, Strand, Identity and cross-panel Correspondence legends together at the bottom, freeing space in the children. Inside legends cannot be collected; nested cross-panel legends remain managed by their own compositions. Collected legends are measured and wrapped at the actual device size. If one legend is too wide, configure its native guide or enlarge the canvas.
