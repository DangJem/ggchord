This composition follows an observed alignment from a four-sequence overview into its local context. Three overview panels share a subset of HSPs and CDS. The remaining panels show the two full-length references, a focused neighbourhood and a short row of observed bases.

<a id="map"></a>

## Draw the full map {#map}

The figure uses a 16 × 14 in canvas at 120 dpi. The script returns the plot so you can pass it to `ggsave()` at the same size.

{{example:case-fusion}}

<a id="reading"></a>

## Read the figure {#reading}

Start with A (chord), B (collinear) and C (radial). Cross-panel ribbons connect corresponding intervals. In D, the outlined windows locate the region expanded in E. The two `fusion_magnify_manual()` components connect D to E and E to F with dashed leaders to the detail frames. F shows the source bases as outlined squares, with source coordinates below the cells. The within-panel identity scale and the cross-panel correspondence legend describe different relationships.

<a id="inputs"></a>

## Source records and display choices {#inputs}

The overview uses one observed HSP per sequence pair, chosen closest to 4,000 bp, and three CDS per sequence. `flagship_fusion_example` includes the selected endpoints, focus windows, and the 60 source-positioned bases for F. The axis in E uses the original coordinates, even though only the windows are visible.

[Build a two-panel Fusion](fusion.md) · [Register a character panel](../../man/register_ggchord_alignment.Rd)
