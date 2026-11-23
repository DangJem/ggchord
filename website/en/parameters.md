A flexible parameter assigns different settings to references or display instances. Use a named vector to match accessions. For physical placement, numeric values are in millimetres and absolute grid units are also accepted where documented. The **Flexible parameter** links in the API reference lead to the matching sections below.

<a id="named-label"></a>

## Move the label of one phage {#named-label}

Keep both backbone radii at 1. Use a fixed horizontal label layout with a named offset of 18 mm for MT108731.1 and select the inward side with `side = "negative"`. Only the named label moves inward.

{{example:flex-labels}}

<a id="scalar"></a>

## Change one backbone radius {#scalar}

<span id="recommendation"></span>

All four genomes share radius 1 by default. Set `radius` in `layout_sequences()` to bring MT118296.1 inward. Radius is a proportion of the standard radius.

{{example:flex-radius}}

| Setting for `radius` | References affected |
| --- | --- |
| `0.72` | All four |
| `c("MT118296.1" = 0.72)` | Only MT118296.1 |
| `c(1, 0.72, 1, 1)` | The second displayed sequence |
| `list(source = c("MT118296.1" = 0.72))` | Only MT118296.1 |

<a id="reorder"></a>

## Match by name or by display order {#reorder}

<span id="vector"></span>
<span id="compatibility"></span>

Both panels use the same reordered layout. On the left the small radius follows its accession, MT118296.1. On the right it follows the second displayed position, now occupied by MT108731.1. An unnamed vector must contain one value for every displayed instance. Use `list(source = ..., instance = ...)` when instance settings should override a shared source setting.

{{example:flex-radius-order}}

<a id="strand"></a>

## Set strand offsets and arrow thickness {#strand}

<span id="genes"></span>
<span id="lists"></span>
<span id="width"></span>

In the Chord example the two CDS from OR222515.1 use signed offsets of 5 mm and −9 mm for the plus and minus strands. A scalar of 4 mm for OQ646790.1 gives its two strands opposite offsets of the same magnitude. The script reads the original strand and uses separate layers for plus and minus arrow thicknesses of 2 and 5 mm.

{{example:flex-strand}}

<a id="invalid-radius"></a>

## What happens if a radius uses strand names? {#invalid-radius}

A backbone belongs to a reference. Using `+` and `-` as radius keys asks for accessions that are absent. The example captures the layout error. To move strand-specific features, use [position_strand()](../../man/position_strand.Rd).

{{example:flex-radius-error}}

<a id="scope"></a>

## Find an argument by task {#scope}

<span id="labels"></span>
<span id="check-result"></span>

| Task | Arguments to look up | Matching rule |
| --- | --- | --- |
| Shape arc references | [layout_sequences()](../../man/layout_sequences.Rd) · `radius`, `curvature`, `angular_gap` | Radius is a positive proportion, curvature is signed and angular gap is in degrees. |
| Adjust straight references | [layout_sequences()](../../man/layout_sequences.Rd) · `orientation`, `row_gap`, `column_gap` | Orientation is 1 or −1. Collinear and Sankey row gaps and Sankey column gaps use mm. |
| Place a fixed label | [layout_label_fixed()](../../man/layout_label_fixed.Rd) · `offset`, `along` | Specify perpendicular distance and displacement along the reference in mm. Use a named setting where the receiving layer supports it. |
| Move features or separate lanes | [position_strand()](../../man/position_strand.Rd), [position_offset()](../../man/position_offset.Rd) and [position_feature_stack()](../../man/position_feature_stack.Rd) | Use strand-aware offsets, shared offsets or lane spacing. |
| Set feature thickness | [geom_feature()](../../man/geom_feature.Rd) and [geom_gene()](../../man/geom_gene.Rd) · `thickness` | Use positive values in mm. Name values by reference or use an explicit instance setting. |
| Size a quantitative band | [layout_track()](../../man/layout_track.Rd) · `thickness`, `gap` | Use a positive thickness and a nonnegative gap. Both are physical sizes. |

For data-driven colours, use [aesthetic mappings and scales](styling.md). The [API reference](reference.md) gives the exact formals and defaults.
