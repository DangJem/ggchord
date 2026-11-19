A flexible parameter lets one layer give different geometry to individual sequences or strands. The **Flexible parameter** link beside an argument in the function reference leads back to the matching rules on this page.

<a id="scalar"></a>

## Change one backbone radius {#scalar}

<span id="recommendation"></span>

The four phage genome sequences in `seq_data_example` share a radius of 1 by default. Give `MT118296.1` a radius of 0.72 to bring its backbone inward. The name matches the sequence's `accver`. The other three retain radius 1.

{{example:flex-radius}}

| Value of `seq_radius` | Which backbone changes? |
| --- | --- |
| `0.72` | All four |
| `c("MT118296.1" = 0.72)` | Only MT118296.1 |
| `c(1, 0.72, 1, 1)` | The second displayed sequence |
| `list("MT118296.1" = 0.72)` | Only MT118296.1 |

<a id="reorder"></a>

## Match by name or by display order {#reorder}

<span id="vector"></span>
<span id="compatibility"></span>

Both panels display the same four sequences in the same reordered layout. On the left the smaller radius follows its name, `MT118296.1`. On the right it follows the **second position**, now occupied by `MT108731.1`. The sequence lengths, colours and sequence labels are identical across the panels.

{{example:flex-radius-order}}

An unnamed vector must contain one value per displayed sequence. For a radius list, each element holds one number. A named vector is the shortest way to change a few sequences and keep their settings when the plot is rearranged.

<a id="strand"></a>

## Set strand offsets and arrow widths {#strand}

<span id="genes"></span>
<span id="lists"></span>
<span id="width"></span>

Two CDS on `OR222515.1` lie on opposite source strands. Give the plus CDS an outward offset of 0.08 and the minus CDS an inward offset of −0.14. A single value of 0.06 gives the two strands of `OQ646790.1` symmetric offsets. The example uses dplyr to select the CDS and read each recorded `source_strand`.

{{example:flex-strand}}

An explicit minus-strand offset already includes its sign. Missing sequences or strands have offset 0 in `position_strand()`. Width changes arrow thickness. Here the plus strand uses 0.03 and the minus strand uses 0.08; the endpoints still come from `start` and `end`. When specifying widths by sequence **and** strand, supply a strand vector for each named entry, for example `list("OR222515.1" = c("+" = 0.03, "-" = 0.08))`.

<a id="invalid-radius"></a>

## What happens if a radius uses strand names? {#invalid-radius}

A backbone belongs to one sequence. Passing `+` and `-` as radius names asks for two sequence IDs absent from the input. This example catches the error when the plot layout is built. Use an `accver` for a backbone radius, or [position_strand()](../../man/position_strand.Rd) for strand offsets.

{{example:flex-radius-error}}

<a id="scope"></a>

## Find an argument by task {#scope}

<span id="labels"></span>
<span id="check-result"></span>

| Task | Arguments to look up | Matching rule |
| --- | --- | --- |
| Shape sequence backbones | [geom_seq()](../../man/geom_seq.Rd) · `seq_radius`, `seq_curvature` | Both default to 1. Use finite positive radii and finite signed curvatures. |
| Set direction and spacing | [layout_sequences()](../../man/layout_sequences.Rd) · `orientation`, `gap` | Each sequence has a direction of 1 or −1. Gap units and defaults depend on the view. |
| Place a sequence name | [geom_seq_label()](../../man/geom_seq_label.Rd) · `labels`, `seq_label_radius`, `seq_label_rotation`, `seq_label_hjust`, `seq_label_vjust` | Match the label to the accession. Radius is relative to its backbone. |
| Clear a connection endpoint | [geom_link_ribbon()](../../man/geom_link_ribbon.Rd) or [geom_link_line()](../../man/geom_link_line.Rd) · `link_gap` | Match a gap to the sequence at each end. |
| Set feature thickness | [geom_feature()](../../man/geom_feature.Rd) · `width` | Match by sequence, strand, or both. |
| Move features | [position_strand()](../../man/position_strand.Rd) or [position_offset()](../../man/position_offset.Rd) · `offset` | Use the first for differences between strands and the second for a shared displacement on each sequence. |
| Place fixed labels | [geom_feature_label()](../../man/geom_feature_label.Rd) or [geom_gene_label()](../../man/geom_gene_label.Rd) · `label_radial_offset`, `label_circum_offset`, `label_circum_limit`, `label_rotation` | Match by sequence or strand for fixed placement. Automatic labels choose their own positions. |
| Size a quantitative band | [track_band()](../../man/track_band.Rd) · `width`, `gap` | Supply one physical unit or native numeric values covering every displayed sequence. Width is positive and the gap may be zero. |

For data-driven colours use [aesthetic mappings and scales](styling.md). For the formals, defaults and flexible-parameter flags, see the [parameter index](parameter-atlas.md).
