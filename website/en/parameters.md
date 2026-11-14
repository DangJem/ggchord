Flexible parameter formats refine shared settings into local control for parameters that support them. These examples use the same four bundled sequences; the first three compare only `seq_radius` on a 7 × 7 inch canvas. Requires ggchord and optional dplyr for CDS preparation.

{{code:courseblockac}}

## One value for every sequence {#scalar}

{{code:parameteruniform}}
{{figure:parameteruniform}}

Every radius is 1. A scalar is sufficient when all sequences share a setting.

## Change one sequence by name {#recommendation}

{{code:parameters}}
{{figure:parameterexample}}

Only `MT118296.1` uses 0.72; the rest retain this parameter's default of 1. Keys match actual `accver` values after display reordering. Keep complete named vectors when every sequence has a different setting.

## Values in display order {#vector}

<span id="compatibility"></span>

{{code:parametervector}}
{{figure:parametervector}}

Values match display order after `seq_order` and must cover all displayed sequences. The second sequence is currently `MT118296.1`, making this plot equivalent to the preceding one.

### Matching after reordering {#reorder}

Keep the vector and change `seq_order`. Named values still change `MT118296.1`; the unnamed vector changes the second displayed sequence, now `MT108731.1`.

{{code:parameterreorder}}

```r
view_ggchord(p_parameter_named_order, width = 7, height = 7, units = "in")
view_ggchord(p_parameter_ordered_values, width = 7, height = 7, units = "in")
```

![Named values after reordering](../assets/parameternamedorder.png)
![Ordered values after reordering](../assets/parameterorderedvalues.png)

## Prepare observed CDS on both strands {#genes}

{{code:parametergenes}}
{{output:parametergenes}}

`OR222515.1` has a plus-strand record at 2624–4195 and a minus-strand record at 19246–20523; the integrase on `OQ646790.1` is plus-strand. Strands come from `source_strand`, with original coordinates retained.

## Set strand positions {#strand}

{{code:parameterstrand}}
{{figure:parameterstrand}}

Explicit strand values are already signed: 0.08 for plus and −0.08 for minus, without another reversal. A scalar `position_strand(0.08)` expands to plus 0.08 and minus −0.08.

### Combine sequence and strand differences {#lists}

Place the two `OR222515.1` strands asymmetrically, while the other sequence uses a smaller symmetric offset. A list expresses these two levels of variation.

{{code:parameterlist}}

![Offsets by observed sequence and strand](../assets/parameterlist.png)

`OR222515.1` uses explicit 0.08/−0.14; the scalar 0.06 for `OQ646790.1` expands symmetrically. Omitted sequences or strands become 0 in this Position family. Do not copy this mixed list directly into `width`.

## Set strand widths {#width}

{{code:parameterwidth}}
{{figure:parameterwidth}}

Plus uses width 0.03 and minus 0.08. Width changes arrow thickness; start/end still determine interval length. Per-sequence width lists require elements following that family's strand-value rules; see [geometry parameters](reference/geom_feature.html).

## Place names by sequence {#labels}

{{code:parameterlabels}}
{{figure:parameterlabels}}

Text radii for `MT108731.1` and `OQ646790.1` are 75% and 115% of their backbone radii. Other sequences keep their defaults. Explicit `seq_label_radius` is incompatible with Collinear.

## Families and omitted values {#scope}

| Interface | Formats and matching | Boundaries and omitted values |
| --- | --- | --- |
| `geom_seq()` | Scalar, accver names, display order and corresponding scalar lists | `seq_radius`, `seq_curvature`, `seq_gap`, `seq_orientation` have individual defaults and Coord support |
| `geom_seq_label()` | Labels, radius, rotation, size and hjust/vjust by sequence | `seq_label_orientation` is a single mode; explicit radius is incompatible with Collinear |
| `geom_link_ribbon()` / `geom_link_line()` | `link_gap` by endpoint sequence | `NULL` retains automatic behavior; control points are not sequence scalars |
| `geom_gene()` / `geom_feature()` | Width by sequence or strand; list elements follow strand-value rules | Omitted entries retain the parameter default; arbitrary mixed Position lists do not apply |
| `position_strand()` | Scalar, sequence values, explicit strand values and mixed sequence lists | Scalar x becomes +x/−x; explicit strands are not reversed; omitted sequences/strands become 0 |
| `position_plasmid()` | Scalar or per-sequence scalar offsets | Both strands share the offset; strand-specific values are unsupported; omitted sequences become 0 |
| Fixed `geom_gene_label()` / `geom_feature_label()` | `label_radial_offset`, `label_circum_offset`, `label_circum_limit`, `label_rotation` by sequence or strand | Follow this family's strand rules; not every text parameter accepts lists |
| `position_track()` | One absolute `unit()` or per-sequence numeric width/gap | Numeric settings cover every displayed sequence; absolute sizes retain their physical meaning |

`seq_order` specifies ordering and `label_layout` specifies a layout mode. Automatic label layers do not accept arbitrary fixed offsets. Axes/themes, aesthetics, scales and guides retain their own interfaces. Manual colour scales match **data categories**; see [styling](styling.html). Consult the [parameter atlas](parameter-atlas.html) for exact defaults and ranges.

## Choose a format and check results {#check-result}

Use scalars for shared settings, named values for local differences, equal-length vectors for display-order settings, and supported lists when two levels of differences require them. Check IDs, order, omitted entries and actual results. Other sections use the simplest expression that serves their task; vector/list counts are not a quality measure.
