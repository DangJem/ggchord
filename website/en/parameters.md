Use one value for every sequence, name individual sequences, or configure the two strands. Each figure uses the same four bundled sequences.

## Apply one value to every sequence {#scalar}

{{code:parameteruniform}}
{{figure:parameteruniform}}

All four backbone gaps use 0.08. A scalar applies to every sequence.

## Change one sequence by name {#recommendation}

{{code:parameters}}
{{figure:parameterexample}}

The named vector sets MT118296.1 to radius 0.72 and leaves the others at 1. Values follow IDs after reordering.

## Set values in display order {#vector}

<span id="compatibility"></span>

{{code:parametervector}}
{{figure:parametervector}}

The four radii follow the four sequence-table rows. An unnamed vector must match the number of sequences. Named values are easier to check when changing order.

## Prepare two observed CDS on opposite strands {#genes}

{{code:parametergenes}}
{{output:parametergenes}}

The integrase CDS is on the plus strand of OQ646790.1 and Int is on the minus strand of OR222515.1. Preserve these source strands for the next examples.

## Set strand positions independently {#strand}

{{code:parameterstrand}}
{{figure:parameterstrand}}

The plus-strand offset is 0.08 and the minus-strand offset is −0.08. `position_strand()` changes placement while retaining source positions and strands.

## Set strand widths independently {#width}

{{code:parameterwidth}}
{{figure:parameterwidth}}

Plus-strand arrows use width 0.03 and minus-strand arrows use 0.08. Width controls visual thickness; start/end still determine interval length.

## Place labels by sequence {#labels}

{{code:parameterlabels}}
{{figure:parameterlabels}}

MT108731.1 uses a text radius at 75% of its backbone radius and OQ646790.1 uses 115%. The other two keep their default positions. Each supported parameter accepts its own named values.

## Supported parameters {#scope}
| Interface | Parameters | Matching scope |
| --- | --- | --- |
| `geom_seq()` | `seq_radius`, `seq_curvature`, `seq_gap`, `seq_orientation` | Sequence; radius/curvature effects depend on coordinates |
| `geom_seq_label()` | `labels`, `seq_label_radius`, `seq_label_rotation`, `size`, `seq_label_hjust`, `seq_label_vjust` | Sequence; `seq_label_orientation` is one mode choice |
| `geom_link_ribbon()` / `geom_link_line()` | `link_gap` | Endpoint sequence; Cartesian control points are not per-sequence scalars |
| `geom_gene()` / `geom_feature()` | `width` | Sequence and/or `+`/`-` strand |
| `position_strand()` / `position_plasmid()` | `offset` | Sequence/strand; prefer Positions over compatibility offset inputs |
| Fixed `geom_gene_label()` / `geom_feature_label()` | `label_radial_offset`, `label_circum_offset`, `label_circum_limit`, `label_rotation` | Sequence/strand, not every text setting |
| `position_track()` | Numeric `width`, `gap` | Native numeric placement by sequence; physical `unit()` is a scalar absolute size |

`seq_order` is itself an ordering. Automatic `label_layout` is a mode; manual angles
and offsets belong to fixed label layers rather than automatic-label `...`.
Axes/themes, aesthetics, scales and guides do not inherit the list extension.
Manual scales do accept vectors named by **data categories**, under ordinary
ggplot2 scale semantics; see [styling](styling.html).

## Check your result {#check-result}
Named-vector keys must match actual accver values or +/− strands. Recheck matching after reordering; list rules do not extend to themes, scales or aesthetics.

<div class="next-step">Continue: <a href="styling.html">Mappings, scales, guides and themes →</a></div>
