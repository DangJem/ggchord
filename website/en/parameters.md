## The capability remains supported {#recommendation}
Sequence-local geometry often needs one value per sequence or strand; duplicating
whole layers would obscure that intent. Keep scalars as the default, use vectors
named by `accver` when values differ, and name strand pairs `+`/`-`.

{{code:parameters}}

Run the [shared setup](index.html#quick-start) first. Named radii survive sequence
reordering; unnamed entries follow current positions. These are layout/geometry
values, not a universal extension to colours or themes.

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

## Compatibility and design tradeoffs {#compatibility}
Equal-length unnamed vectors, sequence-named lists, position-named lists and
unnamed lists remain compatible. Gene/feature entries can additionally contain
strand pairs. Each parameter still validates its own type, length, range and units.
This does not make a list a valid container for any argument.

Avoid expanding low-frequency equivalent forms: positional lists can silently
follow a new order after filtering; nested lists are harder to read. This version
does not remove them. The planned v0.16.6 audit will evaluate convergence;
any breaking decision needs explicit migration and removal timing. This design
recommendation is not an announced deprecation.

Numeric values do not all share one unit: `seq_gap` is a proportion, rotations
are degrees, and physical track sizes use explicit `unit()` values. `seq_backbone_gap`
and `seq_backbone_width` are scalar style controls, not members of the list rule.
Automatic collinear orientation belongs to `coord_collinear(orientation = "auto")`;
`geom_seq(seq_orientation = ...)` supplies numerical +1/-1 sequence directions.
