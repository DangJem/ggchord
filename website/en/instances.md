## Preparation {#input}

{{prepare:foundationdata}}

## Named A–B–A routes {#instances}

The three displayed rows refer to two biological sequences. `seq_order` values
identify sources; names identify instances. Filter each layer's original HSPs
before routing. `between` accepts two names or two declaration-order indices.
Ambiguous routes fail. The two subsets below are disjoint and preserve every HSP.

{{code:foundationinstances}}

{{figure:foundationvertical}}

The same `foundation_plot` can use `coord_chord()` or `coord_collinear()`.
Vertical Collinear changes the global view while keeping source strands and
positions. Positive local normal points left. Explicit text styles retain
priority. `map_ggchord_position()` accepts an instance name; its inverse input
is final y in vertical views.

`export_ggchord_layout()` supplies `display_instances`, `instance_routing` and
`instance_omissions`. Explicit `source_hsp_id` detects repeated display across
layers; use `allow_repeated=TRUE` only when that repetition is intended. Without
an ID, identity stays layer-local; equal row values do not establish identity.
