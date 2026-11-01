## Scope and timing {#scope}
These are current candidate replacements, not a v1.0.0 compatibility promise.
The 0.16.0–0.16.4 phases were delivered together in v0.16.5. The current
0.16.6.9000 checkout is unreleased; the formal installation tag remains v0.16.5.
The RC audit retains flexible parameter formats, spelling aliases, reexports and
getter/export entry points without a new deprecation. Old calls in this
table explain migration and must not be executed as current examples. Exact
first-removal tags are not guessed where historical evidence gives only an upper bound.

## Replacement table {#table}
{{migration}}

## Run the replacement {#replacement}
In a fresh session, run this preparation before the current replacements:

{{prepare:setup,plasmiddata}}

{{code:migration}}

Old minus-strand offset magnitudes need a negative sign in explicit Position
values. Mapping still uses `aes(feature_shape = type)` with a shape scale;
the fixed parameter is `shape`. Automatic labels use `label_layout`; manual
geometry uses the corresponding fixed-label layer. `label_layout = "aligned"`
becomes `"auto"`. Removed names are errors, not silently retained execution paths.

## What remains supported {#retained}
The five colour/color spelling aliases, all upstream reexports, gene/feature/plasmid
shortcuts and both geometry getter/export interfaces remain. Coordinate-specific
defaults and already documented result fields are retained. See the
[API contract](reference/ggchord_api_contract.html) and current reference for
NULL/NA/missing semantics; NULL does not universally mean disabled.
