## Scope and timing {#scope}
These are current candidate replacements, not a v1.0.0 compatibility promise.
Earlier 0.16.x versions are locally completed but unreleased. Old calls in this
table explain migration and must not be executed as current examples. Exact
first-removal tags are not guessed where historical evidence gives only an upper bound.

## Replacement table {#table}
{{migration}}

## Run the replacement {#replacement}
In a fresh session, run this preparation before the current replacements:

{{prepare:setup,docstyle,plasmiddata}}

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
