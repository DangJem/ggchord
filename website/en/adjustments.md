<div class="task-intro"><p>Put common adjustments in the appropriate ggplot2 component: mappings, geometry, positions, coordinates and themes have distinct roles.</p></div>

<span id="prepare"></span>

## Assign each responsibility {#grammar}
[Plotting concepts](modes.html#prepare) explains layers, mappings, scales and coordinates. The example below applies that distinction to alignment identity.

## Map link identity {#mapping}

{{code:adjustments}}

{{figure:mappingexample}}

## Place features {#position}
Use `position_strand(c("+" = .1, "-" = -.2))` for signed strand offsets.
`position_feature_stack()` handles overlapping features; `position_track()`
allocates a named quantitative lane. Geometric sequence parameters accept one
value, a matching unnamed vector, or a vector named by `accver`; named vectors
are recommended. This syntax does not extend to scales, guides or themes.

Each layer owns its own `data`, `mapping` and parameters. Repeated feature or
link layers do not replace one another. Use ordinary ggplot2 `show.legend` and
`inherit.aes`; role input mappings precede geometry, while staged visual mappings
follow ggplot2 evaluation. See the [API contract](reference/ggchord_api_contract.html).

## Check your result {#check-result}
Use a coordinate for global rotation, named values for per-sequence geometry, and mappings plus scales for data-dependent colour.

<div class="next-step">Continue: <a href="parameters.html">Set values by sequence and strand →</a></div>
