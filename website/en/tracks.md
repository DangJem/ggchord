## Declare position and numerical domain {#domain}
After the [shared setup](index.html#quick-start):

{{code:tracks}}
{{figure:tracks}}

The signal is a synthetic sine wave in arbitrary units, not measured coverage or
GC skew. ggchord renders supplied values; it does not compute sequence composition.
Replace it with your measured bins and declare their actual unit.

## Share an ID deliberately {#sharing}
`position_track("signal", ...)` chooses placement and physical width.
`scale_track_value_continuous("signal", ...)` defines the numerical domain;
it is not spatial placement. Layers sharing an ID train one joint domain,
including the bar baseline, and must declare consistent track settings.
Fill/colour scales are separately scoped to that ID. `geom_track_tile()` encodes
categories, `geom_track_line()` connects numeric observations, and
`geom_track_bar()` shows intervals relative to a baseline.

## Inspect retention {#diagnostics}
Use `export_ggchord_layout()` to inspect `track_registry` and `track_diagnostics`.
Their original value limits, physical mm contributions and native offsets have
different units. Focus tracks with `focus_ggchord_track()` alongside the same loci;
check clipped or omitted rows rather than treating an absent bar as zero.
