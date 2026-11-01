## Compute a track from actual DNA {#domain}
After the [shared setup](index.html#quick-start), compute GC fraction in 100 bp
windows of bundled pUC19. The final shorter window uses its actual base count.


{{code:tracks}}
{{figure:tracks}}

The domain is 0–1, in fractions. Replace it with the actual unit/range of other
measurements. `position_track()` controls placement/physical width;
`scale_track_value_continuous()` controls numerical mapping.

## Share an ID deliberately {#sharing}
Same-ID layers share a numerical domain, including bar baselines, and must agree
on width/placement. Colour scales are scoped independently by track ID. Tiles
show interval categories, lines connect observations, bars encode interval values
relative to a baseline. Sharing a domain does not require sharing a spatial band.

## Inspect retention and limits {#diagnostics}
The export's `track_registry` and `track_diagnostics` record domains, physical
sizes, native offsets and clipping. Omitted/clipped rows do not mean zero values.
Focus tracks with `focus_ggchord_track()` using the same windows.

## Combine bars, a moving mean and annotation bands {#multiple}
Run the [circular example](modes.html#circular) to obtain real `features` first.

{{code:multitrack}}
{{figure:multitrack}}

Bars and line use the same actual GC measurements; the line is a three-window
moving mean. Incomplete end windows are omitted, not filled with invented values.
Line positions use integer bin midpoints. Annotation bands have an independent ID
and palette. Bands are explicitly ordered outside; overly wide inner bands at
small radii can fold and are rejected.
