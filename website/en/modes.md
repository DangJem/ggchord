Run the [shared setup](index.html#quick-start) once. Each object below is a normal
composable ggplot; print it or use `ggsave()`.

## Chord: compare multiple sequences {#chord}
{{code:chord}}
{{figure:chord}}
`coord_chord()` controls global rotation and fitting. Sequence-local radius,
curvature, gap and orientation belong to `geom_seq()`. A NULL link control point
uses the historical origin-based chord route; it does not mean a universal
straight link. Ribbon avoidance defaults to `uniform`, line avoidance to `none`.

## Circular: annotate one sequence {#circular}
{{code:circular}}
{{figure:circular}}
This example explicitly hides the axis to leave the feature labels clear.
`coord_circular()` requires a single sequence. General features can carry unknown
strand. `geom_feature_plasmid()` is a convenience preset with its own defaults;
choose it explicitly when those defaults fit the task. Primer and restriction
layers are described in [recognition](recognition.html).

## Collinear: align local comparisons {#collinear}
{{code:collinear}}
{{figure:collinear}}
`anchor` aligns specified intervals; `orientation = "auto"` chooses supported
orientations using available evidence. Inspect orientation diagnostics rather
than assuming biological inversion. Ribbons connect adjacent displayed rows by
default; use `adjacent_only = FALSE` or explicit `pairs` to retain others.
Non-adjacent outlines default to dashed; explicit styles take priority. These
rules differ from chord mode. [Focus and segments](views.html) preserve source
coordinates for cropped or locally reversed views.
