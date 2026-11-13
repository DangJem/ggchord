Prepared annotations need only lengths and tables. DNA is used here to prepare the example annotation, not as a plotting requirement.

{{code:courseblockaj}}

## Input and provenance {#input}

Run independently in a fresh R session with optional dplyr. `finish_features` retains the first six recognition rows without changing intervals or strands. `finish_labels` appends actual coordinates to source names for label-fitting demonstrations; the original table remains separate.

{{code:comparefinishcircle}}
## Backbone styles and double-line spacing {#backbones}

In single circular mode, `seq_style = "auto"` selects double. Single draws one backbone and double draws two. Keep the same sequence and 9 × 8 inch size, changing only style or double-line spacing. `seq_backbone_gap` defaults to 0.025; band width `seq_backbone_width` defaults to 0.035. Neither is stroke `linewidth`.

{{code:comparebonedouble}}
{{figure:comparebonedouble}}
{{code:comparebonesingle}}
{{figure:comparebonesingle}}
{{code:comparebonegap}}
{{figure:comparebonegap}}
## Band backbone {#band-rendering}

Band represents the complete backbone with a filled polygon. seq_backbone_width controls its width, default 0.035. This is actual band geometry; only seq_style changes from the double baseline.

{{code:compareboneband}}
{{figure:compareboneband}}

From the band baseline, double only band width while retaining other settings.

{{code:comparebandwidth}}
{{figure:comparebandwidth}}
## Circular direction, rotation, opening and fit {#coordinates}

The zero origin is at 12 o’clock with clockwise coordinates. Positive rotation turns the complete plot counterclockwise. Relative to double above, each example changes only direction or gap. See the [rotation contrast](circular-layout.html#rotation).

| Parameter | Default | Type and boundary |
| --- | --- | --- |
| `rotation` | 0 | Finite degrees; 90 places origin at 9 o’clock |
| `direction` | clockwise | clockwise/counterclockwise; not biological strand |
| `gap` | 0 | Degrees in [0,360); does not remove bases |
| `fit` | labels | labels/geometry/manual; controls inclusion of labels in fitting |
| `ratio` | 1 | Global aspect ratio, not a change in units |
| `clip` | off | Coordinate clipping can cut exterior text |

Automatic fitting may change overall scale and the readable subset of tick labels at fixed device size. Inspect text and leaders before fixing xlim/ylim; clipping is not loss of source data. Circular mode owns gap/direction and rejects geom_seq seq_gap/seq_orientation.

{{code:comparecirclecounter}}
{{figure:comparecirclecounter}}
{{code:comparecirclegap}}
{{figure:comparecirclegap}}
## Style centre name and length separately {#center}

Centre text is separate from labels along the backbone. `name = TRUE` and `length = TRUE` are separated by a newline by default; name also accepts a display string. Empty name_style/length_style lists accept fixed colour, size, family, fontface and lineheight. First hide only length. Then return to baseline and change only name_style size, retaining the length text style.

{{code:comparecenter}}
{{figure:comparecenter}}
{{code:comparecenternolength}}
{{figure:comparecenternolength}}
{{code:comparecenterstyle}}
{{figure:comparecenterstyle}}

## Switch to a single backbone {#style}

{{prepare:firstcircle,firstfeatures}}

{{code:frontcirclestyle}}

{{figure:frontcirclestyle}}

`seq_style = "single"` replaces the double backbone with one line. The centre label and bp ticks remain visible.

## Rotate the map {#rotation}

{{code:firstfeatures,frontcirclerotation}}

{{figure:frontcirclerotation}}

The default `rotation = 0` places the sequence origin at 12 o’clock, with coordinates increasing clockwise. Here `rotation = 90` moves it to 9 o’clock, rotating the ticks and replication origin together. The annotation still covers the same source bp interval. Next, [prepare your own sequence and annotations](circular-data.html#tables).
