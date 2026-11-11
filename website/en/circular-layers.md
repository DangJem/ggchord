Adjust backbones, centre text and annotation placement on the same real pUC19. Establish the origin and direction before adding information layers; use one complete DNA sequence and source-backed annotations.

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

The zero origin is at 12 o’clock with clockwise coordinates. Positive rotation turns the complete plot counterclockwise. Relative to double above, each example changes only direction or gap. See the [rotation contrast](plasmid-tuning.html#rotation).

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


## Generic shapes and the plasmid preset {#feature-shapes}

On the same six source annotations, change only arrow to block. A block omits the directional arrow without changing strand. The third figure shows the combined geom_feature_plasmid preset: width=0.058, arrow_head_length=0.068, arrow_head_width=1.42 and plasmid Position. It changes several defaults and is not a single-parameter contrast. See [shape choices](annotations.html#shapes) and [arrow/short-feature geometry](feature-geometry.html).

{{code:comparefeaturearrow}}
{{figure:comparefeaturearrow}}
{{code:comparefeatureblock}}
{{figure:comparefeatureblock}}
{{code:comparefeaturepreset}}
{{figure:comparefeaturepreset}}
## Share placement between overlapping features and labels {#stacking}

Compare identity versus a shared feature stack on the same six annotations. Pass the Position to both geometry and labels rather than moving text alone. `position_feature_stack(spacing = 0.1, side = "strand", base_position = NULL)` applies an optional base Position before allocating overlap lanes. `position_plasmid()` is a local preset; `position_strand()` applies signed offsets by biological strand. Inside/outside is not a substitute for strand evidence.

{{code:comparefinishstack}}
{{code:comparepositionidentity}}
{{figure:comparepositionidentity}}
{{code:comparepositionstack}}
{{figure:comparepositionstack}}


## Exercises, replacement input and next steps {#adapt}

Change only double-line spacing to 0.04, then change the centre length_style size. Keep complete source lengths; do not infer them from maximum feature ends. For FASTA/GenBank/GFF3, check IDs and coordinates in [formats](formats.html). Continue with [circular labels](circular-labels.html), [primers](primers.html), [restriction sites](restriction-sites.html) and [tracks](tracks.html). Preview at the figure dimensions and use the same width/height in ggsave.
