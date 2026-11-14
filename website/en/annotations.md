<div class="task-intro"><p>Draw real annotations as independent layers, choose suitable shapes and share stacking positions between overlapping intervals and text.</p></div>

<span id="prepare"></span>

## Preserve feature semantics {#stack}
Use the first six pUC19 annotations to see how overlapping intervals share stacked positions. Source strand is retained.

{{code:featurestack}}
{{figure:featurestack}}

`position_feature_stack()` resolves overlapping annotation lanes. The label layer
uses the same Position, so text anchors follow those lanes. Prefer unique
`feature_id` values for subset/reordered label inputs. The strict gene shortcut
requires known `+/-`; general features retain unknown/bidirectional semantics.

## Choose shapes by meaning {#shapes}
| Geometry | Typical use | Important distinction |
| --- | --- | --- |
| `arrow`, `compact_arrow` | Directional coding/feature interval | Known direction is evidence, not a display guess |
| `block`, `capped_line` | Region or unfilled interval | Shape choice does not infer biological function |
| `promoter_arrow`, `primer_arrow`, `primer_arc` | Declared promoter/primer-style annotation | Computational matching and drawing are separate |
| `marker`, `chevron`, `lollipop` | Custom interval/point-like emphasis | Keep real interval metadata; display shape is a mapping |

`geom_gene()` and `geom_feature_plasmid()` are useful presets, not aliases with
identical defaults. Use `position_strand()` for signed offsets and
`position_plasmid()` as a circular placement preset; the latter is not a fourth
coordinate mode.

## Highlight an interval without rewriting observations {#decoration}
{{code:decorations}}
{{figure:decorations}}

The highlighted interval is the observed pUC19 replication origin. A region overlay uses its own table and `region_fill` scale. `geom_seq_label()`
names sequences, `geom_seq_center_label()` supplies name/length in a single circle,
and `geom_seq_break()` marks cropped/segmented ends. `geom_seq_scale()` is an
independent sequence-position reference bar for collinear plots; it differs from
`scale_seq_position_continuous()`, which controls genomic axes on sequences.
An ordinary ggplot2 annotation should use [position mapping](views.html#mapping)
instead of treating bp as Cartesian x/y.

## Result checks {#check-result}
Check intervals and strands before shapes and labels. Unknown strands remain in general features; stacking changes display placement, not genomic coordinates.

<div class="next-step">Continue: <a href="feature-geometry.html">Gene and feature arrow geometry →</a></div>

Compare shouldered/flush/triangle arrows on one observed CDS in [arrow geometry](feature-geometry.html#arrow-heads).

## Region side, width and offset {#region-parameters}

Keep the original query ID and interval of the first observed HSP. Defaults are width 0.08, offset 0 and inside. Each subsequent figure changes only outside placement, width or offset. These are local geometry quantities, not bp. Set fixed fill in geom parameters or map categories to region_fill with its dedicated scale. Highlighting does not establish a new functional annotation.

{{code:comparefinishchord}}
{{code:compareregioninside}}
{{figure:compareregioninside}}
{{code:compareregionoutside}}
{{figure:compareregionoutside}}
{{code:compareregionwidth}}
{{figure:compareregionwidth}}
{{code:compareregionoffset}}
{{figure:compareregionoffset}}
