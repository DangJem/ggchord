<div class="task-intro"><p>Draw real annotations as independent layers, choose suitable shapes and share stacking positions between overlapping intervals and text.</p></div>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed.

{{prepare:setup,docstyle,plasmiddata}}

## Preserve feature semantics {#stack}
Run this page’s preparation code below. and [circular example](modes.html#circular).
This table retains actual overlapping pUC19 annotations and unknown strands. Geometry and fill map separately
from the same category; an explicit fixed `shape` would override shape mapping.

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

A region overlay uses its own table and `region_fill` scale. `geom_seq_label()`
names sequences, `geom_seq_center_label()` supplies name/length in a single circle,
and `geom_seq_break()` marks cropped/segmented ends. `geom_seq_scale()` is an
independent bp reference bar for collinear plots; it differs from
`scale_seq_position_continuous()`, which controls genomic axes on sequences.
An ordinary ggplot2 annotation should use [position mapping](views.html#mapping)
instead of treating bp as Cartesian x/y.

## Check your result {#check-result}
Check intervals and strands before shapes and labels. Unknown strands remain in general features; stacking changes display placement, not genomic coordinates.

<div class="next-step">Continue: <a href="labels.html">Labels and local placement →</a></div>
