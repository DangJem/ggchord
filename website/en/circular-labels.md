Compare fixed and automatic labels on the same six source-backed pUC19 annotations. Placement, size, wrapping and device dimensions jointly affect readability; a hidden label is not a missing annotation.

## Prepare this page independently {#input}

Requires optional dplyr. Original names remain in finish_features; display text merely appends actual coordinates. Intervals and strands stay fixed. All layers consume the same finish_labels, with 9 × 8 inch figures.

{{code:comparefinishcircle}}


## Fixed labels: orientation, side and rotation {#fixed-labels}

The feature orientation follows each feature. Each subsequent plot changes only horizontal orientation, inside placement or 30-degree rotation. Return to the common baseline instead of accumulating settings. Default overlap is hide, so crowded text may be omitted.

{{code:comparefixedfeature}}
{{figure:comparefixedfeature}}
{{code:comparefixedhorizontal}}
{{figure:comparefixedhorizontal}}
{{code:comparefixedinside}}
{{figure:comparefixedinside}}
{{code:comparefixedrotate}}
{{figure:comparefixedrotate}}


## Fixed-label parameters and conflicts {#fixed-parameters}

| Parameter | Default | Meaning |
| --- | --- | --- |
| label_orientation | feature | feature/horizontal/radial/tangent |
| label_side | outside | outside/inside/auto |
| label_overlap | hide | hide/nudge/allow; allow can overlap |
| label_rotation | NULL | Extra degrees, resolving to 0 |
| label_radial_offset | 0.04 | Local-normal displacement, not bp |
| label_circum_offset | NULL | Resolves to 0 along the circumference |
| label_circum_limit | NULL | Resolves to TRUE for circumferential limits |
| label_wrap | NULL | Optional character wrapping width |

Do not pass fixed placement controls to automatic repel layers. Match the feature Position; prefer unique feature_id values for subsets/reordering and never reuse an ID for a different interval. In linear mode radial/tangent mean vertical/horizontal rather than circular letter placement.



## Compare four automatic layouts on identical text {#automatic-labels}

Default feature first tries internal, then adjacent/external placement. Radial, auto and callout provide alternative placement strategies. Layout changes presentation, not source annotations. Keep identical text, requested size and dimensions across all four.

{{code:comparefeaturefit}}
{{figure:comparefeaturefit}}
{{code:compareradialfit}}
{{figure:compareradialfit}}
{{code:compareautofit}}
{{figure:compareautofit}}
{{code:comparecalloutfit}}
{{figure:comparecalloutfit}}


## Wrapping, shrinking, hiding and leaders {#fitting}

Starting from feature, change only label_fit to shrink; then return to baseline and set only label_wrap=14. Shrink explicitly permits smaller text down to the default minimum 1.5. It does not shrink indefinitely to guarantee all names. Character wrapping width is not a fixed physical width.

{{code:compareshrinkfit}}
{{figure:compareshrinkfit}}
{{code:comparewrapfit}}
{{figure:comparewrapfit}}

From the 14-character wrap example, change only label_fit to ellipsis. Then return to default feature placement and disable only external fallback. Omission/hiding changes display, not source names or features.

{{code:compareellipsisfit}}
{{figure:compareellipsisfit}}
{{code:compareexternaloff}}
{{figure:compareexternaloff}}


## Important automatic-label defaults {#automatic-parameters}

| Parameter | Default | Choices and limits |
| --- | --- | --- |
| label_layout | feature | feature/radial/auto/callout; gene layouts differ |
| label_fit | wrap | wrap/none/ellipsis/auto/shrink; omission does not edit source text |
| label_max_lines | 2 | Inspect together with wrapping |
| label_min_size | 1.5 | Lower bound for shrinking |
| label_side | outside | outside/inside/auto |
| external | TRUE | FALSE prevents exterior fallback; linear mode rejects this constraint |
| max_overlaps | Inf | Does not guarantee zero collisions |
| label_segment_overlap | fade | Appearance of conflicting leaders |
| label_segment_overlap_alpha | 0.18 | Faded leader opacity |
| label_segment_linetype | auto | Automatic leader line type |
| label_segment_layer | below | below/above; text retains the label layer position |

Read exported label placement and fit diagnostics before changing dimensions, label subsets, layout or wrapping. Keep original full names rather than retaining only truncated figure text.




Feature layout uses direct leaders, allowing them to cross features without
multi-bend detours. Filled features cover leaders drawn below; transparent
fills keep their transparency. Internal single-line names follow circular
tracks. Wrapped sources keep one annotation and a continuous glyph on closed
circles; a visible coordinate gap retains separate display pieces. A shifted adjacent name receives a short attachment. To restore the
label layer's drawing order, use `geom_feature_label_repel(label_segment_layer = "above")`.

## Exercises and your annotations {#adapt}

Change only label_max_lines to 3 and explain any improvement. At fixed dimensions, compare exterior fallback on/off. Fixed and automatic layers have different controls. Supply explicit label text or anno, checking IDs, intervals, strands and unique feature_id values. Continue with [primer text](primers.html#primer-labels) and the [dense pETDuet-1 case](plasmids.html#complex).


Neighbouring names compare radial travel with angular displacement; an empty
inner track can preserve association better than squeezing names onto the
nearest track. Dotted arrow boundaries can come from source segments or
cleavage markers. Inspect `segments`, `segment_name` and `cleavage_arrows`;
a boundary is not necessarily an independent feature. Top-level feature
labels do not automatically add text for every named subsegment.
