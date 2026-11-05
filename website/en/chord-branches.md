Ribbons encode intervals; point links encode selected positions. This page compares branches and paths on explicitly synthetic geometry examples, not observed alignments or evidence of homology.

## Prepare an illustrative shared endpoint {#input}

{{code:comparebranchdata}}

Three illustrative sequences each span 1000 bp. Two links share query interval 200–600 bp and connect to different target intervals. Coordinates are one-based inclusive. Every figure uses 9 × 7 inches; endpoints stay fixed across branch settings.

## Independent ribbons and a shared query trunk {#branches}

{{code:comparebranchnone}}
{{figure:comparebranchnone}}

With default `link_branch = "none"`, the ribbons travel independently and may overlap near the query end.

{{code:comparebranchquery}}
{{figure:comparebranchquery}}

Only change to `"query"`: equal query endpoints can share a trunk. Overlapping intervals are not necessarily equal endpoints. Real data without suitable groups may show no change when branching is enabled. `"subject"` groups equal subject endpoints instead.

## Change the trunk fraction {#fraction}

{{code:comparebranchfraction}}
{{figure:comparebranchfraction}}

Relative to the query-branch figure, only `link_branch_fraction` changes from default 0.15 to 0.35. The shared trunk becomes longer relative to the shortest member. It is not identity, a probability or a newly measured alignment length.

| Parameter | Default | Effect and constraints |
| --- | --- | --- |
| `link_branch` | `"none"` | none/query/subject; selects the shared endpoint |
| `link_branch_fraction` | `0.15` | In (0, 0.5]; none does not create a shared trunk |
| `link_ctrl_point` | `NULL` | Shapes far ends when branching; the trunk determines its junction |
| `pairs` | `NULL` | qaccver/saccver table selecting existing pairs, not generating alignments |

## Replace intervals with points {#point-links}

Point links require explicit `qpos/spos`; an interval table is not a point table. Here choose interval midpoints. This loses interval length encoded by ribbon width, so explain the selection rule first.

{{code:comparepointdata}}
{{code:comparelinkcurve}}
{{figure:comparelinkcurve}}

Default `link_type = "curve"` connects the selected positions with curves.

{{code:comparelinkstraight}}
{{figure:comparelinkstraight}}

Only change to `"straight"`, keeping input, anchors and dimensions fixed. Add `arrow = grid::arrow()` separately for direction. It represents the query-to-subject drawing convention, not independent evidence of gene strand or transcription.

## Real data, modes and exercises {#adapt}

Preserve real query/subject endpoints. Do not edit unequal endpoints to manufacture shared trunks. See [link preparation](ribbons.html) for filtering and bundling; geometric branching does not perform those statistics.

Chord and collinear modes accept these link layers; single circular mode rejects ribbon/link. Linear adjacency filtering still applies. Check `adjacent_only` or explicit `pairs` before showing cross-row ribbons; see [cross-row links](collinear-routing.html).

Exercise: change only the trunk fraction to 0.25. Return to the default fraction and compare curved versus straight point links. With your own data, first decide whether interval lengths or individual positions matter, then consult [`geom_link_ribbon()`](reference/geom_link_ribbon.html) and [`geom_link_line()`](reference/geom_link_line.html).
