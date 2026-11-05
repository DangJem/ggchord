Searching cuts and choosing a display subset are separate steps. Compare leaders, text content, order and side on the same three pUC19 hits, without editing DNA or cut coordinates.

## Retain cut coordinates and search assumptions {#input}

Run independently with explicit EcoRI, BamHI and HindIII and methylation="none". This is a search assumption, not a measurement of sample methylation. enzyme_site_count refers to the complete input sequence and survives display filtering.

{{code:comparefinishdnabase}}
{{code:comparefinishsites}}
{{code:comparefinishsiterows}}
{{output:comparefinishsiterows}}


## Compare leaders at the same cut sites {#restriction-labels}

Default radial uses text bounds and circular order, with radial stubs and fan connections where needed. Elbow uses one bend; straight connects the tick tip to the final label edge. Only leader changes across the three figures; cut sites and label grouping are fixed.

{{code:comparerestrictionradial}}
{{figure:comparerestrictionradial}}
{{code:comparerestrictionelbow}}
{{figure:comparerestrictionelbow}}
{{code:comparerestrictionstraight}}
{{figure:comparerestrictionstraight}}


## Change content, order and side separately {#text-content}

Each figure returns to the radial baseline, changing only label_style, label_order or label_side. Default bold enzyme names indicate one cut in the complete target; explicit fontface replaces that automatic distinction. Positions come from computed cleavage, not motif midpoints.

{{code:comparerestrictionenzyme}}
{{figure:comparerestrictionenzyme}}
{{code:comparerestrictionposition}}
{{figure:comparerestrictionposition}}
{{code:comparerestrictionorder}}
{{figure:comparerestrictionorder}}
{{code:comparerestrictioninside}}
{{figure:comparerestrictioninside}}


## Key parameters and filtering meaning {#parameters}

| Parameter | Default | Effect |
| --- | --- | --- |
| label_style | enzyme_position | enzyme_position/enzyme/position |
| label_order | auto | auto/enzyme_position/position_enzyme; automatic mirrored ordering |
| label_side | outside | outside/inside/auto |
| leader | radial | radial/elbow/straight; never changes cuts |
| min_label_gap | NULL | Device-derived default; explicit genomic-fraction spacing |
| tick_length | 0.04 | Local-normal tick length |
| label_offset | 0.18 | Local-normal text offset |
| linewidth / label_size | 0.28 / 3.2 | Stroke width versus text size |
| fontface / family | NULL / NULL | Explicit fontface overrides unique-cutter bolding |

find_restriction_sites returns matches and cuts; filter_restriction_sites selects a display set. unique, six_plus, unique_6plus and reference have different conditions. Reference uses saved sets only for exact reference targets. methylation=auto applies audited backgrounds only to exact built-in references; state assumptions for unknown samples.



## Exercises, dense maps and limitations {#adapt}

Increase only min_label_gap and inspect outward label/leader movement; then restore defaults and set label=FALSE for ticks only. Retain full search results separately from the display subset. Increase export dimensions instead of deleting source observations to fit a small image. Recheck crowding after combining features and primers in [dense pETDuet-1](plasmids.html#complex). Computed cuts do not validate a digestion experiment.
