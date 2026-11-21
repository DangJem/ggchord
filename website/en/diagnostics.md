<a id="source"></a>

## Recover the input to one named layer {#source}

Select the layer by its ID with `layers = "observed"` and request `stages = "source"`. The source table shows the three records that entered this layer. To compare those records with the display geometry, request `stages = c("source", "geometry")` and supply the target device and dimensions.

{{example:inspect-source}}

<a id="checks"></a>

## Choose the evidence for a symptom {#checks}

| Question | What to inspect |
| --- | --- |
| Where did a ribbon go? | Check source IDs, endpoints and window bounds. For a repeated reference, inspect or set the display instance on each endpoint. |
| Which intervals contributed to a coverage value? | Export the `statistic` stage and inspect `source_rows`. A [worked example](tracks-statistics.md#members) includes the actual denominator. |
| Why are labels hidden or crowded? | Run `diagnose_ggchord()` with the target device and dimensions. Inspect the returned `issues`, measurement status and the rendered figure. |
| Which panel contains this position? | When projections differ across facets, call `map_ggchord_position()` with the intended `panel`. |

Read `stage_status` alongside the returned tables. An available stage can contain zero rows. A failed or unavailable stage gives a reason. Check `measurement` before comparing positions measured at different sizes.
