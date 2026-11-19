<a id="source"></a>

## Recover the input to one named layer {#source}

Give the ribbon layer a stable ID and request `data_stage = "source"`. The three resulting records show the accessions and original endpoints that entered this layer. To see both source and ribbon geometry, also request `include = "ribbon"`.

{{example:inspect-source}}

<a id="checks"></a>

## Choose the evidence for a symptom {#checks}

| Question | What to inspect |
| --- | --- |
| Where did a ribbon go? | Compare source IDs and endpoints, active sequence selection and adjacency rules. A repeated reference may need a named display instance. |
| Which intervals contributed to a coverage value? | Look up the track's `track_statistics` table and its `source_rows`. A [worked example](tracks-statistics.md#members) includes the actual denominator. |
| Why are labels hidden or crowded? | Run `diagnose_ggchord()` with the target device and dimensions. Inspect `label_diagnostics`, measurement status and the rendered figure. |
| Which panel contains this position? | When projections differ across facets, call `map_ggchord_position()` with the intended `panel`. |

When querying calculation stages, read `layer_stage_status` first. A successful empty statistic has status `available`. A failed stage has status `failed` and a reason. These cases require different interpretations.
