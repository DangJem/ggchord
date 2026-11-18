{{code:courseblockbg}}

{{prepare:trackdata,tracks,foundationdata}}

<a id="tracks"></a>

## Count observations or display supplied values {#tracks}

`bar` counts observations at source positions; `col` displays supplied values.
Each accession/window computes its Stat and data Position before view projection
or display replication. Layers share the default band when `track` is omitted. Add `track_band("depth")` before or after the layers and refer to it with `track = "depth"`. Use `position` for source-data adjustments such as stack or dodge.

{{code:foundationcount}}

{{figure:foundationcount}}

The shortcut below counts the same observations using the same physical band.

{{code:courseblockbh}}


| Geometry | Default Stat | Input meaning |
| --- | --- | --- |
| bar | count | x/pos observations, optional weight |
| col | identity | supplied y/track_value; interval or position |
| point, line, path, step, text | identity | source x and local y; text requires label |
| tile, rect, area, ribbon, segment, linerange, errorbar | identity | corresponding native ggplot2 extents |

Supported statistics include identity/count/bin/summary, interval coverage and the [documented estimation combinations](tracks-statistics.md#estimates). Bin requires explicit binwidth or breaks. A custom ggplot2 Stat can read the layer’s source indices from source_row and return the members of each result in a source_rows list column. Omitted membership is exported as unknown.

Source position uses x orientation; native text size currently uses mm and polygon outlines are full. Other orientations/size units/partial outlines fail explicitly. Native text styles and path arrows are forwarded.

Tracks use ggchord coordinates. Sample facets with fixed scales partition observations by facet columns in the source table. Samples compute statistics independently and share the value domain of each named band. Use `coord_collinear(seq_axis = "vertical")` for vertical projection. Source membership is recorded before data positions such as dodge.

<a id="interval-coverage"></a>

## Measure interval coverage in each window {#interval-coverage}

`stat_coverage(data = intervals, binwidth = 4)` draws the fraction of each window covered by at least one interval. Map `aes(y = after_stat(mean_depth))` for mean depth. Coordinates start at 1 and include both endpoints. An interval that crosses a window boundary contributes its intersection to each window. Mean depth counts the length of each record’s intersection separately.

Consider a synthetic sequence of length 10. Rows 1, 2 and 3 span 3–6, 4–5 and 10–10, respectively. Positions 4 and 5 are each covered twice. The last window contains just two positions, so its denominator is 2.

| Window | Distinct positions covered | Coverage fraction | Mean depth | Source rows |
| --- | --- | --- | --- | --- |
| 1–4 | 2 | 2/4 | 3/4 | 1, 2 |
| 5–8 | 2 | 2/4 | 3/4 | 1, 2 |
| 9–10 | 1 | 1/2 | 1/2 | 3 |

Empty windows in observed groups return zero coverage, zero depth and an empty source list. References or samples without observations produce no statistical rows. Each reference or focused window starts binning at local position 1. Alternatively, pass explicit intervals in `windows`. Each independently intersects the observations and windows may overlap. Choose either `binwidth` or explicit windows. Interval and window units must agree with the reference. By default out-of-bounds intervals raise an error. Use `bounds = "clip"` to clip intervals to the reference. Split compound intervals and origin-crossing records into single segments first.

<a id="independent-bands"></a>

## Two independent bands and display order {#independent-bands}

Encode the same GC by bar height and colour using explicit gc_height/gc_colour IDs. Set height and colour domains 0–1 separately, with distinct physical widths. The second figure changes only colour order from 2 to 0, moving that band nearer the backbone. Order is a finite priority, not a measurement or required one-based index; ties use IDs. Layers sharing an ID use its declared side, order, width and gap. Lines do not implicitly bridge independent windows.

{{code:comparemultitrackdata}}
{{code:comparetracktwo}}
{{figure:comparetracktwo}}
{{code:comparetrackorder}}
{{figure:comparetrackorder}}

<a id="window-policy"></a>

## DNA windows and line boundaries {#window-policy}

GC here is `mean(bases %in% c("G", "C"))`: ambiguous bases remain in the denominator and contribute zero to the numerator. The final window uses its actual number of bases. This pUC19 input has no ambiguous bases; for other inputs explicitly choose and report a policy rather than assuming they are omitted. Focusing clips intervals and retains the original summary; it does not recompute GC for a shortened window. Line tracks connect the supplied samples within their grouping; `closed` explicitly controls circular closure and `gap` belongs to physical track placement. Do not use a closed line to imply continuity across omitted source windows.

<a id="estimates"></a>

## Native source-space estimates {#estimates}

`geom_track()` also delegates `smooth`, `density` and `ecdf` to native ggplot2 statistics **before projection**. These create estimates of supplied observations, not additional observations. They run independently by accession/window and explicit group. Estimation `source_rows` records the eligible group submitted to the estimator after native missing-value handling; it is not per-point influence, positive read support or a custom model's additional internal filtering report.

| Statistic | Documented geometry | Input and explicit choices | Data Position |
| --- | --- | --- | --- |
| `smooth` | `line`, `path`, `point`, `ribbon` | Source `x` and raw numerical `y`; choose method/formula; ribbon requires SE bounds | `"identity"` or native PositionIdentity |
| `density` | `line`, `path`, `point`, `area` | Positions in `x`, optional weights; choose bandwidth; `after_stat(density)` is density per source-position unit | `"identity"` or native PositionIdentity |
| `ecdf` | `line`, `path`, `point`, `step` | Positions in `x`; `after_stat(ecdf)` is cumulative fraction; explicit `pad = FALSE` avoids infinite source anchors | `"identity"` or native PositionIdentity |

For a density area, use `position = "identity"` explicitly because area defaults to stacking. The `bounds` argument controls clipping at source-sequence boundaries. The table above lists geometries that can be projected.

Reuse the real `bins` and `track_seq` prepared above. This explicit linear fit is a display summary of the GC fractions, with a native confidence ribbon and mean line sharing one Track:

{{code:courseblockbi}}

Statistical fitting uses raw supplied source x/y. A later transformed Track value scale maps the estimates; it does not silently refit a log-space model, unlike placing a transform before a Stat in a native ggplot. If transformed-space fitting is intended, prepare and name that transformed variable explicitly and state its unit. Kernel density of event positions is not measured sequencing depth. ECDF output is a cumulative distribution, not a basewise alignment.

The directed source-backed checks compare native ggplot2 results and group membership for real pUC19 GC windows and G/C positions, plus curved/vertical projection and an explicitly missing-value control. `check-track-estimates.R` retains the numerical report and actual CI/mean image. Finite source bounds and physical-band constraints still apply; automatic placement does not change the fitted model or manufacture observations.

<a id="track-focus"></a>

## Tracks and windows share source coordinates {#track-focus}

Focus pUC19 to 500–1500 bp and trim the same observed GC windows as intervals. `focus_ggchord_track` does not recompute averages, invent boundary samples or bridge omitted windows. Return to actual DNA and declare new windows/denominators if finer GC summaries are needed.

{{code:comparetrackfocus}}
{{console:comparetrackfocus}}

## Compare event counts with supplied totals {#count-and-col}

The preceding `bar` counts HSPs at each query start. Here the totals are computed first and drawn with `col`. The `track` argument names the physical band and `position` controls how values are arranged within it.

```r
start_totals <- foundation_events |>
  dplyr::count(accver, pos, name = "hsp_count")
p_start_totals <- ggchord(data = foundation_seqs) +
  track_band("query_starts") +
  geom_seq() +
  geom_track(
    aes(x = pos, y = hsp_count), data = start_totals,
    geom = "col", track = "query_starts",
    position = "identity"
  ) +
  coord_collinear()
```
