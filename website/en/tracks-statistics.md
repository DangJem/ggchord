{{code:courseblockbg}}

{{prepare:trackdata,tracks,foundationdata}}

## Count observations or display supplied values {#tracks}

`bar` counts observations at source positions; `col` displays supplied values.
Each accession/window computes its Stat and data Position before view projection
or display replication. `track=position_track(...)` selects the physical band;
`position` independently selects a native data adjustment such as stack or dodge.
The earlier `position=position_track(...)` band spelling remains accepted.

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

Available statistics are identity/count/bin/summary, plus the documented [smooth/density/ecdf estimation combinations](tracks-statistics.html#estimates). Bin needs explicit binwidth
or breaks. Native group and after_stat mappings are supported; numerical axes
retain the statistic's value units. Arbitrary third-party Geom/Stat plugins and
unsupported combinations are rejected. Existing
track line/tile shortcuts retain their documented behaviour.

Source position uses x orientation; native text size currently uses mm and polygon outlines are full. Other orientations/size units/partial outlines fail explicitly. Native text styles and path arrows are forwarded.

Native tracks currently require one panel and a ggchord coordinate system. Facets, `coord_flip()`, raster, boxplot and violin are outside this first supported set and are rejected explicitly. Use `coord_collinear(seq_axis = "vertical")` for vertical projection. Statistical source membership is preserved before dodge or other data positions. Rectangles and tiles use their supplied bounds without an artificial zero baseline.
## Two independent bands and display order {#independent-bands}

Encode the same GC by bar height and colour using explicit gc_height/gc_colour IDs. Set height and colour domains 0–1 separately, with distinct physical widths. The second figure changes only colour order from 2 to 0, moving that band nearer the backbone. Order is a finite priority, not a measurement or required one-based index; ties use IDs. Layers sharing an ID must agree on side/order/width/gap. Lines do not implicitly bridge independent windows.

{{code:comparemultitrackdata}}
{{code:comparetracktwo}}
{{figure:comparetracktwo}}
{{code:comparetrackorder}}
{{figure:comparetrackorder}}
## DNA windows and line boundaries {#window-policy}

GC here is `mean(bases %in% c("G", "C"))`: ambiguous bases remain in the denominator and contribute zero to the numerator. The final window uses its actual number of bases. This pUC19 input has no ambiguous bases; for other inputs explicitly choose and report a policy rather than assuming they are omitted. Focusing clips intervals and retains the original summary; it does not recompute GC for a shortened window. Line tracks connect the supplied samples within their grouping; `closed` explicitly controls circular closure and `gap` belongs to physical track placement. Do not use a closed line to imply continuity across omitted source windows.
## Native source-space estimates {#estimates}

`geom_track()` also delegates `smooth`, `density` and `ecdf` to native ggplot2 statistics **before projection**. These create estimates of supplied observations, not additional observations. They run independently by accession/window and explicit group. Estimation `source_rows` records the eligible group submitted to the estimator after native missing-value handling; it is not per-point influence, positive read support or a custom model's additional internal filtering report.

| Statistic | Documented geometry | Input and explicit choices | Data Position |
| --- | --- | --- | --- |
| `smooth` | `line`, `path`, `point`, `ribbon` | Source `x` and raw numerical `y`; choose method/formula; ribbon requires SE bounds | `"identity"` or native PositionIdentity |
| `density` | `line`, `path`, `point`, `area` | Positions in `x`, optional weights; choose bandwidth; `after_stat(density)` is density per source-position unit | `"identity"` or native PositionIdentity |
| `ecdf` | `line`, `path`, `point`, `step` | Positions in `x`; `after_stat(ecdf)` is cumulative fraction; explicit `pad = FALSE` avoids infinite source anchors | `"identity"` or native PositionIdentity |

Other combinations of these statistics fail explicitly. Use `position = "identity"` explicitly for density area because area otherwise defaults to stacking. `bounds` remains the Track source-boundary policy and does not forward the native density boundary-correction parameter. Facets, boxplots, violins, raster and arbitrary third-party Geom/Stat objects remain unsupported; a numerical name in a function is not a projection contract.

Reuse the real `bins` and `track_seq` prepared above. This explicit linear fit is a display summary of the GC fractions, with a native confidence ribbon and mean line sharing one Track:

{{code:courseblockbi}}

Statistical fitting uses raw supplied source x/y. A later transformed Track value scale maps the estimates; it does not silently refit a log-space model, unlike placing a transform before a Stat in a native ggplot. If transformed-space fitting is intended, prepare and name that transformed variable explicitly and state its unit. Kernel density of event positions is not measured sequencing depth. ECDF output is a cumulative distribution, not a basewise alignment.

The directed source-backed checks compare native ggplot2 results and group membership for real pUC19 GC windows and G/C positions, plus curved/vertical projection and an explicitly missing-value control. `check-track-estimates.R` retains the numerical report and actual CI/mean image. Finite source bounds and physical-band constraints still apply; automatic placement does not change the fitted model or manufacture observations.

## Tracks and windows share source coordinates {#track-focus}

Focus pUC19 to 500–1500 bp and trim the same observed GC windows as intervals. `focus_ggchord_track` does not recompute averages, invent boundary samples or bridge omitted windows. Return to actual DNA and declare new windows/denominators if finer GC summaries are needed.

{{code:comparetrackfocus}}
{{console:comparetrackfocus}}
