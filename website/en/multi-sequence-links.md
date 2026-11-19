<a id="intervals"></a>

## Draw a ribbon from each observed alignment {#intervals}

These three local alignments connect MT108731.1 and MT118296.1. Each record supplies four endpoints. The ribbon meets each backbone across that endpoint interval. A decreasing pair of coordinates can record a reverse hit. Keep the recorded order when preparing a ribbon table.

{{example:paired-ribbons}}

<a id="points"></a>

## Choose endpoints and visible pairs {#points}

Use `geom_link_line()` to connect two individual positions. Supply `qaccver, qpos, saccver, spos` when each end is a point. Use `geom_link_ribbon()` for a pair of intervals. The fixed fill above gives the three records the same colour. To encode percent identity, keep the input `pident` and use the default ribbon identity mapping.

In a collinear plot, adjacent rows are the default comparison. Set `pairs` to a table of the sequence pairs of interest, or use `adjacent_only = FALSE` to include nonadjacent comparisons. Check the [source rows](diagnostics.md) before interpreting a missing link.
