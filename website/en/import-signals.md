<a id="coordinates"></a>

## Inspect the first window {#coordinates}

This pUC19 bedGraph records the GC percentage of each 100 bp window, with a shorter final window. The first file interval is `[0, 100)`. The reader preserves those boundaries in `chromStart, chromEnd` and returns `start = 1, end = 100` for plotting. Its value 60 means 60% GC.

{{example:installed-signal}}

<a id="display"></a>

## Choose a geometric representation {#display}

To draw bars that span each measured interval, map its midpoint to `x`, its inclusive length to `width`, and the measurement to `y` in `geom_track_col()`. For a colour strip, use `geom_track_tile()` with `fill = value` and specify the position and width of each tile. Keep the value's measurement units attached to the scale. [Overlay a line and bars](tracks.md) for a complete plot using these data.

`read_bigwig()` returns the same core columns from a binary signal file. Read the signal description to determine its units. For the fraction of a window occupied by a collection of intervals, supply the original intervals to the [coverage statistic](tracks-statistics.md).
