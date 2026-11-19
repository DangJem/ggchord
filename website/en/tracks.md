Start with the GC percentages already measured in the bundled pUC19 bedGraph. A quantitative track gives those values a height axis alongside the sequence. The examples use the optional dplyr package for table preparation.

<a id="domain"></a>

## Draw one bar per measured window {#domain}

Read the installed file with `read_bedgraph()`. The first three 100 bp windows have GC contents of 60%, 57% and 59%. Declare a band named `gc`, draw one column from each window’s value and give the band a 0–100 percentage axis. The first bar begins at base 1 and ends at base 100.

{{example:gc-bars}}

<a id="shared-band"></a>

## Add a line through the same window values {#shared-band}

Compute the midpoint of each window. Keep the column layer and add a line using the midpoint for x and the measured value for y. Both layers name the same track, so the line passes the middle of each bar’s top edge.

{{example:gc-shared-band}}

The band’s `width` reserves 15 mm of drawing space. Its numerical limits set the heights on that space. Increasing the band width makes a one-percentage-point change easier to see while preserving every measured percentage. Next, [compute coverage from alignment intervals](tracks-statistics.md) or [inspect the original bedGraph columns](import-signals.md).
