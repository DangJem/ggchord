Coverage answers a different question from a precomputed measurement. Given a set of alignments to a subject, what fraction of each query window lies inside at least one of those intervals? Start with a single subject before comparing several.

<span id="interval-coverage"></span>

<a id="coverage"></a>

## Compute window coverage for one subject {#coverage}

Select the three MT108731.1-to-MT118296.1 HSPs. For each HSP take the lower and upper query coordinates as interval bounds. Set `stat = "coverage"` and `binwidth = 10000` on the column layer. Overlaps count once toward the covered fraction. The rightmost window is only 4,323 bp long because the query ends at 64,323 bp.

{{example:coverage-pair}}

<a id="domains"></a>

## Compare subjects on the same height scale {#domains}

Keep all subjects for the same query, retain their IDs in the interval table and facet by subject. The MT118296.1 panel reproduces the single-subject result. The other panels show where the same query overlaps each of the two other subjects. Shared limits from 0 to 1 let the same bar height mean the same fraction in every panel.

{{example:coverage-facets}}

<span id="count-and-col"></span>

<a id="members"></a>

## Check the denominator and contributing HSPs {#members}

Name the coverage layer `windows`, export its statistics and inspect the rows for MT118296.1. The 20,001–30,000 bp window has 8,415 covered bases out of 10,000. The `source_rows` column identifies the contributing rows of the input interval table.

{{example:window-members}}

For point observations, `geom_track_bar()` counts records per position. For precomputed window measurements, [draw columns directly](tracks.md). To compare mean depth with coverage, map `after_stat(mean_depth)` to height on a separately labelled track. Mean depth includes multiple contributions from overlapping intervals.
