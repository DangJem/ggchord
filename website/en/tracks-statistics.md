<span id="interval-coverage"></span>

<a id="coverage"></a>

## Compare window coverage across groups {#coverage}

Keep MT108731.1 as the query reference and split its observed HSPs by subject accession. Each panel shows the fraction of query bases covered by that subject in successive 10,000 bp windows. Overlapping HSPs count once toward the union. The common 0–1 axis makes the three panels directly comparable.

{{example:coverage-facets}}

<span id="count-and-col"></span>

<a id="members"></a>

## Read the denominator and contributing rows {#members}

A full window divides by 10,000 bases. The last window ends at base 64,323 and divides by 4,323. For subject MT118296.1, the window 20,001–30,000 contains 8,415 covered bases. The exported `source_rows` column identifies the input intervals that contributed.

{{example:window-members}}

For single-position observations, `geom_track_bar()` counts records at each position. When the input already has a count column, use `geom_track_col()` and map that column to height.

<a id="domains"></a>

## Choose the height domain {#domains}

For measurements from different biological samples, keep the sample ID in the interval table and facet by that column. The default shared height domain supports comparisons between samples. To examine variation within each panel, request `facet_domain = "panel"` on the track value scale, leave its numerical limits open and keep the local axes visible. For separate plots of the same measure, `shared_track_limits()` computes a range you can pass to each scale.

Mean depth counts overlapping intervals separately. Map `aes(y = after_stat(mean_depth))` on a second track to show it alongside coverage. The two axes have different units.
