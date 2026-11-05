<div class="task-intro"><p>Compute GC fractions from real DNA, then compare bars with a moving mean. Each track has an explicit unit, window and baseline.</p></div>

<span id="prepare"></span>

## Compute a track from actual DNA {#domain}
Compute GC fraction in 100 bp
windows of bundled pUC19. The final shorter window uses its actual base count.


{{code:trackdata,frontgctable}}

{{output:frontgctable}}

Each row of `bins` is an actual sequence window. Check its endpoints and `gc`, then map values to track height and colour.

### Draw GC fractions {#draw-gc}

{{code:tracks}}
{{figure:tracks}}

The domain is 0–1, in fractions. Replace it with the actual unit/range of other
measurements. `position_track()` controls placement/physical width;
`scale_track_value_continuous()` controls numerical mapping.

## Share an ID deliberately {#sharing}
Same-ID layers share a numerical domain, including bar baselines, and must agree
on width/placement. Colour scales are scoped independently by track ID. Tiles
show interval categories, lines connect observations, bars encode interval values
relative to a baseline. Sharing a domain does not require sharing a spatial band.

## Inspect retention and limits {#diagnostics}
The export's `track_registry` and `track_diagnostics` record domains, physical
sizes, native offsets and clipping. Omitted/clipped rows do not mean zero values.
Focus tracks with `focus_ggchord_track()` using the same windows.

## Combine bars and a moving mean {#multiple}
Use the same GC windows and add a three-window moving mean.

{{code:multitrack}}
{{figure:multitrack}}

Bars and line use the same actual GC measurements; the line is a three-window
moving mean. Incomplete end windows are omitted, not filled with invented values.
Line positions use integer bin midpoints. Both layers use the same track ID.

## Check your result {#check-result}
GC in bins lies between 0 and 1; each row is an actual window and the last may be shorter than 100 bp. Bars and line share the same 0–1 fraction scale.

<div class="next-step">Continue: <a href="output.html">Export, preview and diagnostics →</a></div>

## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

```r
ggsave("pUC19-tracks.png", p_multitrack, width = 9, height = 8, units = "in", dpi = 144)
```
