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

## Encode values with a tile track {#tile}

Reuse `track_seq` and `bins`, computed from real DNA on this page. Tiles require `accver, start, end` and do not require numeric heights. Here the same `gc` values map to `track_fill`.

{{code:tracktile}}
{{figure:tracktile}}

Use the same 9 × 8 inch dimensions and 100 bp windows as the bars. Bars encode value by height; tiles keep a fixed band width and encode GC fraction from 0 to 1 by colour. This is not another experimental observation and does not average overlapping intervals. It compares encodings, **not a single geometry parameter**.

| Parameter or mapping | Default | Choice and constraints |
| --- | --- | --- |
| `position` | Required | `position_track("gc_colour")` declares an independent track |
| `width` | 4 mm in `position_track()` | 10 mm here; physical band width |
| `side` / `order` | outside / 1 | Lower order places a band nearer the backbone |
| `gap` | 1 mm | Physical spacing from the preceding band, not a sequence gap |
| `track_fill` | Never inferred | Map gc inside `aes()`; use `fill` for a fixed colour |
| `bounds` | `"error"` | Out-of-range coordinates error; clip reports clipping |
| `alpha` / `linewidth` | 1 / 0.25 | Fixed opacity and border width; border colour defaults to NA |

Use `scale_track_fill_gradientn()` for continuous colours or `scale_track_fill_manual()` for categories, with the matching track ID. Do not label a fraction as a percentage. To combine with the bars, add this tile layer to `p_tracks` and retain separate IDs for separate bands. Layers sharing an ID must agree on side/order/width/gap.

Change only tile width from 10 to 4 mm: colours stay the same while physical band width changes. Your own data must declare windows, units, NA handling and overlap policy. Check coordinate errors rather than hiding them with clip. See [`geom_track_tile()`](reference/geom_track_tile.html) and [`position_track()`](reference/position_track.html).

## Change only physical tile width {#tile-width}

Compare with the 10 mm band above. Keep the same 9 × 8 inch device, windows, colours and value domain, changing only `position_track(width)` to 4 mm:

{{code:comparetrackthin}}
{{figure:comparetrackthin}}

Colours still encode the same GC fractions. Only physical band width shrinks. The geom draws interval tiles and the Position controls band width; neither changes base composition. Try 6 mm while retaining export dimensions. Bar heights are separately controlled by `track_value` and its value scale; reducing band width is not numerical normalization.

## Check your result {#check-result}
GC in bins lies between 0 and 1; each row is an actual window and the last may be shorter than 100 bp. Bars and line share the same 0–1 fraction scale.

<div class="next-step">Continue: <a href="output.html">Export, preview and diagnostics →</a></div>

## Save this result {#save-result}
Save the plot just created to the working directory; use matching dimensions for preview.

```r
ggsave("pUC19-tracks.png", p_multitrack, width = 9, height = 8, units = "in", dpi = 144)
```

## Two independent bands and display order {#independent-bands}

Encode the same GC by bar height and colour using explicit gc_height/gc_colour IDs. Set height and colour domains 0–1 separately, with distinct physical widths. The second figure changes only colour order from 2 to 0, moving that band nearer the backbone. Order is a finite priority, not a measurement or required one-based index; ties use IDs. Layers sharing an ID must agree on side/order/width/gap. Lines do not implicitly bridge independent windows.

{{code:comparemultitrackdata}}
{{code:comparetracktwo}}
{{figure:comparetracktwo}}
{{code:comparetrackorder}}
{{figure:comparetrackorder}}

## DNA windows and line boundaries {#window-policy}

GC here is `mean(bases %in% c("G", "C"))`: ambiguous bases remain in the denominator and contribute zero to the numerator. The final window uses its actual number of bases. This pUC19 input has no ambiguous bases; for other inputs explicitly choose and report a policy rather than assuming they are omitted. Focusing clips intervals and retains the original summary; it does not recompute GC for a shortened window. Line tracks connect the supplied samples within their grouping; `closed` explicitly controls circular closure and `gap` belongs to physical track placement. Do not use a closed line to imply continuity across omitted source windows.
