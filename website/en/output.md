<a id="prepare"></a>

## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed.

{{code:frontlinearbase}}
{{figure:frontlinearbase}}

<a id="size"></a>

## Save the intended size {#size}
The example uses `p_linear_base` from this page. If you already have a plot, use your plot object in its place.

{{code:frontpreview}}

For circular plots, omit both dimensions to estimate the preview size. For chord and collinear plots, supply a width and use `height = NULL` to estimate height. A webpage also scales a thumbnail to its content column, so increasing dpi alone does not enlarge its text.

After checking text, labels and legends, save at the same dimensions.

{{code:frontsave}}

`view_ggchord()` returns a temporary image path; it does not keep a final export.
The website, README and vignettes share PNG previews generated with
`view_ggchord()` at sizes listed in the [gallery](gallery.md). They do not capture an
interactively resized panel. Use your intended PDF/SVG device explicitly for
final vector output; font availability and embedding remain device-dependent.

<a id="check-result"></a>

## Result checks {#check-result}
Use matching dimensions for preview and saving. Inspect text, guides and warnings; a successful measurement does not certify fonts or release quality across devices.

<a id="formats"></a>

## Choose a file format {#formats}

PNG serves web pages and slides; PDF/SVG serve vector-output workflows. The extension selects the default device unless a device is explicitly supplied. Reuse the target dimensions above and check fonts and embedding in the final file.

```r
ggsave("comparison.pdf", p_linear_base,
  width = output_width, height = output_height, units = "in")
```

<a id="layout-results"></a>

## Choose an output by question {#layout-results}

| Question | Call and result | Rows and units |
| --- | --- | --- |
| Where did an input link go? | `connection_report` and `original_data` from `export_ggchord_layout(p, original_data = TRUE)` | Look up the row in its layer input. A statistic can have several members in `source_rows`. Polygon rows describe vertices. Count observations in the source table. |
| Where is a base or residue drawn? | `map_ggchord_position()` returns `position`, `x`, `y` and `sequence_unit` | Source positions are in bp or aa. Display coordinates use layout units. Select an axis for a dot plot. |
| Is the plot legible at the intended size? | `measurement` and diagnostics from `diagnose_ggchord()` | Check panel measurement status on the target device before reading bounds in mm. In ordinary exports, `render_context` identifies layout estimates or existing measurements. |
| How do the panels line up? | `panels` and child layouts from `export_ggchord_fusion()` | Global composition geometry uses mm. Look up source positions in the corresponding child panel. |
