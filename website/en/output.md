<a id="size"></a>

## Save a figure at its final width {#size}

This example saves a four-sequence figure at 160 × 110 mm. It writes a PDF named `sequence-backbones.pdf` in the current R session's temporary directory. Set `output_file` to a permanent path for a manuscript. Pass the plot explicitly so saving does not depend on which figure was drawn last.

{{example:export-size}}

<span id="layout-results"></span>

<a id="preview"></a>

## Inspect text and spacing at that size {#preview}

For a static screen preview, pass the same plot, width, height and units to `view_ggchord()`. Its `dpi` controls the raster resolution. Use `ggsave()` for the final file and choose the extension for the intended device. When exporting PNG for print, set `dpi` explicitly and retain the required physical size.

Check the last labels around the edges, the space between panel titles and axes, and the legend keys. Font metrics and available space determine which labels fit. For a recorded measurement, see [layout diagnostics](diagnostics.md).
