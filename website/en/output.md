## Prepare this page {#prepare}
In a fresh R session, run the code below, then follow this page in order. Colours and layout are optional display choices; inputs select only the fields needed.

{{code:frontlinearbase}}
{{figure:frontlinearbase}}

## Save the intended size {#size}
The example uses `p_linear_base` from this page. If you already have a plot, use your plot object in its place.

{{code:frontpreview}}

For circular plots, omit both dimensions to estimate the preview size. For chord and collinear plots, supply a width and use `height = NULL` to estimate height. A webpage also scales a thumbnail to its content column, so increasing dpi alone does not enlarge its text.

After checking text, labels and legends, save at the same dimensions.

{{code:frontsave}}

`view_ggchord()` returns a temporary image path; it does not keep a final export.
The website, README and vignettes share PNG previews generated with
`view_ggchord()` at sizes listed in the [gallery](gallery.html). They do not capture an
interactively resized panel. Use your intended PDF/SVG device explicitly for
final vector output; font availability and embedding remain device-dependent.

## Check your result {#check-result}
Use matching dimensions for preview and saving. Inspect text, guides and warnings; a successful measurement does not certify fonts or release quality across devices.
