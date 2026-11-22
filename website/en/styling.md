The same five pUC19 features can answer two questions. A single fill shows where the intervals lie. Mapping fill to the recorded annotation type also tells you which class each interval belongs to.

<a id="styles"></a>

## Use a fixed fill to show the intervals {#styles}

Select the two CDS, two primer-binding intervals and the multiple cloning site. Set `fill` directly on the feature layer. The long and short intervals all have the same teal fill, including those separated into lanes near the origin.

{{example:fixed-features}}

<a id="categories"></a>

## Map the type column to a named palette {#categories}

Keep the same feature selection and placement. Put `fill = type` inside `aes()`, then name the colours in `scale_fill_manual()` with the three source values. The two CDS share one colour, as do the two primer-binding intervals. The third colour picks out the multiple cloning site.

{{example:category-fill}}

The palette names keep each colour attached to its type after sorting or filtering. To change typography or move the legend, use `theme()`. For differences between sequences in radius, offset or width, continue to [flexible parameters](parameters.md). To display measured values by height, try [a GC track](tracks.md).

Use `scale_colour_sequence()` for reference IDs, `scale_fill_strand()` for strand, `scale_fill_feature_type()` for annotation classes and `scale_fill_pident()` for percentage identity. If genes and ribbons both use `fill`, give the gene layer `id = "genes"` and add `scale_for("genes", scale_fill_strand())`. The gene legend shows strand directions while the ribbons keep their percentage scale.

The default canvas is warm ivory. Use `theme_ggchord_publication()` for a white panel. A ribbon has a fine green outline by default. Set `colour = NA` on its layer to omit the outline.
