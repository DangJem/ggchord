The same five pUC19 features can answer two questions. A single fill shows where the intervals lie. Mapping fill to the recorded annotation type also tells you which class each interval belongs to.

<a id="styles"></a>

## Use a fixed fill to show the intervals {#styles}

Select the two CDS, two primer-binding intervals and the multiple cloning site. Set `fill` directly on the feature layer. The long and short intervals all have the same teal fill, including those separated into lanes near the origin.

{{example:fixed-features}}

<a id="categories"></a>

## Map the type column to a named palette {#categories}

Keep the same feature selection and placement. Put `fill = type` inside `aes()`, then name the colours in `scale_fill_manual()` with the three source values. The two CDS share one colour, as do the two primer-binding intervals. The third colour picks out the multiple cloning site.

{{example:category-fill}}

The palette names keep each colour attached to its type after sorting or filtering. To change typography or move the legend, use `ggplot2::theme()`. For differences between sequences in radius, offset or width, continue to [flexible parameters](parameters.md). To display measured values by height, try [a GC track](tracks.md).
