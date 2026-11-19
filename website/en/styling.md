<a id="categories"></a>

## Map a feature class to colour {#categories}

Each category in the source `type` column gets its own fill. Naming the palette entries keeps the same colour attached to each category after sorting or filtering. The named scale supplies the legend title and keys.

{{example:category-fill}}

<a id="styles"></a>

## Style a layer or the complete figure {#styles}

Use a fixed layer argument such as `fill = "#2B7580"` for one uniform colour. Inside `aes()`, map the source column to `fill` and choose a scale for its values. Each starting aesthetic reads from the layer's source table. If sequence placement uses `aes(accver = ref)`, map `fill = ref` as well to colour by that displayed reference.

Set `show.legend = FALSE` on a layer to hide all its guide contributions. On a point track, `show.legend = c(fill = FALSE)` suppresses just that layer's fill contribution. For plot-wide typography or legend position, use `ggplot2::theme()`. Use [track value scales](tracks.md) for numerical heights and a fill or colour scale for the palette.
