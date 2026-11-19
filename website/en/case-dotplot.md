Four full-length references form a 4 × 4 matrix of possible ordered sequence pairs. The figure preserves all 32 HSPs from the bundled comparison.

<a id="map"></a>

## Draw the full map {#map}

The figure uses a 13 × 12 in canvas at 144 dpi. The script returns the plot so you can pass it to `ggplot2::ggsave()` at the same size.

{{example:case-dotplot}}

<a id="reading"></a>

## Read the figure {#reading}

Read the query accession from the column and the subject accession from the row. Within the cell, line colour shows identity, line width shows alignment length in bp, and the slope follows the order of the paired endpoints. The cell label reports how many input HSPs it contains.

<a id="inputs"></a>

## Source records and display choices {#inputs}

Each of the six unordered pairs has records in one query-to-subject direction. The dataset contains no self-comparisons or reciprocal copies. `flagship_dotplot_example$cell_labels` supplies the labels for all 16 cells at the displayed 4,000 bp gap. Recompute label positions if you change that gap or reorder an axis.

[Draw one HSP](dotplot.md) · [Dotplot axes](../../man/coord_dotplot.Rd)
