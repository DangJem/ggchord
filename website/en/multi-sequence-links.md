An HSP supplies one interval on each of two sequences. Choose a ribbon if both spans matter. When you only need the location of each match, a line can join a representative position from each interval.

<a id="intervals"></a>

## Keep both interval spans {#intervals}

The same directed pair has three HSPs. Map each one to a ribbon using its original endpoints. A fixed fill gives them all one colour, so it is easy to compare the widths. The shortest ribbon meets MT118296.1 close to base 1.

{{example:paired-ribbons}}

<a id="points"></a>

## Join the midpoints of the same three records {#points}

Compute one midpoint on the query and one on the subject, rounding down to whole bases. Pass these as `qpos` and `spos` to `geom_link_line()`. The rightmost line connects the midpoints of the longest HSP. With only the two midpoints the line no longer shows the interval widths or endpoint order.

{{example:link-midpoints}}

For more than two rows, the collinear ribbon layer connects adjacent rows by default. Set `adjacent_only = FALSE` to show selected nonadjacent pairs. Next, [inspect the input records of a layer](diagnostics.md) or [route records to repeated copies of a reference](multi-sequence-layout.md#instances).
