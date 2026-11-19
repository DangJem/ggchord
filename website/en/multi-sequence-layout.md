Use sequence layout to decide how a reference is displayed. Begin with the two-row comparison of MT108731.1 and MT118296.1, then try a change of order, a second copy of a reference, or equal visual lengths.

<a id="direction"></a>

## Reverse one displayed sequence {#direction}

The [original arrangement](collinear.md#adapt) has MT108731.1 on top, with both arrows pointing right. Place MT118296.1 first and assign a named orientation of −1 to MT108731.1. Now the longer backbone is below the shorter one and its arrow points left. All three HSPs keep their original coordinates.

{{example:order-orientation}}

<a id="length"></a>

## Compare relative positions on equal-length rows {#length}

Look at the first HSP on its own and set `length_mode = "equal"`. The query and subject rows now occupy the same visual width. Read the terminal ticks to see that the rows still represent 64,323 bp and 32,090 bp. The subject interval covers a larger proportion of its sequence.

{{example:layout-equal}}

<a id="instances"></a>

## Choose which displayed copy receives each ribbon {#instances}

Name three display rows `top`, `middle` and `bottom` in the order vector. The top and bottom rows are copies of MT108731.1. The first layer sends the longest HSP to the top copy. The second sends the remaining two to the bottom copy. In each layer, `between` names the two rows that receive the selected records.

{{example:layout-instances}}

Next, [compare interval ribbons and midpoint lines](multi-sequence-links.md). For spacing, orientation and per-sequence geometry, see [flexible parameters](parameters.md) and [the layout reference](../../man/layout_sequences.Rd).
