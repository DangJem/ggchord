<a id="direction"></a>

## Reverse one displayed sequence {#direction}

Place MT118296.1 first and reverse MT108731.1. The accession names in the orientation vector let you reorder the rows later without moving the reversal to the other sequence. Compare the ends of the three ribbons with the [original arrangement](multi-sequence-links.md).

{{example:order-orientation}}

<a id="instances"></a>

## Display the same reference twice {#instances}

For an A–B–A comparison, name the two copies in the order vector, for example `layout_sequences(order = c(left = "A", middle = "B", right = "A"))`. In each ribbon layer, select the relevant source records and use `between = c("left", "middle")` or `between = c("middle", "right")`. Source coordinates and biological strand still belong to the original A. Both copies retain its complete length.

<a id="length"></a>

## Preserve or standardise visual lengths {#length}

Keep the default proportional lengths when comparing physical distances between sequences. Use `layout_sequences(length_mode = "equal")` to compare the relative positions of features on equal-length backbones. Read the local sequence axis for base or residue distances.
