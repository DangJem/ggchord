Verify the shared grammar before tuning linear composition. These objects replace only the coordinates on the chord teaching objects: no table edits or new alignments.

## Prepare this page {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) first.
{{prepare:sequencepalette,docstyle,showcasedata}}
{{prepare:learnseq,learnlinks,learngenes,learnlabels}}

## Replace coord {#switch}
{{code:learnlinear}}

## Read backbone lengths {#backbones}
{{figure:linearseq}}
Rows share one bp scale. Different sequence lengths remain different rather than being stretched to equal width.

## Add comparisons {#links}
{{figure:linearlinks}}
The same observed HSPs and cross-row connections remain. This also exposes the limitations of default long routes, addressed on the next page.

## Add strand-separated genes and labels {#labels}
{{figure:linearlabels}}
Arrows and labels use the same selected table. Display flips do not rewrite source strands. Continue to [2.3 Route cross-row links](collinear-routing.html).
