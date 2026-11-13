## Complete example: a dense pETDuet-1 map {#complex}

Use the pUC19 example above to learn the layers, then move to the richer pETDuet-1 map. Inputs come from the bundled sequence, feature and reference-primer catalogues. The restriction subset retains unique six-base-or-longer recognition sites from commercial_all; it is not every cut made by every enzyme.

This example runs independently. Put the preparation script in your working directory.

{{recipe:plasmid}}

{{code:frontduetsummary}}

{{output:frontduetsummary}}

{{code:showcaseplasmid}}

{{figure:showcaseplasmid}}

The adjacent code renders this 16.5 × 8 inch preview. Restriction, primer and selected feature labels share the outer region, and dense sites can push their labels farther from the backbone. The coordinate system fits the labels automatically. Check the leaders and text again if you reduce the preview size.

Inputs contain 20 features, five reference primers and 58 grouped cut positions. The feature preset retains plum coding regions, mint resistance genes, yellow origins and blue regulatory features; primers are purple, consistently across this site. Dense labels depend on export size; open the full image to read them.

The explicit viewport keeps the dense map at the documented export size. `label_fit = "shrink"` allows short-feature labels to fit down to 2.5 mm; labels that do not fit remain callouts. Preview and save at the same 16.5 × 8 in size.

{{code:courseblockal}}
