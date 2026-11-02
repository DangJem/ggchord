Focus on the visual information added at each step. Preparation computes inputs; independent layers still control drawing.

## Prepare this page {#prepare}
{{prepare:featurepalette,docstyle,plasmidinput}}
{{code:learnplasmid}}

## Step 1: double backbone and centre name {#backbone}
{{figure:plasmidbackbone}}
Check name and full length first. `coord_circular()` accepts one sequence; the geom owns the double-line style.

## Step 2: feature shapes {#features}
{{figure:plasmidfeatures}}
`position_feature_stack()` separates overlaps; `geom_feature_plasmid()` selects shapes and colours by real category. Compact short-feature shapes do not imply missing data.

## Step 3: fitted feature labels {#labels}
{{figure:plasmidtext}}
Geometry and labels share a Position. Allow shrinking only above an explicit minimum and external placement when necessary; recheck after changing dimensions.

## Step 4: primers and selected restriction sites {#details}
{{code:plasmid}}
{{figure:plasmid}}
Only three selected enzyme sites are displayed; the full search table remains available. See [3.3](recognition.html) for recognition provenance and assumptions. Compare annotation densities across [13 references](plasmid-catalog.html), then finish the [dense pETDuet-1 showcase](plasmids.html#complex).
