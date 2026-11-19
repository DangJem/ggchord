<a id="membership"></a>

## Read existing group memberships {#membership}

The installed file contains group memberships from public plant data. Each output row identifies a group and one member.

```r
groups <- system.file("extdata", "plant-feature-groups.txt", package = "ggchord")
members <- read_feature_relations(groups, format = "named_groups")
head(members$membership[c("group_id", "feature_id")])
```

<a id="links"></a>

## Check both source IDs {#links}

The pair table retains source IDs from Arabidopsis and soybean. Supply feature coordinates to obtain intervals for each end.

```r
pairs <- system.file("extdata", "plant-feature-pairs.homology", package = "ggchord")
positions <- system.file("extdata", "plant-pair-features.gff", package = "ggchord")
anchors <- read.table(positions, col.names = c("accver", "feature_id", "start", "end"))
joined <- read_feature_relations(pairs, feature_data = anchors)
head(joined$ribbon_data[c("qfeature_id", "sfeature_id", "qstart", "sstart")])
```

<a id="focus"></a>

## Choose local windows {#focus}

Before plotting, obtain complete sequence lengths independently, then choose windows around both anchors. See [local views](views.md).

<a id="projection"></a>

## Protein and genome coordinates {#projection}

Existing membership records identify related features. For precise amino-acid to nucleotide coordinates, see the [CDS projection example](projection.md#cds).

<a id="check-result"></a>

## Check the relationships {#check-result}

Check source ID matches and the coordinate units at both ends. The score in the third input column retains its original source meaning.

<a id="save-result"></a>

## Export matched relationships {#save-result}

The returned ribbon_data is an ordinary table that can be saved with both endpoint intervals.
