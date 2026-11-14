## Complete pUC19 inputs {#input}

This page explains the composition of one complete reference map. pETDuet-1 and its dense-label choices are described [separately](circular-dense-case.html).

{{code:courseblockg}}

{{code:referencebuilder}}

```r
reference <- make_reference_map(plasmid_example_pUC19)
p <- reference$plot
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpuc}}

## Read the layers {#layers}

The double backbone carries source positions. All reference features retain their original intervals, strands and segments. Feature labels share the stacking Position; independent primers and restriction sites use their own tables and layers. The centre gives identity and complete length.

## Inspect the complete tables {#sites}

```r
reference$policy
nrow(reference$features)
nrow(reference$primers)
nrow(reference$sites)
```

The chosen source profile contains 9 features and 39 restriction rows at 33 distinct positions. It has no independent reference-primer rows. Raw enzyme matches remain in `reference$all_sites`; the declared source profile is not a three-site display subset.

## Use your own annotations {#adapt}

Prepared annotation, primer and site tables can go directly to independent layers. Sequence recognition is optional when the tables already exist. Full lengths and source units must be supplied independently.

## Save at the same size {#save-result}

```r
ggsave("pUC19-complete.png", p, width = 10, height = 8, units = "in", dpi = 144)
```


<span id="check-result"></span><span id="draw"></span><span id="evidence-check"></span><span id="limits"></span><span id="question"></span><span id="source-rights"></span>
