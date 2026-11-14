## Dense, complete annotation inputs {#complex}

This independent pETDuet-1 example demonstrates the interaction of complete annotations, primers and enzyme labels at a normal output size. It does not depend on objects from the pUC19 page.

{{code:courseblockg}}

{{code:referencebuilder}}

```r
reference <- make_reference_map(plasmid_example_pETDuet_1)
p <- reference$plot
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogduet}}

## What complete means here {#source-profile}

All 20 features, 5 independent reference primers and 58 site rows in the selected profile are supplied. These sites occupy 54 distinct positions. The source saves an empty enzyme profile, so the explicit policy is unique six-or-more-base commercial nonredundant cutters. Enzyme-profile selection and arbitrary display sampling are different operations.

## Label and geometry responsibilities {#layout}

Stacking determines feature bands. The label layer first uses available internal space and can move text outside; shared source positions remain fixed. A shared enzyme label can describe several rows at one position. Do not infer source counts from the number of text strings.

## More room without deleting data {#size}

The same plot can be exported on a wider canvas. This changes physical space; it does not alter features, selected sites, strands or complete length. Check label diagnostics when selecting an output size.

```r
ggsave("pETDuet-complete.png", p, width = 16.5, height = 8, units = "in", dpi = 144)
```

## Related functions {#next}

[Feature stacking](circular-features.html) · [Labels](circular-labels.html) · [All 13 reference maps](plasmid-catalog.html)
