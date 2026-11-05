These 13 maps use distinct repository FASTA records. The sequence and annotation sources remain available for each map.

## Prepare all cases once {#prepare}
Run `library(dplyr)` before using the shared plotting script.

{{recipe:catalog}}

Save and run the shared script, then copy the three lines beside a figure. All features are drawn. Their names are placed inside when possible and may move outside when space runs out. Restriction sites nearest three evenly spaced target positions supply the outer enzyme and position labels.

## Inputs and recognition results {#counts}
{{plasmidcatalog}}

## pSB1C3: Small vector and basic structure {#catalogpsb}
```r
sequence <- plasmid_example_pSB1C3
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpsb}}

## pUC19: Common features and replication origin {#catalogpuc}
```r
sequence <- plasmid_example_pUC19
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpuc}}

## pBR322: Multiple resistance regions and direction {#catalogpbr}
```r
sequence <- plasmid_example_pBR322
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpbr}}

## pBluescript II SK(+): Local regulation and short features {#catalogblue}
```r
sequence <- plasmid_example_pBluescript_II_SK_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogblue}}

## pET-28a(+): Expression regions and tags {#catalogpet}
```r
sequence <- plasmid_example_pET_28a_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpet}}

## pETDuet-1: Dual expression units; transition to the dense case {#catalogduet}
```r
sequence <- plasmid_example_pETDuet_1
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogduet}}

## pDONR221: Recombination-related architecture {#catalogdonr}
```r
sequence <- plasmid_example_pDONR221
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogdonr}}

## pcDNA3.1(+): Annotation hierarchy in an expression vector {#catalogpcdna}
```r
sequence <- plasmid_example_pcDNA3_1_plus
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogpcdna}}

## pTRE-Tight-BI: Bidirectional regulatory structure {#catalogtre}
```r
sequence <- plasmid_example_pTRE_Tight_BI
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogtre}}

## pTRIPZ: Display density on a longer sequence {#catalogtripz}
```r
sequence <- plasmid_example_pTRIPZ
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogtripz}}

## pCAMBIA1300: Architecture of a larger plant vector {#catalogcambia}
```r
sequence <- plasmid_example_pCAMBIA1300
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogcambia}}

## pEarleyGate 201: Category comparison with another plant vector {#catalogearley}
```r
sequence <- plasmid_example_pEarleyGate_201
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogearley}}

## pSpCas9(BB)-2A-GFP (PX458): Large expression vector and label fitting {#catalogcas}
```r
sequence <- plasmid_example_pSpCas9_BB_2A_GFP_PX458
p <- make_catalog_plot(sequence, find_common_features(sequence))
view_ggchord(p, width = 10, height = 8, units = "in")
```

{{figure:catalogcas}}

## Next {#next}
After comparing densities, continue to [The complete case](plasmids.html#complex), adding separate primer and restriction layers.
