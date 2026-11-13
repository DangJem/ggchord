These 13 maps are an optional reference gallery, not a mandatory sequence of lessons. Complete [first circular maps](circular-start.html#input) and [annotations/labels](circular-layout.html#input), then choose a vector and density similar to your own.

These 13 maps use distinct repository FASTA records. The sequence and annotation sources remain available for each map.

## Prepare all cases once {#prepare}
Run `library(dplyr)` before using the shared plotting script.

{{code:courseblockg}}

{{recipe:catalog}}

Save and run the shared script, then copy the three lines beside a figure. All features are drawn. Their names are placed inside when possible and may move outside when space runs out. Restriction sites nearest three evenly spaced target positions supply the outer enzyme and position labels.

## Inputs and recognition results {#counts}
{{plasmidcatalog}}

## pSB1C3: Small vector and basic structure {#catalogpsb}
{{code:courseblockh}}

{{figure:catalogpsb}}

## pUC19: Common features and replication origin {#catalogpuc}
{{code:courseblocki}}

{{figure:catalogpuc}}

## pBR322: Multiple resistance regions and direction {#catalogpbr}
{{code:courseblockj}}

{{figure:catalogpbr}}

## pBluescript II SK(+): Local regulation and short features {#catalogblue}
{{code:courseblockk}}

{{figure:catalogblue}}

## pET-28a(+): Expression regions and tags {#catalogpet}
{{code:courseblockl}}

{{figure:catalogpet}}

## pETDuet-1: Dual expression units; transition to the dense case {#catalogduet}
{{code:courseblockm}}

{{figure:catalogduet}}

## pDONR221: Recombination-related architecture {#catalogdonr}
{{code:courseblockn}}

{{figure:catalogdonr}}

## pcDNA3.1(+): Annotation hierarchy in an expression vector {#catalogpcdna}
{{code:courseblocko}}

{{figure:catalogpcdna}}

## pTRE-Tight-BI: Bidirectional regulatory structure {#catalogtre}
{{code:courseblockp}}

{{figure:catalogtre}}

## pTRIPZ: Display density on a longer sequence {#catalogtripz}
{{code:courseblockq}}

{{figure:catalogtripz}}

## pCAMBIA1300: Architecture of a larger plant vector {#catalogcambia}
{{code:courseblockr}}

{{figure:catalogcambia}}

## pEarleyGate 201: Category comparison with another plant vector {#catalogearley}
{{code:courseblocks}}

{{figure:catalogearley}}

## pSpCas9(BB)-2A-GFP (PX458): Large expression vector and label fitting {#catalogcas}
{{code:courseblockt}}

{{figure:catalogcas}}

## Next {#next}
After comparing densities, continue to [The complete case](circular-dense-case.html#complex), adding separate primer and restriction layers.

## When comparing gallery maps {#gallery-limits}

Each map uses its own complete FASTA and source annotations; do not transfer another reference’s features or cuts to different DNA. Sites nearest three evenly spaced targets are a display selection, not evidence that other cuts are absent. Names may be internal, exterior or hidden while feature geometry remains. Inspect dense cases at their declared dimensions rather than editing source intervals to fit a reduced screenshot.


## Comparing reference selections {#reference-selection}

This gallery deliberately draws three site rows per map. For a saved enzyme
profile, apply `filter_restriction_sites(sites, set = "reference",
parent_set = "commercial_nonredundant")` to complete recognition results.
Saved profiles do not recover every old PNG export setting or enzyme-library
version: pETDuet-1 saves None despite showing cutters in the PNG. Independent
primers require `find_primer_bindings(sequence, set = "reference")` and primer
layers. Named subsegments remain in `segments`; they are not extra top-level
features, and dotted boundaries do not automatically receive extra labels.
Inspect dense maps again at their actual export dimensions.
