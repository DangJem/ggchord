Each map retains all source features, independent reference primers and every site in its declared enzyme profile. No three-site sampling is used. Maps use their own complete sequence, the same layers and a 10 × 8 in canvas. Open the full-resolution image for annotation text.

## Preparation and complete display policy {#prepare}

{{code:courseblockg}}

{{code:referencebuilder}}

`reference$features`, `reference$primers` and `reference$sites` are the complete tables submitted to the plot. `reference$all_sites` retains raw matches; `reference$policy` identifies the enzyme profile. Feature geometry is not thinned and the site profile is not sampled to a fixed count.

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

## Enzyme profiles and provenance {#reference-selection}

Saved reference profiles take priority within `commercial_nonredundant`. The pETDuet-1 file saves None while its reference PNG shows sites, so `unique_6plus` is selected explicitly. This does not claim to recover the original PNG export settings. Every row of the chosen profile is retained. Complete means that declared profile, rather than every candidate from every enzyme database on one map.

Independent primers use `set = "reference"`. A source without independent reference primers has zero such rows; no records are invented. Named pieces remain in `segments` rather than becoming new top-level features. Enzymes at a shared position can share a label while their input records remain available.

## Related views {#next}

[pUC19 layers](plasmids.html) · [pETDuet-1 dense layout](circular-dense-case.html) · [Labels and output size](circular-labels.html)

<span id="gallery-limits"></span>
