See [Flagship data sources and construction](flagship-data.md) for provenance, display selections and rebuilding. This page loads bundled inputs directly and focuses on plotting settings.

Each map retains all source features, independent reference primers and every site in its declared enzyme profile. No three-site sampling is used. Maps use their own complete sequence, the same layers and a 10 × 8 in canvas. Open the full-resolution image for annotation text.

<a id="prepare"></a>

## Preparation and complete display policy {#prepare}

{{code:courseblockg}}

{{code:referencebuilder}}

`reference$features`, `reference$primers` and `reference$sites` are the complete tables actually drawn. `reference$all_site_count` records the number of original candidates, with complete search results reproducible through the generator; `reference$policy` describes the enzyme set. Feature geometry is not reduced and site sets are not sampled to a fixed count.

<a id="counts"></a>

## Inputs and recognition results {#counts}
{{plasmidcatalog}}

<a id="catalogpsb"></a>

## pSB1C3: Small vector and basic structure {#catalogpsb}
{{code:courseblockh}}

{{figure:catalogpsb}}

<a id="catalogpuc"></a>

## pUC19: Common features and replication origin {#catalogpuc}
{{code:courseblocki}}

{{figure:catalogpuc}}

<a id="catalogpbr"></a>

## pBR322: Multiple resistance regions and direction {#catalogpbr}
{{code:courseblockj}}

{{figure:catalogpbr}}

<a id="catalogblue"></a>

## pBluescript II SK(+): Local regulation and short features {#catalogblue}
{{code:courseblockk}}

{{figure:catalogblue}}

<a id="catalogpet"></a>

## pET-28a(+): Expression regions and tags {#catalogpet}
{{code:courseblockl}}

{{figure:catalogpet}}

<a id="catalogduet"></a>

## pETDuet-1: Dual expression units; transition to the dense case {#catalogduet}
{{code:courseblockm}}

{{figure:catalogduet}}

<a id="catalogdonr"></a>

## pDONR221: Recombination-related architecture {#catalogdonr}
{{code:courseblockn}}

{{figure:catalogdonr}}

<a id="catalogpcdna"></a>

## pcDNA3.1(+): Annotation hierarchy in an expression vector {#catalogpcdna}
{{code:courseblocko}}

{{figure:catalogpcdna}}

<a id="catalogtre"></a>

## pTRE-Tight-BI: Bidirectional regulatory structure {#catalogtre}
{{code:courseblockp}}

{{figure:catalogtre}}

<a id="catalogtripz"></a>

## pTRIPZ: Display density on a longer sequence {#catalogtripz}
{{code:courseblockq}}

{{figure:catalogtripz}}

<a id="catalogcambia"></a>

## pCAMBIA1300: Architecture of a larger plant vector {#catalogcambia}
{{code:courseblockr}}

{{figure:catalogcambia}}

<a id="catalogearley"></a>

## pEarleyGate 201: Category comparison with another plant vector {#catalogearley}
{{code:courseblocks}}

{{figure:catalogearley}}

<a id="catalogcas"></a>

## pSpCas9(BB)-2A-GFP (PX458): Large expression vector and label fitting {#catalogcas}
{{code:courseblockt}}

{{figure:catalogcas}}

<a id="reference-selection"></a>

## Enzyme profiles and provenance {#reference-selection}

Saved reference profiles take priority within `commercial_nonredundant`. The pETDuet-1 file saves None while its reference PNG shows sites, so `unique_6plus` is selected explicitly. This does not claim to recover the original PNG export settings. Every row of the chosen profile is retained. Complete means that declared profile, rather than every candidate from every enzyme database on one map.

Independent primers use `set = "reference"`. A source without independent reference primers has zero such rows; no records are invented. Named pieces remain in `segments` rather than becoming new top-level features. Enzymes at a shared position can share a label while their input records remain available.

<a id="next"></a>

## Related views {#next}

[pUC19 layers](plasmids.md) · [pETDuet-1 dense layout](circular-dense-case.md) · [Labels and output size](circular-labels.md)

<span id="gallery-limits"></span>
