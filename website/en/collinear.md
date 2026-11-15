<div class="section-intro">The complete multi-sequence Collinear result: four real genomes, strand-separated features and outward-routed cross-row links.</div>


<a id="question"></a>

## Research question, alignment and tradeoffs {#question}

This view asks how selected CDS and observed intervals compare on one shared bp scale across four complete sequences. FASTA, GFF3 CDS and original BLAST HSPs are the same sources as the [chord case](chord.md#question). Auto orientation is a link-supported display decision. Center aligns complete sequences without asserting homologous centres.

Adjacent rows use explicitly filtered observations; each cross-row pair retains only its longest HSP. Omitting two crowded features is a display decision, with complete source tables retained. Right-side controls are derived from the current layout and must be recomputed after order/direction/gap changes. Local clearance avoids arrows, not every intervening row along a long path.

<a id="overview"></a>

## Complete result {#overview}

{{figure:showcasecollinear}}

Compare four backbones on a shared base-pair scale with 16 selected CDS placed on either side of the rows. Read gene positions along one row first, then follow a band to its matching interval. A solid outline marks adjacent rows, while a dashed outline marks rows with another displayed row between them. Start with the [four-backbone plot](multi-sequence-start.md#collinear) and [alignment and CDS example](multi-sequence-start.md#collinear) if you are new to the workflow.

<a id="prepare"></a>

## Prepare bundled data {#prepare}
Use bundled `genome_comparison_example` directly, without extracting files; original provenance remains in the [input bundle](https://github.com/DangJem/ggchord/archive/refs/heads/main.zip). The example requires dplyr and purrr:

Retrieve the display subsets selected by the original rules, then prepare shared colours and styles; source-file filtering does not need to be repeated.

{{code:sequencepalette,docstyle,showcasedata,chordshowcasegenes}}

Retain four full FASTA records, 266 CDS and 32 original HSPs. Adjacent rows use the explicitly selected observations; each non-neighbour pair retains only its longest HSP to reduce overlap.

<a id="gc"></a>

## Reserve space for connections and labels {#gc}
To compare sequence composition, see [tracks](tracks.md) for a GC track example.

<a id="compose"></a>

## Compose the layers {#compose}
{{code:showcasecollinear}}

Dashed outlines identify non-adjacent displayed rows, with sufficient linewidth to read. Right-side controls separate long routes; smooth local clearance retreats near arrows. Identity has the same left-side colourbar as the chord view, and the bottom strand keys use arrows. The ruler's horizontal offset is 0; its length still uses the shared bp conversion.

<a id="read"></a>

## Adjustment and limits {#read}
Rebuild controls when changing seq_gap, order, orientation or link selection. Local arrow clearance is not an automatic route around intervening whole sequences; ribbons may still cross legitimately. Use [focus/segments](views.md) for local questions. Explicit circular radius/curvature cannot currently be passed to linear coordinates; see [conversion boundaries](multi-sequence-views.md#input).

<a id="save"></a>

## Save at the intended size {#save}
{{code:courseblockf}}

<a id="evidence-check"></a>

## Evidence checks before saving {#evidence-check}

Preview with view_ggchord at the declared width/height and save with ggsave at matching dimensions. Check IDs, source strands, endpoint directions, selection reports, legend units and hidden-label reasons. A successful PNG does not certify font portability across PDF/SVG viewers. Keep complete source tables, display subsets, parameters, version and input hashes. Download the [asset/provenance record](../ASSET-SOURCES.md).

<a id="source-rights"></a>

## Source rights and reproducibility {#source-rights}

The package code is MIT-licensed; source database records and reference catalogues retain their original terms. The package license does not relicense third-party data. Keep accession identities and provenance with exported figures. The [asset/source record](../ASSET-SOURCES.md), input archive hashes and displayed preparation script document this example; algorithmic matches remain hypotheses for experimental verification.
