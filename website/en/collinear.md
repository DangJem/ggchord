<div class="section-intro">The complete multi-sequence Collinear result: four real genomes, strand-separated features and outward-routed cross-row links.</div>


## Research question, alignment and tradeoffs {#question}

This view asks how selected CDS and observed intervals compare on one shared bp scale across four complete sequences. FASTA, GFF3 CDS and original BLAST HSPs are the same sources as the [chord case](chord.html#question). Auto orientation is a link-supported display decision. Center aligns complete sequences without asserting homologous centres.

Adjacent rows use explicitly filtered observations; each cross-row pair retains only its longest HSP. Omitting two crowded features is a display decision, with complete source tables retained. Right-side controls are derived from the current layout and must be recomputed after order/direction/gap changes. Local clearance avoids arrows, not every intervening row along a long path.

## Complete result {#overview}

{{figure:showcasecollinear}}

Compare four backbones on a shared base-pair scale with 16 selected CDS placed on either side of the rows. Read gene positions along one row first, then follow a band to its matching interval. A solid outline marks adjacent rows, while a dashed outline marks rows with another displayed row between them. Start with the [four-backbone plot](multi-sequence-start.html#collinear) and [alignment and CDS example](multi-sequence-start.html#collinear) if you are new to the workflow.

## Prepare original records {#prepare}
Extract the [input bundle](../assets/tutorial-inputs.zip) and install dplyr, purrr and tidyr first:

Put the preparation script in your working directory. It contains the input, selection and display settings.

{{recipe:collinear}}

Retain four full FASTA records, 266 CDS and 32 original HSPs. Adjacent rows use the explicitly selected observations; each non-neighbour pair retains only its longest HSP to reduce overlap.

## Reserve space for connections and labels {#gc}
To compare sequence composition, see [tracks](tracks.html) for a GC track example.

## Compose the layers {#compose}
{{code:showcasecollinear}}

Dashed outlines identify non-adjacent displayed rows, with sufficient linewidth to read. Right-side controls separate long routes; smooth local clearance retreats near arrows. Identity has the same left-side colourbar as the chord view, and the bottom strand keys use arrows. The ruler's horizontal offset is 0; its length still uses the shared bp conversion.

## Adjustment and limits {#read}
Rebuild controls when changing seq_gap, order, orientation or link selection. Local arrow clearance is not an automatic route around intervening whole sequences; ribbons may still cross legitimately. Use [focus/segments](views.html) for local questions. Explicit circular radius/curvature cannot currently be passed to linear coordinates; see [conversion boundaries](multi-sequence-views.html#input).

## Save at the intended size {#save}
{{code:courseblockf}}

## Evidence checks before saving {#evidence-check}

Preview with view_ggchord at the declared width/height and save with ggsave at matching dimensions. Check IDs, source strands, endpoint directions, selection reports, legend units and hidden-label reasons. A successful PNG does not certify font portability across PDF/SVG viewers. Keep complete source tables, display subsets, parameters, version and input hashes. Download the [asset/provenance record](../assets/ASSET-SOURCES.md).

## Source rights and reproducibility {#source-rights}

The package code is MIT-licensed; source database records and reference catalogues retain their original terms. The package license does not relicense third-party data. Keep accession identities and provenance with exported figures. The [asset/source record](../assets/ASSET-SOURCES.md), input archive hashes and displayed preparation script document this example; algorithmic matches remain hypotheses for experimental verification.
