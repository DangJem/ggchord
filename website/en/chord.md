<div class="section-intro">Bring four genomes into one connected view. Combine straight and curved backbones, strand-aware genes and the default stepped identity scale without hiding where the data came from.</div>


<a id="question"></a>

## Research question and evidence limits {#question}

This view asks where observed alignments lie across four complete source records and which original CDS occur near their endpoints. It is not a phylogenetic tree and does not independently establish inversions or function. Check the accession records at NCBI: [MT108731.1](https://www.ncbi.nlm.nih.gov/nuccore/MT108731.1), [MT118296.1](https://www.ncbi.nlm.nih.gov/nuccore/MT118296.1), [OQ646790.1](https://www.ncbi.nlm.nih.gov/nuccore/OQ646790.1), [OR222515.1](https://www.ncbi.nlm.nih.gov/nuccore/OR222515.1).

Length/identity thresholds, at most three HSPs per pair and strand-aware CDS sampling are display choices. Original FASTA/GFF3/BLAST records stay unchanged. Rotation, curvature, radius and reversal affect display only; labels retain source names. Consult complete GFF3 before reporting gene presence or absence.

<a id="overview"></a>

## What this view shows {#overview}
{{figure:showcasechord}}

Compare the dark blue `MT118296.1` and teal `OQ646790.1` backbones first. Their radii are 0.82 and 1.18. Orange `OR222515.1` and red `MT108731.1` use 1.02 and 1.08. Radius separates backbones visually. Read sequence lengths from the bp ticks. [Compare radius on its own](chord-geometry.md).

CDS are sampled separately for each available source strand, without a length threshold. Neither MT record has a negative-strand CDS. The sole negative-strand CDS on OQ646790.1 is 207 bp. An explicit Position places positive-strand arrows inside and negative-strand arrows outside. Arrow direction, colour and data retain the source strands. Positive-strand sampling favours CDS that overlap displayed link endpoints so the view demonstrates local ribbon clearance. Every arrow uses its original annotation as the label, including hypothetical protein.

<a id="prepare"></a>

## Prepare bundled data {#prepare}
This example uses bundled `genome_comparison_example` directly, without extracting files; original records remain traceable through the [input bundle](https://github.com/DangJem/ggchord/archive/refs/heads/main.zip). The example requires dplyr and purrr.

Retrieve the display subsets selected by the original rules, then prepare shared colours and styles; source-file filtering does not need to be repeated.

{{code:sequencepalette,docstyle,showcasedata,chordshowcasegenes}}

The input contains four unchanged FASTA records, their GFF3 annotations and 32 original HSPs. The display keeps HSPs at least 900 bp long and at least 75% identical, with at most three per sequence pair. This reduces clutter without merging or inventing alignments.

<a id="compose"></a>

## Compose the layers {#compose}
Endpoint-overlap sampling was completed when generating the data; the following code directly composes the layers.

{{code:showcasechord}}

| Setting | Purpose in this plot |
| --- | --- |
| `seq_orientation`, `seq_curvature`, `seq_radius` | Adjust display direction, bend and radius per sequence. Zero curvature draws a straight backbone. |
| `seq_gap = 0.06`, `rotation = 50` | Add space between backbones and turn the whole layout to position the elements. |
| `link_avoid = "smooth"`, `alpha = 0.46` | Leave smooth clearance near gene arrows and show overlaps through translucent ribbons. |
| `position = chord_position`, `width = 0.025` | Place arrows on opposite sides by source strand with offsets of ±0.030. Narrower arrows sit close to the backbones, and labels use the same position. |
| `label_layout = "radial"`, `label_wrap = 24` | Arrange gene labels radially and wrap longer names. |
| `pident` and the default stepped scale | Map percent identity to ribbon colour and show the default stepped scale. |

<a id="read"></a>

## Read and adapt the result {#read}
Start with the sequence colours and gene directions, then trace the larger aligned intervals. A ribbon crossing is not evidence of a rearrangement by itself. Read alignment identity using the default stepped colour scale.

Use [link filtering](ribbons.md) to change the display threshold, [labels](labels.md) to choose annotation density, or [sequence parameters](parameters.md) to change radius and curvature. For a file-by-file introduction, start with [original sequence comparisons](real-data.md).

<a id="save"></a>

## Save at its intended size {#save}
{{code:courseblockd}}

<a id="evidence-check"></a>

## Evidence checks before saving {#evidence-check}

Preview with view_ggchord at the declared width/height and save with ggsave at matching dimensions. Check IDs, source strands, endpoint directions, selection reports, legend units and hidden-label reasons. A successful PNG does not certify font portability across PDF/SVG viewers. Keep complete source tables, display subsets, parameters, version and input hashes. Download the [asset/provenance record](../ASSET-SOURCES.md).

<a id="source-rights"></a>

## Source rights and reproducibility {#source-rights}

The package code is MIT-licensed; source database records and reference catalogues retain their original terms. The package license does not relicense third-party data. Keep accession identities and provenance with exported figures. The [asset/source record](../ASSET-SOURCES.md), input archive hashes and displayed preparation script document this example; algorithmic matches remain hypotheses for experimental verification.
