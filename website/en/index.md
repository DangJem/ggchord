<div class="home-intro">
<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="320" height="320">
<p class="eyebrow">GENOMIC VISUALIZATION · R / GGPLOT2</p>
<p class="home-lead">Genome comparisons and plasmid maps in R</p>
<p>Compose chord comparisons, collinear views and dense plasmid maps. Build reproducible genomic figures from real sequences, alignments and annotations.</p>
<div class="home-actions"><a class="primary-action" href="installation.html">Get started →</a><a href="about.html">Browse the chapters</a></div>
</div>

## Choose your visualization {#choose}
- [Chord plot](chord-start.html) for corresponding intervals across sequences.
- [Collinear plot](collinear-start.html) for comparing gene positions and order at a shared scale.
- [Plasmid map](plasmid-start.html) for annotating features on one sequence.

## Chord · Four genomes, one connected view {#chord}

{{figure:showcasechord}}

Positive-strand arrows sit inside the backbones, where ribbons leave clearance around them. Negative-strand arrows sit outside. Ribbon colour encodes alignment identity. The two MT records have only positive-strand CDS in the source annotations. The other two display both strands. [Data, selection and complete code](chord.html).

## Collinear · Four rows, many layers {#collinear}

{{figure:showcasecollinear}}

Compare gene order and aligned positions on a common base-pair scale. Arrows above and below the backbones distinguish strands. Dashed outlines mark links across rows. [First collinear lesson](collinear-start.html) · [Complete case](collinear.html).

## Plasmid · A richly annotated pETDuet-1 {#circular}
{{figure:showcaseplasmid}}

A double backbone, layered features, five reference primers and restriction-enzyme names grouped by shared cut position. The input includes 58 cut positions and 20 features; labels are laid out for the actual export size.
[From a first circular view to complete plasmid annotation →](plasmid-start.html) · [Complete case](plasmids.html#complex)

## Continue from where you are {#learn}
- [Import your own data](formats.html)
- [Adjust colours, legends and themes](styling.html)
- [Browse all workflows and capabilities](about.html)
- [Look up functions, arguments and results](reference/index.html)
