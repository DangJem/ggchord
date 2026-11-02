<div class="home-intro">
<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="320" height="320">
<p class="eyebrow">GENOMIC VISUALIZATION · R / GGPLOT2</p>
<p class="home-lead">Reveal the structure. Follow the connections.</p>
<p>Compose chord comparisons, collinear views and dense plasmid maps. Build reproducible genomic figures from real sequences, alignments and annotations.</p>
<div class="home-actions"><a class="primary-action" href="installation.html">Draw your first plot →</a><a href="articles.html">Find your workflow</a></div>
</div>

## Choose your visualization {#choose}
{{showcases}}

## Chord · Four genomes, one connected view {#chord}
{{figure:showcasechord}}

The original ggchord view, refined with signed curvature, varying radii, optimized display direction, local ribbon clearance and individually labelled selected features from real records.
[Data, complete code and interpretation →](chord.html)

## Collinear · Four rows, many layers {#collinear}
{{figure:showcasecollinear}}

A linear expression of the chord grammar, reusing real sequences, observed alignments and the identity palette. Features separate by strand; long cross-row links route to the right. The tutorial demonstrates coordinate-only conversion and explains current circular-only parameter limits.
[Complete workflow and the meaning of each layer →](collinear.html)

## Plasmid · A richly annotated pETDuet-1 {#circular}
{{figure:showcaseplasmid}}

A double backbone, layered features, five reference primers and restriction-enzyme names grouped by shared cut position. The input includes 58 cut positions and 20 features; labels are laid out for the actual export size.
[From a first circular view to complete plasmid annotation →](plasmids.html#complex)

## Continue from where you are {#learn}
- [Install and draw a first plot](installation.html)
- [Import your own data](formats.html)
- [Adjust colours, legends and themes](styling.html)
- [Browse all workflows and capabilities](articles.html)
- [Look up functions, arguments and results](reference/index.html)

Examples prefer dplyr, purrr, tidyr and the native `|>` pipe, with ggplot2 layer composition. Updating older code? Read the[migration guide](migration.html).
