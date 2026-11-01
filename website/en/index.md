<img class="package-logo" src="../assets/logo.webp" alt="ggchord" width="150" height="150">

Create **chord comparisons, circular maps and collinear plots** with ggplot2.
Add sequences, links, genes, labels and quantitative data as independent layers.

## Chord: compare several sequences {#chord}
{{figure:realchord}}

Links show existing alignments and outer arrows show selected annotations from
four original repository records. [Follow the complete workflow](real-data.html#chord).

## Circular: annotate a plasmid {#circular}
{{figure:plasmid}}

Combine pUC19 genes, replication origin, primers and selected restriction sites.
[Build a complete plasmid map](plasmids.html).

## Collinear: inspect relationships by row {#collinear}
{{figure:reallinear}}

Compare positions, directions and aligned intervals along sequence rows.
[Start from original files](real-data.html#linear), or try [three short examples](modes.html).

## Get started {#learn}
- [Install and load](installation.html), then inspect the [minimal tables](data.html#tables).
- Browse [tutorials](articles.html), [function documentation](reference.html) and the [gallery](gallery.html).
- Adjust [labels](labels.html), [tracks](tracks.html) and [colours](styling.html).
- Read the [migration guide](migration.html), [FAQ](faq.html) and [changelog](news.html).

## A short example {#quick-start}
Select only the columns used by this example from the bundled data before plotting.
The original table does not need to be passed with every metadata column.

{{code:setup}}
{{code:chord}}

Print `p_chord`, then [preview and save at explicit dimensions](output.html).
These are local 0.16.5 development docs; the formal release remains 0.15.0.
See [installation](installation.html#source) for the distinction.
