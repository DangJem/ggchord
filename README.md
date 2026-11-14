<!-- Generated from website/readme-template.md and website/readme-text.json. -->
<p align="center">
  <img src="man/figures/logo.png" alt="ggchord" width="280">
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord: sequences, annotations and relationships

Build layered sequence figures in R with ggplot2: compare multiple sequences, annotate a circular map, inspect HSP alignments and connect independent views.

[Website](https://dangjem.github.io/ggchord/) · [Usage guide](website/en/modes.md) · [Figure gallery](website/en/gallery.md) · [Report an issue](https://github.com/DangJem/ggchord/issues)

<a id="modules"></a>
## Four plotting modules

| Module | Use it to | Views |
| --- | --- | --- |
| [I · Multi-sequence relationships](website/en/multi-sequence.md) | Compare sequence intervals, gene order and observed relationships. | Chord · Collinear · Radial · Sankey |
| [II · Single-sequence circular maps](website/en/circular.md) | Place features, primers and restriction sites on one sequence. | Circular |
| [III · HSP Dotplot](website/en/dotplot.md) | Inspect alignment intervals on query and subject axes. | Dotplot |
| [IV · Multi-view composition](website/en/fusion.md) | Connect independent plots through explicit source intervals. | Fusion |

The multi-sequence module has four Coord types. Circular and Dotplot add two more; Fusion is the composition layer. Layers, mappings, scales, positions and themes retain their ggplot2 roles.

<a id="examples"></a>
## Real-data figures

### Chord · relationships across four genomes

Combine mixed-curvature backbones, strand-aware CDS arrows, selected source labels and Identity-coloured ribbons. Features and HSPs are display subsets of the original records.

![Four-genome Chord flagship with labelled genes and alignment ribbons](man/figures/showcasechord.png)

[Complete example](website/en/chord.md) · [Chord / Collinear / Radial / Sankey](website/en/multi-sequence-views.md) · [Plot source](tools/documentation-examples.R)

### Collinear · order, orientation and cross-row links

Compare four genomes on a shared base-pair scale. Selected CDS sit on either side of each backbone; dashed outlines distinguish links across non-adjacent rows.

![Four-genome Collinear flagship with strand-separated genes and cross-row ribbons](man/figures/showcasecollinear.png)

[Complete example](website/en/collinear.md) · [Alignment and arrangement](website/en/collinear-layout.md) · [Plot source](tools/documentation-examples.R)

### Circular · dense pETDuet-1 annotations

A complete circular map with 20 features, five reference primers and 58 grouped cut positions from the declared enzyme set. Original sequence positions and annotation identities are retained.

![pETDuet-1 circular flagship with features, primers and restriction sites](man/figures/showcaseplasmid.png)

[Complete example](website/en/circular-dense-case.md) · [13 complete reference maps](website/en/plasmid-catalog.md) · [Plot source](tools/documentation-examples.R)

### Alignment matrices and linked views

<table>
  <tr>
    <th>III · HSP Dotplot</th>
    <th>IV · Multi-view composition</th>
  </tr>
  <tr>
    <td width="50%"><a href="man/figures/flagdotplot.png"><img src="man/figures/flagdotplot.png" alt="Four-sequence matrix of 32 original HSPs" width="460"></a></td>
    <td width="50%"><a href="man/figures/flagfusion.png"><img src="man/figures/flagfusion.png" alt="Linked Chord, Collinear, Radial and local-detail views" width="460"></a></td>
  </tr>
  <tr>
    <td>32 original HSPs; Identity by colour, reported alignment length by line width.</td>
    <td>Global relationships, ordered comparisons and a magnified source interval in one composition.</td>
  </tr>
  <tr>
    <td><a href="website/en/dotplot-start.md">Axes, direction and mappings →</a></td>
    <td><a href="website/en/fusion-details.md">Local magnification and character detail →</a></td>
  </tr>
</table>

Figures share traceable inputs and plot definitions. See the source records and reproduction notes:
[Input files](examples/) · [Dotplot / Fusion source](tools/site-flagships.R) · [Data sources and attribution](website/ASSET-SOURCES.md)

<a id="installation"></a>
## Installation

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", dependencies = NA)
```

The minimal example below needs only ggchord. Complete data-preparation examples declare optional dependencies such as dplyr, purrr and tidyr alongside their code.

<a id="quick-start"></a>
## Draw a first plot

Use bundled tables to draw sequence backbones, alignment ribbons and sequence names. No file download or alignment program is needed.

```r
library(ggchord)

p <- ggchord(seq_data_example, ribbon_data_example) +
  geom_seq() +
  geom_link_ribbon() +
  geom_seq_label() +
  coord_chord()

p
```

This small example shows the layer grammar. The figures above add their own source annotations, display choices and styling.

<a id="output"></a>
## Preview and save

Preview and save at the same explicit dimensions so labels are arranged for the intended figure size.

```r
view_ggchord(p, width = 10, height = 8, units = "in")
ggplot2::ggsave("ggchord.png", p, width = 10, height = 8, units = "in", dpi = 144)
```

<a id="documentation"></a>
## Documentation and resources

- [Plotting concepts](website/en/modes.md) · [Input tables](website/en/data.md) · [File formats](website/en/formats.md)
- [Labels](website/en/labels.md) · [Appearance](website/en/styling.md) · [Figure output](website/en/output.md)
- [Diagnostics](website/en/diagnostics.md) · [Function reference sources](man/) · [Changes](NEWS.md)

<details>
<summary>Design references and acknowledgements</summary>

ggchord draws on the following projects and public tools for design and comparison:

- [ggplot2](https://ggplot2.tidyverse.org/) — layered graphical grammar.
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — query/subject field naming.
- [Circos](https://circos.ca/) · [circlize](https://github.com/jokergoo/circlize) — circular hierarchy and tracks.
- [clinker](https://github.com/gamcil/clinker) · [gggenes](https://github.com/wilkox/gggenes) · [genoPlotR](https://genoplotr.r-forge.r-project.org/) — gene arrows and synteny comparisons.
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) · [gggenomes](https://thackl.github.io/gggenomes/) — annotation placement and local comparisons.
- [ggview](https://github.com/idmn/ggview) — previewing figures at export dimensions.

</details>

[Report an issue](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
