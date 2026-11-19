<!-- Generated from website/readme-template.md and website/readme-text.json. -->
<p align="center">
  <img src="man/figures/logo.png" alt="ggchord" width="280">
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord: sequences, annotations and relationships

Build layered sequence figures in R with ggplot2: compare multiple sequences, annotate a circular map, inspect HSP alignments and connect independent views.

> [!IMPORTANT]
> Use v0.29.0 for code that relies on the previous API, because v0.50.0 and later are incompatible with that API.

[Website](https://dangjem.github.io/ggchord/) · [Usage guide](website/en/modes.md) · [Figure gallery](website/en/gallery.md) · [Report an issue](https://github.com/DangJem/ggchord/issues)

<a id="modules"></a>
## Four plotting modules

| Module | Use it to | Views |
| --- | --- | --- |
| [I · Multi-sequence relationships](website/en/multi-sequence-views.md) | Compare sequence intervals, gene order and observed relationships. | Chord · Collinear · Radial · Sankey |
| [II · Single-sequence circular maps](website/en/circular.md) | Place features, primers and restriction sites on one sequence. | Circular |
| [III · HSP Dotplot](website/en/dotplot.md) | Inspect alignment intervals on query and subject axes. | Dotplot |
| [IV · Multi-view composition](website/en/fusion.md) | Connect independent plots through explicit source intervals. | Fusion |

The multi-sequence module has four Coord types. Circular and Dotplot add two more; Fusion is the composition layer. Layers, mappings, scales, positions and themes retain their ggplot2 roles.

<a id="examples"></a>
## Real-data figures

### I · Multi-sequence relationships

<table>
  <tr>
    <td width="50%"><strong>Chord · relationships across four genomes</strong><br><a href="man/figures/showcasechord.png"><img src="man/figures/showcasechord.png" alt="Four-genome Chord flagship with labelled genes and alignment ribbons" width="460"></a><br>Combine mixed-curvature backbones, strand-aware CDS arrows, selected source labels and Identity-coloured ribbons. Features and HSPs are display subsets of the original records.<br><a href="website/en/case-chord.md">Complete example →</a></td>
    <td width="50%"><strong>Collinear · order, orientation and cross-row links</strong><br><a href="man/figures/showcasecollinear.png"><img src="man/figures/showcasecollinear.png" alt="Four-genome Collinear flagship with strand-separated genes and cross-row ribbons" width="460"></a><br>Compare four genomes on a shared base-pair scale. Selected CDS sit on either side of each backbone; dashed outlines distinguish links across non-adjacent rows.<br><a href="website/en/case-collinear.md">Complete example →</a></td>
  </tr>
</table>

<table>
  <tr>
    <td width="50%"><strong>Radial · synteny along radial spokes</strong><br><a href="man/figures/flagradial.png"><img src="man/figures/flagradial.png" alt="Radial flagship derived from the four-genome Chord example" width="460"></a><br>Derived from the Chord flagship, using its source preparation and selected features. Radial backbones expose interval relationships while preserving sequence identities and native positions.<br><a href="website/en/case-radial.md">Complete example →</a></td>
    <td width="50%"><strong>Sankey · source intervals across two mirrored columns</strong><br><a href="man/figures/flagsankey.png"><img src="man/figures/flagsankey.png" alt="Sankey flagship derived from the four-genome Collinear example" width="460"></a><br>Derived from the Collinear flagship, retaining its selected genes and HSPs. Mirrored sequence columns show observed interval relationships; ribbon width is not a conserved flow quantity.<br><a href="website/en/case-sankey.md">Complete example →</a></td>
  </tr>
</table>

### II · Single-sequence circular maps

<a href="man/figures/showcaseplasmid.png"><img src="man/figures/showcaseplasmid.png" alt="pETDuet-1 circular flagship with features, primers and restriction sites" width="460"></a>

A complete circular map with 20 features, five reference primers and 58 grouped cut positions from the declared enzyme set. Original sequence positions and annotation identities are retained.

[Code for the dense map above](tools/documentation-examples.R) · [Compare the pETDuet-1 reference map](website/en/plasmid-catalog.md#catalogduet)

### III · HSP Dotplot

<a href="man/figures/flagdotplot.png"><img src="man/figures/flagdotplot.png" alt="Four-sequence matrix of 32 original HSPs" width="460"></a>

32 original HSPs; Identity by colour, reported alignment length by line width.

[Complete example](website/en/case-dotplot.md)

### IV · Multi-view composition

<a href="man/figures/flagfusion.png"><img src="man/figures/flagfusion.png" alt="Linked Chord, Collinear, Radial and local-detail views" width="460"></a>

Global relationships, ordered comparisons and a magnified source interval in one composition.

[Complete example](website/en/case-fusion.md)

Figures share traceable inputs and plot definitions. See the source records and reproduction notes:
[Input files](examples/) · [Dotplot / Fusion source](tools/site-flagships.R) · [Plot source](tools/documentation-examples.R) · [Data sources and attribution](website/ASSET-SOURCES.md)

<a id="get-started"></a>
## Get started

<a id="installation"></a>
### Installation

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

Install the stable release from [CRAN](https://CRAN.R-project.org/package=ggchord):

```r
install.packages("ggchord")
```

Install the v0.29.0 release from GitHub:

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", ref = "v0.29.0", dependencies = NA)
```

The minimal example below needs only ggchord. Complete data-preparation examples declare optional dependencies such as dplyr, purrr and tidyr alongside their code.

<a id="quick-start"></a>
### Draw a first plot

Use bundled tables to draw sequence backbones, alignment ribbons and sequence names. No file download or alignment program is needed.

```r
library(ggchord)

p <- ggchord(prepare_ggchord_data(seq_data_example, ribbon_data_example)) +
  geom_seq() +
  geom_link_ribbon() +
  geom_seq_label() +
  coord_chord()

p
```

This small example shows the layer grammar. The figures above add their own source annotations, display choices and styling.

<a id="output"></a>
### Preview and save

Preview and save at the same explicit dimensions so labels are arranged for the intended figure size.

```r
view_ggchord(p, width = 10, height = 8, units = "in")
ggplot2::ggsave("ggchord.png", p, width = 10, height = 8, units = "in", dpi = 144)
```

<a id="design-inspiration"></a>
## Design inspiration

The entries distinguish documented design influences from projects considered only in capability research. See the [design-source audit](tools/site-audit/design-inspiration.md) for evidence.

### Graphical grammar and sequence relationships

- [ggplot2](https://ggplot2.tidyverse.org/) — Provides the layered grammar and the division of responsibilities between mappings, scales, guides, geometry, coordinates and themes. ggchord follows this model so users can compose layers and supply their own scales and styles.
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — Informs the compact query/subject field names used to trace alignment records, endpoints and sequence identifiers. Plotting consumes existing results rather than running an alignment program.
- [Circos](https://circos.ca/) · [circlize](https://github.com/jokergoo/circlize) — Circos informs circular information hierarchy and ribbon presentation; circlize is a capability comparison for multi-track layouts and explicit track/ring contracts.

### Gene arrows, annotation placement and local comparison

- [clinker](https://github.com/gamcil/clinker) · [gggenes](https://github.com/wilkox/gggenes) · [genoPlotR](https://genoplotr.r-forge.r-project.org/) — Offer clear examples of strand-aware gene arrows, gene-cluster comparisons and linear synteny. They inform how sequence direction, feature placement and connections are presented together.
- [gggenomes](https://thackl.github.io/gggenomes/) — Informs tidy multi-table inputs, focus/sync workflows, feature stacking, local anchoring and local display reversal. Its approach helps keep data preparation distinct from display choices.
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) — Provides a reference for readable feature annotations, crowded-label handling and local feature layouts. It informs the emphasis on connecting labels clearly to their source intervals.

### Local magnification, aligned characters and figure preview

- [ggmagnify](https://github.com/hughjonesd/ggmagnify) — Provides a reference for connecting an overview to a magnified inset. ggchord develops this presentation idea around explicit source intervals and independently registered sequence views.
- [ggmsa](https://bioconductor.org/packages/ggmsa/) — Listed as a character-view research reference in v0.20.0. In v0.28.0, its residue cells and chemistry palettes were inspected directly and informed geom_residue() background sizing and default fill scales. Glyph rendering, exact source registration and user-scale priority remain independently implemented.
[ggseqalign](https://www.bioconductor.org/packages/release/bioc/html/ggseqalign.html), [ggseqlogo](https://omarwagih.github.io/ggseqlogo/), [Jalview](https://www.jalview.org/help/html/features/annotation.html) and [ESPript](https://espript.ibcp.fr/ESPript/ESPript/) are research comparisons for aligned characters, sequence logos, column annotations or static output. Their integrations, logo, consensus and conservation layers are not implemented here. Historical records do not establish that each influenced the current layout algorithms.
- [ggview](https://github.com/idmn/ggview) — Informs the workflow of previewing a static figure at its intended export size. ggchord uses view_ggchord() for sized preview and standard ggsave() for saving.

[Report an issue](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
