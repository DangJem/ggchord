<!-- Generated from website/readme-template.md and website/readme-text.json. -->
<p align="center">
  <img src="man/figures/logo.png" alt="ggchord" width="280">
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# {{title}}

{{intro}}

[{{website}}](https://dangjem.github.io/ggchord/) · [{{manual}}](website/{{lang}}/modes.md) · [{{gallery}}](website/{{lang}}/gallery.md) · [{{issues}}](https://github.com/DangJem/ggchord/issues)

<a id="modules"></a>
## {{modules}}

| {{module}} | {{purpose}} | {{view}} |
| --- | --- | --- |
| [{{multi}}](website/{{lang}}/multi-sequence.md) | {{multi-purpose}} | Chord · Collinear · Radial · Sankey |
| [{{circular}}](website/{{lang}}/circular.md) | {{circular-purpose}} | Circular |
| [{{dotplot}}](website/{{lang}}/dotplot.md) | {{dotplot-purpose}} | Dotplot |
| [{{fusion}}](website/{{lang}}/fusion.md) | {{fusion-purpose}} | Fusion |

{{grammar}}

<a id="examples"></a>
## {{examples}}

### {{chord-title}}

{{chord-caption}}

![{{chord-alt}}](man/figures/showcasechord.png)

[{{case}}](website/{{lang}}/chord.md) · [{{four-views}}](website/{{lang}}/multi-sequence-views.md) · [{{plot-source}}](tools/documentation-examples.R)

### {{linear-title}}

{{linear-caption}}

![{{linear-alt}}](man/figures/showcasecollinear.png)

[{{case}}](website/{{lang}}/collinear.md) · [{{layout}}](website/{{lang}}/collinear-layout.md) · [{{plot-source}}](tools/documentation-examples.R)

### {{radial-title}}

{{radial-caption}}

![{{radial-alt}}](man/figures/flagradial.png)

[{{case}}](website/{{lang}}/radial-synteny.md) · [{{plot-source}}](tools/site-flagships.R)

### {{sankey-title}}

{{sankey-caption}}

![{{sankey-alt}}](man/figures/flagsankey.png)

[{{case}}](website/{{lang}}/sankey.md) · [{{plot-source}}](tools/site-flagships.R)

### {{circular-title}}

{{circular-caption}}

![{{circular-alt}}](man/figures/showcaseplasmid.png)

[{{case}}](website/{{lang}}/circular-dense-case.md) · [{{reference-maps}}](website/{{lang}}/plasmid-catalog.md) · [{{plot-source}}](tools/documentation-examples.R)

### {{more-views}}

<table>
  <tr>
    <th>{{dotplot}}</th>
    <th>{{fusion}}</th>
  </tr>
  <tr>
    <td width="50%"><a href="man/figures/flagdotplot.png"><img src="man/figures/flagdotplot.png" alt="{{dotplot-alt}}" width="460"></a></td>
    <td width="50%"><a href="man/figures/flagfusion.png"><img src="man/figures/flagfusion.png" alt="{{fusion-alt}}" width="460"></a></td>
  </tr>
  <tr>
    <td>{{dotplot-caption}}</td>
    <td>{{fusion-caption}}</td>
  </tr>
  <tr>
    <td><a href="website/{{lang}}/dotplot-start.md">{{dotplot-guide}}</a></td>
    <td><a href="website/{{lang}}/fusion-details.md">{{fusion-guide}}</a></td>
  </tr>
</table>

{{source-note}}
[{{source-inputs}}](examples/) · [{{source-figures}}](tools/site-flagships.R) · [{{source-attribution}}](website/ASSET-SOURCES.md)

<a id="get-started"></a>
## {{get-started}}

<a id="installation"></a>
### {{installation}}

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord", dependencies = NA)
```

{{optional-deps}}

<a id="quick-start"></a>
### {{quick-start}}

{{quick-intro}}

```r
library(ggchord)

p <- ggchord(seq_data_example, ribbon_data_example) +
  geom_seq() +
  geom_link_ribbon() +
  geom_seq_label() +
  coord_chord()

p
```

{{quick-note}}

<a id="output"></a>
### {{output}}

{{output-intro}}

```r
view_ggchord(p, width = 10, height = 8, units = "in")
ggplot2::ggsave("ggchord.png", p, width = 10, height = 8, units = "in", dpi = 144)
```

<a id="design-inspiration"></a>
## {{inspiration-title}}

{{inspiration-intro}}

### {{inspiration-grammar-title}}

- [ggplot2](https://ggplot2.tidyverse.org/) — {{inspiration-grammar}}
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) — {{inspiration-naming}}
- [Circos](https://circos.ca/) · [circlize](https://github.com/jokergoo/circlize) — {{inspiration-circular}}

### {{inspiration-features-title}}

- [clinker](https://github.com/gamcil/clinker) · [gggenes](https://github.com/wilkox/gggenes) · [genoPlotR](https://genoplotr.r-forge.r-project.org/) — {{inspiration-synteny}}
- [gggenomes](https://thackl.github.io/gggenomes/) — {{inspiration-genomes}}
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) — {{inspiration-features}}

### {{inspiration-detail-title}}

- [ggmagnify](https://github.com/hughjonesd/ggmagnify) — {{inspiration-magnify}}
- [ggmsa](https://bioconductor.org/packages/ggmsa/) · [ggseqalign](https://bioconductor.org/packages/ggseqalign/) — {{inspiration-alignment}}
- [ggseqlogo](https://omarwagih.github.io/ggseqlogo/) · [Jalview](https://www.jalview.org/help/html/features/annotation.html) — {{inspiration-columns}}
- [ggview](https://github.com/idmn/ggview) — {{inspiration-preview}}

[{{issues}}](https://github.com/DangJem/ggchord/issues) · [MIT {{license}}](LICENSE)
