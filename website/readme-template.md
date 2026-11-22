<!-- Generated from website/readme-template.md and website/readme-text.json. -->
<p align="center">
  <img src="man/figures/logo.png" alt="ggchord" width="280">
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# {{title}}

{{intro}}

> [!IMPORTANT]
> {{development-status}}

[{{website}}](https://dangjem.github.io/ggchord/) · [{{manual}}]({{tutorial-path}}) · [{{issues}}](https://github.com/DangJem/ggchord/issues)

<a id="modules"></a>
## {{modules}}

| {{module}} | {{purpose}} | {{view}} |
| --- | --- | --- |
| {{multi}} | {{multi-purpose}} | Chord · Collinear · Radial · Sankey |
| {{circular}} | {{circular-purpose}} | Circular |
| {{dotplot}} | {{dotplot-purpose}} | Dotplot |
| {{fusion}} | {{fusion-purpose}} | Fusion |

{{grammar}}

<a id="examples"></a>
## {{examples}}

### {{multi}}

<table>
  <tr>
    <td width="50%"><strong>{{chord-title}}</strong><br><a href="man/figures/showcasechord.png"><img src="man/figures/showcasechord.png" alt="{{chord-alt}}" width="460"></a><br>{{chord-caption}}<br><a href="examples/make_readme_plots.R">{{case}} →</a></td>
    <td width="50%"><strong>{{linear-title}}</strong><br><a href="man/figures/showcasecollinear.png"><img src="man/figures/showcasecollinear.png" alt="{{linear-alt}}" width="460"></a><br>{{linear-caption}}<br><a href="examples/make_readme_plots.R">{{case}} →</a></td>
  </tr>
</table>

<table>
  <tr>
    <td width="50%"><strong>{{radial-title}}</strong><br><a href="man/figures/flagradial.png"><img src="man/figures/flagradial.png" alt="{{radial-alt}}" width="460"></a><br>{{radial-caption}}<br><a href="examples/make_readme_plots.R">{{case}} →</a></td>
    <td width="50%"><strong>{{sankey-title}}</strong><br><a href="man/figures/flagsankey.png"><img src="man/figures/flagsankey.png" alt="{{sankey-alt}}" width="460"></a><br>{{sankey-caption}}<br><a href="examples/make_readme_plots.R">{{case}} →</a></td>
  </tr>
</table>

### {{circular}}

<a href="man/figures/showcaseplasmid.png"><img src="man/figures/showcaseplasmid.png" alt="{{circular-alt}}" width="460"></a>

{{circular-caption}}

[{{case}}](examples/make_readme_plots.R)

### {{dotplot}}

<a href="man/figures/flagdotplot.png"><img src="man/figures/flagdotplot.png" alt="{{dotplot-alt}}" width="460"></a>

{{dotplot-caption}}

[{{case}}](examples/make_readme_plots.R)

### {{fusion}}

<a href="man/figures/flagfusion.png"><img src="man/figures/flagfusion.png" alt="{{fusion-alt}}" width="460"></a>

{{fusion-caption}}

[{{case}}](examples/make_readme_plots.R)

{{source-note}}
[{{source-inputs}}](examples/) · [{{plot-source}}](examples/make_readme_plots.R) · [{{source-figures}}](data-raw/) · [{{source-attribution}}](website/ASSET-SOURCES.md)

<a id="get-started"></a>
## {{get-started}}

<a id="installation"></a>
### {{installation}}

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

{{install-cran}}

```r
install.packages("ggchord")
```

{{install-development}}

```r
install.packages("remotes")
remotes::install_local(".", dependencies = NA)
```

{{optional-deps}}

<a id="quick-start"></a>
### {{quick-start}}

{{quick-intro}}

```r
library(ggchord)
library(ggplot2)

p <- ggchord(seq_data_example) +
  geom_seq() +
  geom_link_ribbon(aes(fill = pident), data = ribbon_data_example) +
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
ggsave("ggchord.png", p, width = 10, height = 8, units = "in", dpi = 144)
```

<a id="design-inspiration"></a>
{{design-references}}

[{{issues}}](https://github.com/DangJem/ggchord/issues) · [MIT {{license}}](LICENSE)
