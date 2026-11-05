<div class="task-intro"><p>Give colours a clear meaning: identify categories or continuous quantities, then choose mappings, scales and guides. Website accents do not dictate data colours.</p></div>

<span id="prepare"></span>

## Build styles from data roles {#roles}

Use four chord backbones and one observed alignment. Map the ribbon fill to the query ID and assign that category one named colour. Sequence colours retain their default scale.

{{code:styling}}
{{figure:styling}}

| Role | Mapped data | Useful scale family |
| --- | --- | --- |
| `seq_colour`, `seq_ring` | Sequence/category | Sequence colour or ring manual scales |
| `ribbon_fill`, `ribbon_alpha`, `ribbon_colour`, `ribbon_linetype` | Alignment metadata or statistics | Manual, gradient/steps/identity and alpha scales |
| `feature_fill`, `feature_shape`, `gene_fill`, `region_fill` | Annotation category or strand | Role-specific manual/preset scales |
| `link_colour`, `link_alpha`, `link_linewidth`, `link_linetype` | Point-link metadata | Continuous or manual link scales |
| `track_value`, `track_fill`, `track_colour` | Values/categories in one named track | Named value domain and colour/fill scales |

## Fixed styles and user scale priority {#priority}
Put a constant style outside `aes()`. Map variable styles inside it and supply an
appropriate scale. An explicitly fixed geom value takes precedence; omitting that
argument keeps mapping active. Your scale takes priority over automatic defaults.
Repeated same-kind layers keep independent data/mapping/parameters. A geometric
vector named by `accver` does not mean a scale also accepts per-sequence lists.

## Make legends readable {#guides}
`guide_ggchord_legend()` controls discrete role legends; `guide_ggchord_colourbar()`
controls continuous colour guides. The colour/color spelling aliases are retained,
with `colour` preferred in this documentation. Use `show.legend` on a layer,
`guides()` for guide choices and `theme(legend.position = ...)` for placement.
A wide guide may need multiple rows or a wider export, not smaller biological data.

## Select a theme for the scene {#themes}
This page uses the default theme. Add `theme()` to adjust text, legends or background; for example, `theme(text = element_text(size = 12))` sets the base text size. `theme_ggchord_publication()` provides a publication style and `theme_ggchord_plasmid()` suits circular maps. Themes control non-data appearance without changing coordinates or mapped colours.

## Check your result {#check-result}
Links from the same query share a colour; continuous identity needs a continuous scale. Put fixed colours outside aes() and make each guide explain its mapping.

<div class="next-step">Continue: <a href="adjustments.html">Common adjustments and grammar →</a></div>

See [scale families](scale-families.html) for runnable manual/continuous/steps/identity, line type, opacity and theme examples.
