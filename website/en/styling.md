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

## Target-size layout {#adaptive}

Add `adapt_ggchord()` to opt into target-size layout. The bundled example uses
recorded sequence lengths and HSPs; the reference canvas is in inches and the
automatic text floor is in mm. The report describes the requested 8 × 6 inch
PNG; use the same size with standard `ggsave()` to reproduce that arrangement.

```r
library(ggchord)
library(dplyr)
seqs <- seq_data_example |>
  select(accver, length)
plot <- ggchord(seqs, ribbon_data_example) +
  geom_seq() +
  geom_link_ribbon() +
  adapt_ggchord(
    mode = "adaptive",
    reference_size = c(8, 6),
    text_scale_range = c(0.85, 1.35),
    min_text_size = 2.5
  )
report <- diagnose_ggchord(plot, width = 8, height = 6, device = "png")
report$guide_layout$guides
report$adaptive_layout$styles
```

`fixed` holds default theme and guide furniture at reference physical sizes.
`scale` adjusts theme-driven defaults within the chosen multiplier range.
`adaptive` adds a bounded joint search for native guide grids, sides and
directions, then remeasures the result. It reuses existing label layout;
`feasible` is a guide-space result, not a promise that every label is collision-free.
No policy preserves the existing plotting behaviour.

Explicit geom styles, theme sizes, guide positions/directions/grids/key sizes,
scales and `show.legend` remain authoritative. Continuous/stepped bars remain
intact; categories and source records are retained. When constraints exceed the
canvas budget, diagnostics report `insufficient_space` and rendering warns.
Increase the export size or revise the explicit guide settings. The bounded
search does not guarantee a global optimum. Layout-only exports record the
policy but require a target device for final guide measurements.

Apply the policy to each desired Fusion panel before composing it. Each panel
uses its allocated size; outer panel allocations and the cross-panel guide are
separate. Facets, arbitrary custom grobs, blanket scaling of geometric constants
and font embedding certification are outside this policy.

## Check your result {#check-result}
Links from the same query share a colour; continuous identity needs a continuous scale. Put fixed colours outside aes() and make each guide explain its mapping.

<div class="next-step">Continue: <a href="adjustments.html">Common adjustments and grammar →</a></div>

See [scale families](scale-families.html) for runnable manual/continuous/steps/identity, line type, opacity and theme examples.
