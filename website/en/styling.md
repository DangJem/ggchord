## Build styles from data roles {#roles}
Run the [shared setup](index.html#quick-start). The example declares categories
explicitly: ribbon fill maps query ID; sequence colour maps sequence ID; gene fill
maps strand. Manual scale names match those values exactly.

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
The example stores minimal, dark and plasmid alternatives in `theme_variants`.
Add the appropriate object to the matching plot. `theme_ggchord_publication()`
is used above; the plasmid preset belongs with the circular workflow. These
objects style non-data appearance. They do not choose a coordinate system or
replace data-driven fill scales. Use `theme(text = element_text(...))` for text
appearance; removed `base_size` constructor-style shortcuts have a
[migration route](migration.html).
