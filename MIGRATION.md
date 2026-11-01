[English](MIGRATION.md) · [现代汉语](MIGRATION-Hans.md)

# Migrate to the ggchord candidate API

## Scope and timing {#scope}
These are current candidate replacements, not a v1.0.0 compatibility promise.
The 0.16.0–0.16.4 phases were delivered together in v0.16.5. The current
0.16.6.9000 checkout is unreleased; the formal installation tag remains v0.16.5.
The RC audit retains flexible parameter formats, spelling aliases, reexports and
getter/export entry points without a new deprecation. Old calls in this
table explain migration and must not be executed as current examples. Exact
first-removal tags are not guessed where historical evidence gives only an upper bound.

## Replacement table {#table}
| Old (do not run) | Current | Removal / guard | Meaning |
| --- | --- | --- | --- |
| `geom_ribbon()` | `geom_link_ribbon()` | ≤ 0.16.1 | Renamed link layer; not ggplot2::geom_ribbon(). |
| `ggchord(..., rotation = 20)` | `ggchord(...) + coord_chord(rotation = 20)` | ≤ 0.16.1 | Global coordinates are composed with +. |
| `geom_feature(feature_shape = "block", feature_width = .05)` | `geom_feature(shape = "block", width = .05)` | 0.16.1 | Mapped feature_shape is retained. |
| `geom_gene(gene_offset = c("+" = .1, "-" = .2))` | `geom_gene(position = position_strand(c("+" = .1, "-" = -.2)))` | 0.16.1 | Explicit Position values are signed. |
| `geom_feature(feature_offset = .1)` | `geom_feature(position = position_strand(.1))` | 0.16.1 | Use a Position object. |
| `geom_gene(gene_width = .05)` | `geom_gene(width = .05)` | 0.16.1 | Fixed geometry uses the shared width name. |
| `geom_link_ribbon(ribbon_gap = .01)` | `geom_link_ribbon(link_gap = .01)` | 0.16.1 | Use link_gap; likewise use link_ctrl_point. |
| `geom_gene_label_repel(gene_label_layout = "aligned")` | `geom_gene_label_repel(label_layout = "auto")` | 0.16.1 | Migrate both the name and value. |
| `geom_gene_label_repel(label_rotation = 90)` | `geom_gene_label(label_rotation = 90)` | 0.16.4 guard | Previously ignored manual control now errors; alternatively remove it. |
| `geom_feature_label_repel(label_orientation = "horizontal")` | `geom_feature_label(label_orientation = "horizontal")` | 0.16.4 guard | Use the fixed feature-label layer for unknown-strand features. |
| `geom_restriction_site(leader = "trunk")` | `geom_restriction_site(leader = "radial")` | 0.16.1 | The radial default remains. |
| `theme_ggchord(base_size = 10)` | `theme_ggchord() + theme(text = element_text(size = 10))` | ≤ 0.16.1 | Text appearance belongs to theme elements. |

## Run the replacement {#replacement}
In a fresh session, run this preparation before the current replacements:

<details class="preparation"><summary>Preparation code — run once on this page</summary>

```r
library(ggchord)
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)
seqs <- seq_data_example[c("accver", "length")]
links <- ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
# The fixture preserves biological strands separately from its old display strands.
genes <- transform(gene_data_example, strand = source_strand)[c("accver", "start", "end", "strand")]

data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19[c("accver", "length")]
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source[c("accver", "start", "end", "strand", "type", "anno")]
```

</details>


```r
p_migrated <- ggchord(plasmid_seq) + geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) + coord_circular()
p_offsets <- ggchord(seqs, gene_data = genes) + geom_seq() +
  geom_gene(position = position_strand(c("+" = .1, "-" = -.2)))
```


Old minus-strand offset magnitudes need a negative sign in explicit Position
values. Mapping still uses `aes(feature_shape = type)` with a shape scale;
the fixed parameter is `shape`. Automatic labels use `label_layout`; manual
geometry uses the corresponding fixed-label layer. `label_layout = "aligned"`
becomes `"auto"`. Removed names are errors, not silently retained execution paths.

## What remains supported {#retained}
The five colour/color spelling aliases, all upstream reexports, gene/feature/plasmid
shortcuts and both geometry getter/export interfaces remain. Coordinate-specific
defaults and already documented result fields are retained. See the
API contract (`help("ggchord_api_contract", package = "ggchord")`) and current reference for
NULL/NA/missing semantics; NULL does not universally mean disabled.
