## 0.19.0: static panel composition

Existing single-view calls keep their defaults. Build independent plots, then
use `ggchord_fusion(center) + fusion_panel(detail, "detail") +
fusion_ribbon(links, qpanel="center", spanel="detail")`. Use ordinary q/s
columns; panel/instance selectors are needed only for composition routes.
Cross-panel aesthetics use ordinary `aes(fill=...)` and `scale_fill_*()`;
role-specific scales still belong to their own panel. Save with `ggsave()`;
inspect physical geometry with `export_ggchord_fusion()`, not `ggplot_build()`.
Fully visible same-unit intervals are required. This additive interface does
not rename or remove any previous API.

[English](MIGRATION.md) · [现代汉语](MIGRATION-Hans.md)

# Migrate to the ggchord candidate API

## Radial synteny in 0.18.0

Switch a supported sequence/feature/Track composition with
`coord_radial_synteny()`. Spokes keep proportional source lengths; use
`geom_seq(length_mode = "equal")` only for explicit normalization. Cyclic
adjacency includes an observed last/first relationship, with no invented
closure. Full interval fronts use ordered `[min-1,max]` boundaries in this new
view; existing coordinate conventions are unchanged. Position axes start at 1.
Radial inverse mapping requires an x/y table on the selected backbone.
Shared straight scale bars and branches, joined window rows and circular-only
annotations are explicitly rejected. See the paired view-expansion tutorial
for source data, physical Track bands, routing and the complete support matrix.


## HSP Dotplots in 0.18.0

Use the same sequence and ribbon tables with `geom_alignment_segment()` and
`coord_dotplot()`; no `geom_seq()` is needed. Keep `geom_link_line()` for point
links. Alignment styles use the existing `link_*` roles and `scale_link_*()`.
Use normalized endpoint directions (`read_paf()` already does this); a raw
ascending target interval plus a minus strand must first be normalized.
Unknown directions are reported and omitted. There is no base-level path
without explicit alignment trace data. Select `axis = "x"` or `"y"` for
`map_ggchord_position()`; source coordinates after focusing use `source = TRUE`.
Default window-axis labels remain relative to the selected window. These new
APIs do not require changes to existing Chord, Circular or Collinear code.


## Unified annotations in 0.17.0

Prefer `ggchord(seq_data = seqs, feature_data = annotations)` for genes and
other features. `geom_feature()` inherits the full table. Without explicit
layer data, gene layers select `gene`/`CDS` from `feature_class`, `type`, or
`role`, in that order; tRNA/rRNA use generic feature layers. Missing strand is
unknown, and generic layers also accept bidirectional annotations.

`gene_data` remains accepted, including the third positional argument. It emits
one message when constructing a plot, with no repeated rendering messages and
no scheduled removal. Untyped legacy tables retain their gene semantics. To
migrate an untyped directional gene table, add `type = "gene"` before passing
it as `feature_data`, or supply it explicitly to `geom_gene(data = genes)`.
Do not pass both constructor inputs. Layer data overrides the global table;
combine disjoint gene/feature subsets rather than drawing the full table twice.

```r
library(ggchord)
library(dplyr)
data(seq_data_example)
data(gene_data_example)
annotations <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand, type = "CDS")
p <- ggchord(seq_data_example, feature_data = annotations) +
  geom_seq() +
  geom_gene()
```

`export_ggchord_layout()` exposes annotation identities in `layer_sources`;
the old geometry `source_row` still refers to the resolved layer input.
Constructor rows share source IDs across default annotation layers. Explicit
layer tables and functions have layer-local scope, without content-based
inference or removal of identical source observations.

## Track values and sequence axes in 0.17.0

`geom_track_bar()` now counts observations, matching native bar semantics;
`geom_track_col()` displays supplied numeric heights with `stat = "identity"`.
Neither function is removed. Migrate numerical bars by replacing the function
name with `geom_track_col()`, retaining the data, value mapping and physical band.

| Previous use | Replacement | Change |
| --- | --- | --- |
| `geom_track_bar(aes(track_value = gc), data = bins, ...)` for supplied values | `geom_track_col(aes(track_value = gc), data = bins, ...)` | Numeric-height behavior moves to col; bar counts rows. |
| `geom_track(..., geom = "bar")` for supplied values | `geom_track(..., geom = "col")` | Use bar with `aes(x = position)` for observation counts. |
| Explicit zero-based sequence tick vectors in examples | `scale_seq_position_continuous(n.breaks = 3)` | Automatic sequence axes start at 1 and prioritize real endpoints. |

`track = position_track(...)` names the physical band; `position` adjusts data.
The earlier `position = position_track(...)` band spelling remains accepted.
Source-coordinate windows retain their original genomic interval. Explicit
breaks and labels remain available; a closed circular axis combines 1/length
at the seam. Automatically managed numeric pident mappings use the 0–100
domain and ticks 0/50/80/90/95/100. Explicit user scales retain priority.

## Cross-row ribbon outlines in 0.16.7

`geom_link_ribbon()` uses `nonadjacent_linewidth = 0.4` for non-adjacent
collinear rows, independently of ordinary fixed `linewidth`. Set
`nonadjacent_linewidth = 0` to make those ribbons inherit the layer
linewidth. An explicit linewidth mapping controls the mapped widths.
The parameter changes width without changing the automatic dashed line type.

## Circular origin in 0.16.7

`coord_circular()` keeps `rotation = 0`, placing the origin at 12 o'clock.
Rotation is counterclockwise from this baseline: 90 places the origin at
9 o'clock, and -90 at 3 o'clock. To preserve an explicit angle from v0.16.6
(or the earlier 0.16.7 development snapshot), subtract 90 degrees. Calls that
omit rotation retain the current development appearance. Increasing coordinates
remain clockwise by default. The chord rotation contract stays unchanged.

## Scope and timing {#scope}
These are current candidate replacements, not a v1.0.0 compatibility promise.
The 0.16.0–0.16.4 phases were delivered together in v0.16.5. The current installation tag is v0.16.6.
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


```r
library(dplyr)
seqs <- seq_data_example |>
  select(accver, length)
links <- ribbon_data_example |>
  select(qaccver, saccver, qstart, qend, sstart, send)
# The fixture preserves biological strands separately from its old display strands.
genes <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand)
```


```r
library(dplyr)
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19 |>
  select(accver, length)
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source |>
  select(accver, start, end, strand, type, anno)
```



```r
p_migrated <- ggchord(plasmid_seq) +
  geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) +
  coord_circular()
p_offsets <- ggchord(seqs, gene_data = genes) +
  geom_seq() +
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

## Default chord rotation in v0.16.7

Chord plots now default to 0 degrees. Use `coord_chord(rotation = 45)` to request the previous 45-degree layout explicitly. Native protein plots declare `sequence_unit = "aa"` and use complete protein lengths.

## Foundation in v0.17.0

`geom_track_bar()` now counts observations at source positions (`stat="count"`).
For previously supplied heights, replace it with `geom_track_col()`; retain the
same mapping and `position = position_track("signal")`, or use
`track = position_track("signal")` with a separate native data `position`.
`geom_track()` supports point/line/path/step/col/bar/tile/rect/area/ribbon/segment/
linerange/errorbar/text and identity/count/bin/summary; bin needs explicit
`binwidth` or `breaks`. Arbitrary Geom/Stat extensions remain outside this version.

Feature workflows accept `feature_data` in prepare/clean/focus/segment; the result
uses `feature_data`, while old gene calls retain `gene_data`. Never supply both.
Source rows, unknown strands and metadata remain available.

`geom_seq(seq_order=c(A_high="A", B="B", A_low="A"))` creates display copies.
Filter real HSP rows on each layer and specify `between=c("A_high","B")` or
`between=c(2,3)`. A route must be unique. Supply stable `source_hsp_id` values to
check duplicate display across layers; `allow_repeated=TRUE` explicitly opts in.
Without that column identity is layer-local, never guessed from equal values.
The same plot can use Chord or Collinear. Sankey owns query-left/subject-right
copies; it does not interpret `between` as a multistage flow.

Use `coord_collinear(seq_axis="vertical")` for vertical backbones. Inverse mapping
uses the final y coordinate. `coord_sankey()` filters entity self links by default
and retains reciprocal input rows; shared branches are outside its first MVP.
Use exported `display_instances` to select an unambiguous instance for mapping.

`geom_seq(length_mode="equal")` changes only visual allocation. Source lengths,
coordinates and local axes remain scientific coordinates; a shared absolute
scale bar is rejected. Default ticks now label true position 1 and the final
position. A closed circular seam combines them; explicit breaks take priority.

Ribbons retain the previous Identity fill and full interval fronts by default.
Use `gradient = TRUE` to opt into endpoint sequence-colour interpolation.
Use `tip="subject"` for an optional decorative taper; its shoulder retains the
translated target endpoint span and its tip is not a single matched position.
Use `tip="none"` for untapered fronts. Map `aes(ribbon_fill=pident)` to select the
Identity mode; it uses the fixed 0–100 domain and proportional intervals
0/50/80/90/95/100. Fixed fill, user mappings and user scales override the gradient.

Gradient fills use a device-resolved raster mesh along the true curve; outlines remain vector paths. This avoids translucent strip seams and preserves curved source fronts.

Automatic Identity guides treat generic legend key sizes as swatch dimensions, rather than compressing the entire numerical ruler. Their length accounts for boundary-label measurements and interval proportions. Set a Ribbon-specific key dimension or an explicit guide to override it. Circular parallel axis text follows the arc glyph by glyph; use horizontal or numeric orientation for straight text.
