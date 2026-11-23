# API version choice

Use v0.29.0 for code that relies on the previous API, because v0.50.0 and later are incompatible with that API.

<details>
<summary>Internal change records for maintainers</summary>

These archived records document removals, replacement calls and verified behaviours for the versions named in each section. Maintainers use them to trace compatibility decisions and reproduce earlier audits.

## Migrating to 0.50.0

The table covers removed calls and changes in behaviour in 0.50.0. Later sections identify the versions to which historical examples apply.

| Former call or argument | 0.50.0 replacement |
| --- | --- |
| Layers shared scales only with `id = NULL` or `id = NA` | Ordinary visual scales now share the aesthetic by default. Scope a scale with `scale_for(scale_fill_discrete(), layers = "features")`. |
| A label layer relied on independent plot data with `anchor_layer` omitted | Set `anchor_layer = NULL` for independent labels. Omit it to follow a single compatible target. With multiple candidates, give the intended layer an ID and use that name as `anchor_layer`. |
| Feature colour previously followed the remapped accession in `aes(accver = ref, fill = accver)` | Each starting mapping reads the layer's source columns. Keep this expression to colour by the original accession. Use `aes(accver = ref, fill = ref)` to colour by the displayed reference. |
| A point track uses categories for shape or needs to hide only its fill guide | Map the categories with `aes(shape = kind)` and set symbols with `scale_shape_manual(values = c(a = 21, b = 24))`. Set `show.legend = c(fill = FALSE)` on the track layer to hide its contribution to the fill guide. Other layers can still contribute to that guide. |
| Mixing coverage fraction and mean depth on a common height axis | Give them separate named bands. Declare custom heights with `measure = "read_depth"`. To share a range across compatible plots, apply the output of `shared_track_limits(p1, p2)` to each track scale. |
| Individual track ranges in sample panels | Use `scale_for(scale_y_continuous(), tracks = "depth", by = "panel")`. |
| Colour whole coverage windows | Use `geom_track_tile(aes(fill = after_stat(coverage)), stat = "coverage", binwidth = 4, data = intervals)`. The final tile uses its actual window length. |
| Assign pairwise observations to panels | Facet original ribbons or point links by a source column, e.g. `facet_wrap(vars(saccver))`. On a Dotplot specify both endpoints, e.g. `facet_grid(rows = vars(qaccver), cols = vars(saccver))`. |
| Repeated keys or two different focused intervals called W | Keep one row per reference ID. Give different intervals distinct `locus_id` values. Use `layout_sequences(order = c(left = "A", right = "A"))` for two displays of one reference. |
| A sample key called `count` in a count track | Copy it to `sample` before using `facet_wrap(vars(sample))`. Computed columns keep their statistical meaning. |
| Source-stage exports that also used default geometry | Pass `include = "track"` alongside `data_stage = "source"` for track geometry, or request both stages. |
| Region highlights specified as display edges | Supply one-based inclusive endpoints. A single-position highlight now occupies one sequence unit. Restriction site coordinates must be in bp on a known reference. |
| Overlaid `geom_seq()` layers give different lengths or units for the same accession | Use one reference table for the layers and set styles separately. Draw separate plots for incompatible references. |
| Track, region or point-link data declared in bp on an aa reference, or vice versa | Supply coordinates in the backbone unit. For DNA intervals on a protein reference, project through the relevant CDS or draw separate plots on the appropriate references. |
| A track expression needs the newly mapped accession, or a filtered input lacks a mapped column | Compute the new accession in a source column first and refer to that column in each mapping. For empty layers, filter a complete input table to zero rows to retain its columns and types. |
| `geom_track_col(aes(fill = kind), fill = "red", stat = "coverage", binwidth = 10)` | Remove the fill mapping to pool intervals. For separate estimates in the same colour, use `aes(group = kind)` and `fill = "red"`. |
| Relying on identical row contents to attach labels to a separately created feature table | Pass the same table to both layers, or supply matching `feature_id` and source intervals. Reordered subsets use the key. |
| Unrecognised arguments to `geom_*()`, such as `geom_seq(linewditth = 2)` | Use `geom_seq(linewidth = 2)`. Starting in 0.50.0, unknown names cause an error. Supply Stat-specific parameters to the selected Stat, for example `geom_track_line(stat = "coverage", binwidth = 4)` with interval data. |
| Partial argument names such as `position_offset(offse = .1)` or `guide_seq_axis(sequence = "A")` | Use `position_offset(offset = .1)` and `guide_seq_axis(sequences = "A")`. Spell out names in refactored calls, including arguments forwarded through `...`. Positional calls such as `position_offset(.1)` still work. |
| One explicit `width` for an external Stat's neighbourhood and the column Geom | Name the Stat argument `radius` and pass both controls, for example `geom_track_col(stat = MyStat, radius = 4, width = 2)`. |
| `ggchord(seq_data = ref, gene_data = f, ribbon_data = h)` | Start with `ggchord(ref)`, then add `geom_feature(data = f, preset = "gene")` and `geom_link_ribbon(data = h)`. For a classified primary table use `ggchord(prepare_ggchord_data(ref, ribbon_data = h, feature_data = f))`. When several prepared inputs are eligible, specify the layer table. |
| `gene_data` in prepare, validate, clean, focus or segment, and `mapping$gene` in preparation | In 0.50.0 the third argument is `feature_data`. Use `prepare_ggchord_data(ref, feature_data = f, mapping = list(feature = c(strand = "source_strand")))`. The primary annotation output, conversion report and validation count use `$feature_data`, `$report$conversions$feature` and `$data_summary$n_feature` respectively. Named tables remain independent in `features = list(cds = cds, sites = sites)`. |
| A saved prepared object with a populated `$gene_data`, or an unclassified third positional argument used with a gene preset | Rebuild from the original inputs. An inherited gene preset selects gene/CDS rows by the first available classification column (`feature_class`, `type`, `role`). If the data source identifies all rows as genes, use `feature_data = dplyr::mutate(f, type = "gene")` on an unclassified table. For a directional gene table without a classification column, use `ggchord(ref) + geom_seq() + geom_feature(data = f, preset = "gene")`. Inherited annotation row IDs use the `feature:` prefix. An explicit layer has layer-local row IDs. |
| `track_band("depth", thickness = w)` | Use `layout_track("depth", thickness = w)` and `track = "depth"` on its layers. For an implicit track, add `layout_track(thickness = w)` after its layer. |
| `seq_order`, `seq_orientation`, `seq_gap`, `length_mode` on `geom_seq()` | Add `layout_sequences(order = ids, orientation = directions, gap = gaps, length_mode = "equal")`. The geom still sets the backbone's appearance. |
| `position_plasmid(offset)` and `position_feature_stack(side = ...)` | Use `position_offset(offset)`. For old `side = "outside"`, use `base_position = position_offset(-0.1), direction = "inward"`. For old `"inside"`, use `position_offset(0.1)` and `"outward"`. The former `"strand"` setting is the default. |
| `geom_gene()`, `geom_feature_plasmid()` | Use `geom_feature(preset = "gene")` or `geom_feature(preset = "plasmid")`, with the original layer data. |
| `geom_feature_label_repel()` | Use `geom_feature_label(placement = "repel")`. Move the old placement policy to `label_layout` and the old `label_fit` to `overflow`. Supply `anchor_layer = "features"` when several feature layers could own the labels. |
| `scale_seq_position()` | Use `scale_seq_continuous()`. Select reference axes with `guide = guide_seq_axis(sequences = ids)`. |
| `scale_track_value("depth", limits = c(0, 1000))` | Add `scale_y_continuous(limits = c(0, 1000))` after the track, or use `scale_for(scale_y_continuous(limits = c(0, 1000)), tracks = "depth")`. |
| Traditional default colours or theme elements supplied to a preset | Add `scales_ggchord()` for the palettes and `theme(...)` for theme elements. |

When inspecting a staged export, check `layer_stage_status$status` before using
track statistics. A successful zero-row statistic is `available`. A failed Stat
sets `status = "failed"` for its statistic and geometry. Read `reason` for the
cause, including when some panels still return rows.
An ordinary build stops on a track Stat error. Use `get_chord_layout(p, build = FALSE)`
to inspect the last successful layout, then rebuild after correcting the input.
For zero input rows, built-in computed fields such as `count` and `prop` remain available for selection.

The [audit example](examples/migrate-v0.50.0.R) records a two-sample comparison and the contributing rows of each window. Each panel uses the same 10 bp reference and numerical scale. A second plot labels selected features in the lanes assigned to the full input.

## 0.29.0 data and mappings

| Workflow | Current usage |
| --- | --- |
| Standalone cleaning and reverse hits | Default sorting preserves relative direction. Use `reversed_interval = "keep"` on the original table when directional endpoints are required. |
| Custom BLAST fields | Untyped custom metadata retains text. Use `colClasses = c(score_custom = "numeric")` to read a measurement, or convert it afterward with `mutate(score_custom = as.numeric(score_custom))`. Other custom identifiers keep their spelling. |
| Track accessions | `aes(accver = id)` must return a scalar or one accession per row. Shorter vectors that were previously recycled now raise an error. |
| Scientific examples | Use `strand = source_strand` for biological interpretation and `unstranded = "keep"` to retain unknown GFF3 strands. |
| Feature direction | When supplying both `strand` and `directionality` (or its `direction` alias), correct the conflicting rows at the source or pass just the authoritative field to the layer. |
| Layer data functions | Each layer reads its function once per fresh render. Orientation inference, geometry and the exported layer input use the same rows. |
| Prepared data | Pass the default tables together as `ggchord(prepared)`. Select a named collection at its layer, for example `geom_feature(data = prepared$features$cds)`. Rename columns in `feature_data` with `mapping = list(feature = c(strand = "source_strand"))`. |
| Quantitative tracks | Use `track = position_track("gc")` for the physical band on line and tile layers. Ordinary lines accept `aes(x = pos, y = gc)`. Set `closed = TRUE` for circular closure. Interval tiles keep their start and end positions. Choose `geom_track_bar()` for counts or `geom_track_col()` for supplied heights. |
| Automatic labels | Set `placement = "callout", overflow = "ellipsis"` on `geom_feature_label_repel()` for external labels with truncation, or `placement = "feature", overflow = "shrink"` to fit text to a feature. Existing `label_layout` and `label_fit` calls retain their interpretation. Use one pair per call. |
| Inspect layer data | `export_ggchord_layout(p, data_stage = c("source", "mapped", "statistic", "geometry"))$layer_data` returns the selected stages by layer ID. Source is the adapted input. Native track mappings execute inside ggplot2, so their mapped entry is NULL and their statistic is available. A data callback runs afresh for this export. |
| Themes and guides | Use `theme(text = element_text(family = "sans"))` for plot-wide text, `labs(ribbon_fill = "Identity")` for a single ribbon legend title, and `theme(ggchord.seq.label = element_text(size = 10))` for sequence label size. |
| Shape a ribbon | Start with `geom_link_ribbon()`. Set `link_gap` for endpoint clearance, `link_avoid = "smooth"` for local avoidance, and `link_ctrl_point` when the middle of the curve needs explicit controls. |
| Align the origin | For twelve, three, six and nine o’clock, use 0, −90, 180 and 90 degrees respectively in `coord_circular(rotation = ...)`, or 90, 0, 270 and 180 in `coord_radial_synteny(start_angle = ...)`. Circular `direction` sets the direction of increasing positions. Radial `clockwise` orders subsequent spokes. |
| Ribbon statistics | `geom_link_ribbon(stat = "ribbon_bundle", bins = 20, weight = "count")` is equivalent to `stat_ribbon_bundle(bins = 20, weight = "count")`. To keep each link and map local density to alpha, use `geom_link_ribbon(aes(ribbon_alpha = after_stat(density)), stat = "ribbon_density", bins = 20, weight = "count")`. Only bundling accepts `min_bundle`. |
| Feature and track positions | Use `geom_feature(position = "feature_stack")` for default lanes or `position_feature_stack(spacing = 0.15)` to set their spacing. A quantitative layer declares its physical band with `track = position_track("depth")` and adjusts source coordinates separately with e.g. `position = ggplot2::position_dodge(width = 10)`. Incompatible local position objects now raise an error. |
| Track statistics and provenance | For example, use `geom_track(aes(x = pos), data = events, geom = "bar", stat = ggplot2::StatBin, binwidth = 10, track = position_track("events"))`. Exported `source_rows` are determined before moving the geometry. Each summary or weighted count lists its eligible observations. |
| Interval coverage | Use `geom_track_bar(data = intervals, stat = "coverage", binwidth = 100, track = position_track("coverage"))` to draw the fraction covered by the union in each window. Map `aes(y = after_stat(mean_depth))` to count overlapping records separately in the depth. |
| Lane direction | `base_position` places the first lane and `direction` sets the growth of successive lanes. Legacy `side = "outside"` corresponds to `base_position = position_plasmid(-0.1), direction = "inward"`. Legacy `side = "inside"` corresponds to `base_position = position_plasmid(0.1), direction = "outward"`. The default follows strand and existing calls retain their placement. |
| Staged colours | On a native track, `aes(fill = stage(start = condition, after_scale = scales::alpha(fill, 0.5)))` groups and trains the track fill scale using condition, then adjusts the resulting colours. Source column names may contain after_stat or after_scale. |

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

The constructor arguments and example in this section apply to v0.17.0.
For the v0.50.0 equivalent, see the [runnable feature and gene examples](#replacement).

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

## Scope of the v0.16.5–v0.16.6 migration {#scope}
The 0.16.0–0.16.4 phases were delivered together in v0.16.5. The v0.16.6
release candidate audit retained flexible parameter formats, spelling aliases,
reexports and the getter/export entry points. The table records the replacements
for those releases. Where historical evidence establishes only an upper bound,
the removal column gives that bound. For this checkout, see the opening v0.50.0
migration table and the examples below.

## Historical replacement table {#table}
| Older call | Replacement in v0.16.5–v0.16.6 | Removal / guard | Meaning |
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

## Run the equivalent examples in v0.50.0 {#replacement}
Run the following in a fresh R session with ggchord v0.50.0 and dplyr.

```r
library(ggchord)
library(dplyr)

seqs <- seq_data_example |>
  select(accver, length)
# The example data stores the original annotation strands in source_strand.
genes <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand)

plasmid_seq <- plasmid_example_pUC19 |>
  select(accver, length)
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source |>
  select(accver, start, end, strand, type, anno)

p_migrated <- ggchord(plasmid_seq) +
  geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) +
  coord_circular()

p_offsets <- ggchord(seqs) +
  geom_seq() +
  geom_feature(
    data = genes,
    preset = "gene",
    position = position_strand(c("+" = .1, "-" = -.2))
  )

p_migrated
p_offsets
```

Old minus-strand offset magnitudes need a negative sign in explicit Position
values. Map categories with `aes(feature_shape = type)` and a shape scale, or
set `shape` for fixed geometry. Automatic label arrangement uses `label_layout`.
For manual geometry use the corresponding fixed-label layer. Replace
`label_layout = "aligned"` with `"auto"`.

## Interfaces retained by the v0.16.6 audit {#retained}
The v0.16.6 audit retained the five colour/color spelling aliases, upstream
reexports, gene/feature/plasmid shortcuts and both geometry getter/export interfaces.
In v0.50.0, use `geom_feature(preset = "gene")` and
`geom_feature(preset = "plasmid")` for the former shortcuts. See
`help("ggchord_api_contract", package = "ggchord")` and the corresponding function
help for how each current interface handles `NULL`, `NA` and omitted arguments.

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

</details>
