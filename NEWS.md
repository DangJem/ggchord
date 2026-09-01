# ggchord 0.10.0

## Highlights

* Dense alignments can now be simplified explicitly with
  `bundle_ggchord_ribbons()`, while `optimize_ggchord_layout()` searches for a
  deterministic sequence order and orientation with fewer weighted crossings.
  `geom_ribbon()` remains conservative: one input row still draws one ribbon
  unless the user requests preprocessing.

* `geom_feature()` supports independent `"arrow"`, `"block"`,
  `"chevron"`, and `"lollipop"` geometries with
  `scale_feature_shape_manual()`. Shapes and legend keys follow each
  sequence's local tangent/normal frame; lollipop heads remain circular on
  curved sequences.

* New `export_ggchord_layout()` provides stable layer geometry with
  `layer_id`, `source_row`, mapped input columns, and coordinate-space
  metadata. New `view_ggchord()` previews the real dimensions produced by
  `ggsave()` without introducing a second saving API or an interactive
  dependency.

## Ribbons, labels, and layout

* The default `ribbon_gap = NULL` resolves clearance per endpoint. Ribbons
  move closer to `geom_seq()` where no gene or feature polygon is present and
  retain clearance where a physical obstacle overlaps the genomic interval.
  Text and leader lines are not treated as ribbon obstacles.

* Ribbon boundaries use a darker neutral outline, and unrequested ribbon or
  gene polygons are no longer computed. Dense-data reports preserve source-row
  mappings, alignment direction, weights, and approximation diagnostics.

* `geom_gene_label()` is the single fixed/manual label layer and now supports
  horizontal, radial, or tangent text; inside, outside, or automatic side
  selection; and hide, nudge, or allow overlap policies. Automatic leader-line
  placement remains in `geom_gene_label_repel()`.

* Coordinate fitting uses independent tight x/y ranges. Text measurement,
  collision boxes, leader clipping, and limits share a device-aware physical
  scale, preventing small exports from being measured as a six-inch canvas.
  The transparent default legend background no longer masks labels extending
  into compact plot margins.

## Visual system and ggplot2 integration

* Default sequence, ribbon, gene, label, axis, and background styles were
  recalibrated for a restrained publication-oriented hierarchy. Strand legend
  keys use the same tapered silhouettes and directions as `geom_gene()`.

* Identity, sequence, strand, and feature guides scale keys, colourbars,
  typography, and margins relative to the output device, with bounded scaling
  for extreme sizes. Explicit guide dimensions remain authoritative.

* Namespace imports are deliberately narrow. Internal ggplot2, grid, and
  grDevices calls are qualified, while ggchord selectively re-exports only
  helpers directly useful for chord plots: `aes()`, `after_stat()`,
  `after_scale()`, `annotate()`, `labs()`, `ggtitle()`, `guides()`,
  `theme()`, the four theme element constructors, `margin()`, `rel()`,
  `ggsave()`, `last_plot()`, `waiver()`, and `expansion()`. Cartesian
  geoms, facets, coordinates, generic scales, and complete theme presets remain
  in ggplot2 so IDE completion stays focused.

* `theme_ggchord()` and its variants no longer repeat Cartesian-axis,
  panel-grid, background, or margin settings already inherited unchanged.

## Breaking changes

* Sequence grouping has been removed: `geom_seq_group_label()`,
  `scale_group_colour_manual()` / `scale_group_color_manual()`, and every
  `seq_group*` argument or aesthetic now produce a clear migration error.

* No separate `geom_gene_label_manual()` is introduced. Use
  `geom_gene_label()` for fixed/manual placement and
  `geom_gene_label_repel()` for automatic layouts.

* Plotly conversion remains removed. v0.10.0 focuses on deterministic static
  ggplot2 output; interactive rendering will be reconsidered only after the
  static API is stable.

# ggchord 0.9.0

## Release audit

* The v0.9.0 public API was audited across constructors, geoms, role-specific
  scales, themes, guides, coordinates, data utilities and import helpers.
  Layer-local data/mappings, same-type layer isolation, old/new scale
  conflicts and plot-owned layout retrieval were rechecked against the final
  v0.9.0 interface.

* The obsolete, unused `ggchord_label_pad()` internal helper and its generated
  help page were removed. Adaptive limits remain the single implementation
  used to fit rendered label boxes.

* `ggchord()` now reuses the structured validator for its always-on safety
  checks instead of maintaining a second validation rule set. Validation also
  stops advanced coordinate and duplicate calculations when malformed numeric
  columns make those calculations unsafe, returning a complete report rather
  than a secondary type error.

## Static rendering focus

* The experimental Plotly conversion method and dependency have been removed.
  v0.9.0 focuses on deterministic ggplot2 output; a future interactive design
  will be considered separately after the static API is stable.

* Default discrete colours now use a colour-vision-friendly palette (with a
  qualitative HCL fallback for larger sets). Strand colours, sequence and gene
  outlines, ribbon separation, highlight colour and legend key glyphs were
  recalibrated for clearer screen, PDF and greyscale output. These defaults
  remain fully replaceable through the role-specific scales and geom styles.

## Data correctness and import fixes

* `clean_ggchord_data(unknown_id = "keep")` now retains unknown gene rows
  without attempting coordinate checks against a missing sequence length.
  Sorting reversed ribbon intervals records their original same/reverse
  direction so drawing does not silently change alignment orientation.

* Ribbon filtering reports every removal reason for a row;
  `deduplicate_ggchord_ribbons(keep = "first")` now means the first input row;
  and `merge_ggchord_ribbons()` no longer leaves stale values in disagreeing
  auxiliary columns. Use `extra_columns = "first"` to request the previous
  first-row behaviour explicitly.

* BLAST outfmt 7 imports now parse and validate `# Fields:` instead of assuming
  a fixed 17-column layout. GFF3 parsing stops at `##FASTA`. All three import
  helpers can add `.source_file` with `source_file = TRUE`.

* Unknown ribbon and gene sequence IDs now have the same severe validation
  level. Skipped duplicate checks for exceptionally large pair groups are
  reported rather than omitted silently. Feature categories, region outlines,
  curved-region side selection and highlight argument validation were fixed.

## Layer-specific data and geometry

* Every ggchord layer now receives a stable `layer_id` and its own geometry
  registry entry. Multiple gene, feature, region, ribbon, highlight, axis or
  label layers no longer reuse the last layer's data and parameters.

* Layer `data` and role mappings such as `aes(seq_id = chromosome, start =
  from)` are evaluated against that layer's input. Original columns are joined
  back to expanded geometry through `source_row`, so ordinary visual mappings
  remain available during the ggplot2 build.

* `get_chord_layout(plot, build = TRUE)` retrieves the layout owned by a
  specific plot. Calling `get_chord_layout()` without a plot still works for
  compatibility but is deprecated because “most recently built plot” is
  ambiguous when plots are built in an interleaved order.

* Sequence reference paths are cached within one build and reused by
  independent same-type layers. The cache is local to that build and cannot
  leak geometry between plots.

## Role-specific scales

* Sequence, group, ribbon, gene, feature and region layers now use independent
  role aesthetics: `seq_colour`, `group_colour`, `ribbon_fill`,
  `ribbon_alpha`, `ribbon_colour`, `ribbon_linetype`, `gene_fill`,
  `feature_fill` and `region_fill`. Their public `scale_*()` constructors can
  coexist in one plot without replacing another layer's fill or colour scale.

* `scale_seq_position_continuous()` controls genomic major/minor breaks and
  labels. It is trained independently against each sequence length.

* `scale_group_color_manual()` and `scale_ribbon_color_manual()` are available
  as American-English aliases of their `colour` counterparts, matching
  ggplot2's spelling convention.

* Old scale-like geom arguments remain functional during v0.9.0 and emit one
  migration warning per session. Supplying both an old argument and the new
  role scale is an error rather than silently choosing one. Principal
  migrations are:

| Old geom argument | New interface |
| --- | --- |
| `seq_colors`, `seq_group_colors` | `scale_seq_colour_manual()`, `scale_group_colour_manual()` |
| `ribbon_colors`, colour limits/breaks/name | `scale_ribbon_fill_*()` |
| `ribbon_*_by` | the corresponding `aes(ribbon_* = ...)` |
| `ribbon_alpha_range` | `scale_ribbon_alpha_continuous(range = ...)` |
| ribbon outline/linetype/direction visual values | ribbon colour/linetype/alpha scales |
| `gene_colors`, `gene_order` | `scale_gene_fill_manual()` |
| `feature_colors`, `feature_order` | `scale_feature_fill_manual()` |
| axis major/minor counts and labels | `scale_seq_position_continuous()` |

## Coordinate, theme and guide helpers

* `coord_chord()` now owns global rotation, aspect ratio, clipping and view
  fitting. `fit = "labels"` measures annotation boxes, `fit = "geometry"`
  fits only geometric marks, and `fit = "manual"` requires explicit limits.
  User-supplied limits take priority, and replacing it with another ggplot2
  coordinate system is no longer silently undone during the build.

* `theme_ggchord()`, `theme_ggchord_minimal()`, `theme_ggchord_dark()` and
  `theme_ggchord_publication()` provide a small set of composable themes.
  Dedicated theme elements are registered for axes, sequence/group labels,
  gene labels and leader segments; data-driven colours remain scales.

* `guide_ggchord_legend()` and `guide_ggchord_colourbar()` are thin wrappers
  around ggplot2 guides with compact chord-diagram defaults. Existing
  `ggchord()` arguments `title`, `rotation`, `panel_margin` and `show_legend`
  remain functional in v0.9.0 but point users to `labs()`, `coord_chord()` and
  `theme()` respectively.

## Geom and annotation interfaces

* Gene, sequence and axis text sizes are now fixed layer values rather than a
  shared `size` scale, so adding one label layer cannot rescale another.
  Registered axis, sequence-label, gene-label and leader-line theme elements
  are resolved before the standard ggplot2 build; an explicit geom style still
  takes priority.

* `geom_axis()` routes shared styles only to compatible child geoms and accepts
  separate `line_params`, `tick_params` and `text_params`. Its former
  `show_legend` argument is removed because axis annotations never participate
  in a legend.

* `geom_seq_group_label()` provides an independent group-label layer, while
  the labels created implicitly by `geom_seq()` remain available for
  compatibility. Group values now train `scale_group_colour_manual()` instead
  of treating already-resolved colour strings as categories.

* `geom_seq_label(labels = ...)` separates displayed sequence text from scale
  labels; `seq_labels` remains a deprecated alias. Label arguments passed to
  `geom_gene()` now fail clearly instead of being warned about and then leaked
  into `geom_polygon()`. Geom-level legend positions and colourbar dimensions
  remain functional during v0.9.0 but direct users to `guides()`.

## Deterministic gene-label layouts

* `geom_gene_label_repel()` now uses the deterministic
  `gene_label_layout = "aligned" | "radial" | "arc"` interface. `"aligned"`
  remains the default and arranges horizontal labels on orderly cardinal
  rails. `"radial"` uses the nearest collision-free local offset track while
  keeping text horizontal. `"arc"` rotates readable text along the sequence
  tangent and draws a short leader only when a label has moved away from its
  first track.

* All modes now use each sequence's real curve and local normal, including
  custom `seq_radius`, `seq_curvature`, `seq_gap`, mixed `seq_orientation`,
  rotation and sequence groups. They share fixed-obstacle avoidance,
  cross-sequence collision handling, order-preserving leader routing and
  device-aware clipping. Rotated labels use oriented-rectangle collision and
  clipping, preventing spurious whitespace and oversized breaks in leaders.

* The layout ideas are informed by orderly multi-sequence callout figures and
  by the external/inside feature-label approaches offered by SnapGene and
  Geneious. ggchord uses its own generic mode names and implementation; it
  does not copy third-party assets or visual designs. See the
  [SnapGene feature-label documentation](https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map)
  and [Geneious label options](https://manual.geneious.com/en/latest/Sequences.html).

## Breaking API simplification

`geom_gene_label_repel()` now has the following focused interface:

```r
geom_gene_label_repel(
  mapping = NULL, data = NULL,
  gene_label_layout = "aligned",
  gene_label_size = NULL, gene_label_wrap = NULL,
  gene_label_side = "outside", max_overlaps = Inf,
  gene_label_segment_linetype = "auto",
  show_legend = FALSE, ...
)
```

Removed arguments fail immediately rather than being silently ignored:

| Removed argument(s) | Migration |
| --- | --- |
| `gene_label_rotation`, `gene_label_radial_offset`, `gene_label_circum_offset`, `gene_label_circum_limit` | Use `geom_gene_label()` for manual rotation or offsets. |
| `box_padding`, `point_padding`, `min_segment_length`, `force`, `seed` | Select an automatic `gene_label_layout`; collision and line settings are managed by the mode. |
| `gene_label_orientation`, `gene_label_segment` | Use `gene_label_layout = "aligned"`, `"radial"`, or `"arc"`. |

This was the v0.9.0 interface. v0.10.0 subsequently strengthens the existing
fixed-position `geom_gene_label()` instead of introducing a separate manual
geom; see the development section above.

# ggchord 0.8.0

## New features: improved label placement and de-overlap

* `geom_seq_label()` now places sequence names on the arc by default
  (`seq_label_radius = 1`) and rotates them along the arc while keeping them
  readable; `seq_label_orientation = "horizontal"` draws every label
  horizontally, extending away from the chord centre.

* `geom_gene_label()` now sits right beside the gene arrows by default
  (`gene_label_radial_offset = 0.04`) and gains `gene_label_wrap` for
  wrapping long annotations into narrower, less overlapping labels.

* `geom_gene_label_repel()` now defaults to `gene_label_orientation =
  "horizontal"`, `gene_label_segment = "elbow"` (an L-shaped leader line that
  adapts to each label's position and text width) and `gene_label_side =
  "outside"`, so labels stay readable and out of the ribbon area. A
  deterministic final de-overlap pass measures the exact rendered text boxes
  and treats the sequence, group and axis labels as hard rectangular
  obstacles; `max_overlaps` hides labels that still collide after repulsion
  (ggrepel-style decluttering).

* The label text-box projection is now shared by the repulsion solver, the
  obstacle boxes and the coordinate limits, so all three agree on where text
  will actually be drawn.

## New features: adaptive plot limits

* Plot limits now fit the rendered text boxes instead of adding one global
  text-width pad on every side. The actual gene/sequence/group/axis label
  boxes are measured and only the sides that need it are expanded, reducing
  empty margins and using the panel area more efficiently.

## New features: sequence grouping

* `geom_seq()` gains sequence-grouping support via `seq_group`,
  `seq_group_gap`, `seq_group_labels`, `seq_group_label_radius` and
  `seq_group_colors`. Groups can come from a `seq_group` column in `seq_data`
  or be supplied as a single value, a named/positional vector, or a list.

* An extra inter-group gap (`seq_group_gap`) is inserted only at boundaries
  between different groups, and optional group labels are drawn at the
  angular midpoint of each group, at a customisable radius.

* Group labels are rendered horizontally and use their own internal
  `zcolour` identity scale, so they never interfere with the Seq ID colour
  legend. `geom_seq()` stays backward compatible and still returns a single
  layer; group labels are appended lazily at build time.

## New features: ribbon visual mappings and direction

* `geom_ribbon()` can now map any numeric column to a continuous fill via
  `ribbon_color_by` (for example `"bitscore"` instead of `pident`), with
  `ribbon_color_limits`, `ribbon_color_breaks` and `ribbon_color_name` to
  control the colourbar.

* `ribbon_alpha_by` / `ribbon_alpha_range` scale ribbon transparency
  continuously from a numeric column.

* `ribbon_outline_by` / `ribbon_outline_colors` and `ribbon_linetype_by` /
  `ribbon_linetypes` map discrete columns to outline colour and linetype
  using internal aesthetics, without disturbing the Seq ID or Identity(%)
  legends.

* `ribbon_direction` (one of `"none"`, `"alpha"`, `"outline"` or
  `"linetype"`) visually distinguishes same- vs reverse-orientation
  alignments, with `ribbon_direction_colors`, `ribbon_direction_linetypes`
  and `ribbon_direction_alpha` for fine control.

* `legend_key_width` / `legend_key_height` control the size of the Identity(%)
  colourbar key.

## New features: highlights and generic features

* New `geom_seq_region()` draws rectangular bands along sequence arcs to mark
  loci, repeats, CRISPR arrays or other user-defined intervals. It accepts
  `seq_id`, `start` and `end` (plus optional `label`, `category` and `color`)
  and exposes `region_fill`, `region_color`, `region_alpha`, `region_width`,
  `region_offset` and `region_side`.

* New `geom_ribbon_highlight()` emphasizes selected ribbons without changing
  the underlying Identity(%) legend. Selection uses safe, explicit filters
  (`ribbon_ids`, query/subject IDs, pident/length ranges, or a predicate
  function) and reuses the computed ribbon geometry.

* New `geom_feature()` is a thin, backwards-compatible convenience layer for
  CDS, tRNA, rRNA, repeat, CRISPR, promoter and custom feature tables; it
  prepares a gene-compatible table and reuses `geom_gene()`'s geometry and
  scales, with `feature_colors`, `feature_width`, `feature_offset` and
  `feature_order` for styling.

## Documentation

* Added man pages and runnable examples for the new layers
  (`geom_seq_region()`, `geom_ribbon_highlight()`, `geom_feature()`) and
  expanded the documentation for the updated `geom_seq()`, `geom_ribbon()`,
  `geom_seq_label()` and `geom_gene_label_repel()` parameters.

# ggchord 0.7.0

## New features: structured data validation and cleaning

* New exported function `validate_ggchord_data()` returns a structured
  `ggchord_validation` object: a `valid` flag, `errors` (severe problems),
  `warnings` (drawable but noteworthy issues), per-category `summary` counts,
  a `data_summary` (sequences/ribbons/genes, unknown IDs, out-of-range rows,
  ...), the original row numbers of every problem (`invalid_rows`) and the
  automatically fixable issues (`cleanable`). `print()` and `summary()`
  methods are provided; `strict = TRUE` stops on severe problems.

* New exported function `clean_ggchord_data()` applies explicit, conservative
  policies (`unknown_id`, `out_of_range`, `reversed_interval`,
  `invalid_pident`, `empty_annotation`) and returns the cleaned tables plus a
  full report of every change (original row number, reason, original/new
  values, action). The input objects are never modified and nothing is dropped
  silently.

* `ggchord()` gains a `validate = c("warn", "error", "none")` argument.
  The default `"warn"` emits a single summary warning (never one warning per
  row) and caches the full report on the plot (`p$ggchord$validation`);
  `"error"` stops on severe problems; `"none"` keeps a fast path. Valid input
  renders exactly as before.

## New features: data import and ribbon preparation

* `read_blast()` parses BLAST `-outfmt 6/7` tabular output (12 or 17 columns,
  auto-detected) into `ribbon_data` format, preserving `evalue`, `bitscore`,
  `qcovs`, `qlen`, `slen`, `sstrand` and `stitle` when present.

* `read_gff3()` parses GFF3 files into `gene_data` format, selecting
  `feature_types` (default `CDS`), extracting `anno` from attribute keys
  (`product`, `Name`, ...), decoding percent-encoding, and mapping
  unstranded features to `+` (or dropping them).

* `read_fasta_lengths()` reads FASTA headers and sequence lengths into
  `seq_data` format, with optional `header_delim` splitting.

* `filter_ggchord_ribbons()` filters ribbons by sequence IDs, pident, length,
  E-value, bitscore, query/subject coverage, undirected sequence pairs and
  self-links, with optional sorting; missing columns produce clear errors.

* `deduplicate_ggchord_ribbons()` removes exact, coordinate-near or highly
  overlapping duplicate blocks (`by = "exact" | "coordinates" | "overlap"`)
  keeping the best pident, longest, or first representative.

* `merge_ggchord_ribbons()` merges adjacent/overlapping blocks of the same
  sequence pair with length-weighted pident. Merging is deliberately
  conservative: blocks with inconsistent spans, large pident differences or
  different orientations are left unmerged.

* All ribbon utilities keep extra columns and the original column order,
  attach the original row numbers as the `source_rows` attribute, and return a
  report of what was removed/merged and why.

## Testing

* Added test files covering validation, cleaning, the `validate` integration,
  data import, ribbon filtering/deduplication/merging and a lightweight
  visual-regression suite (deterministic layout fingerprints plus an opt-in
  PNG md5 baseline behind `GGCHORD_VISUAL_REGRESSION=1`).

# ggchord 0.6.1

(No user-facing changes; internal bug fixes.)

# ggchord 0.6.0

## New features
* Plot objects are now fully self-contained: data and parameters are stored on
  the plot itself instead of in a package-wide environment. Multiple plots can
  be created and built independently in the same session, and plots survive
  `saveRDS()` / `readRDS()`.

* The layout is now computed by `ggplot_build()` rather than by a custom
  `print()` method. As a result `print()`, `ggsave()`, `ggplot_build()` and
  other standard ggplot2 workflows all work directly on ggchord plots, and
  rendering no longer modifies the user's plot
  object.

* New layer `geom_seq_label()`: places sequence labels at the midpoint of each
  sequence arc with control over radial offset (`seq_label_radius`), rotation
  (`seq_label_rotation`) and font size (`seq_label_size`).

* New ribbon color scheme `"subject"`: colors ribbons by the subject sequence
  (`saccver`), complementing the existing `"query"` scheme.

* The layout accessor `get_chord_layout()` is now exported, making the computed
  geometry available for custom layers and annotations.

* Themes, scales and other ggplot2 objects can now be added with `+` (e.g.
  `p + theme(legend.position = "bottom")`), and user-supplied colour/fill
  scales are respected instead of being overwritten.

* `ggchord()` now warns about suspicious input data: reversed or out-of-range
  alignment/gene coordinates, `pident` outside [0, 100], and sequence IDs that
  are not present in `seq_data`.

* `geom_gene_label_repel()` gains `gene_label_side = "auto" | "inside" |
  "outside"`. With `"outside"`, labels that would sit inside the chord (where
  they can overlap the ribbons) are mirrored to the outside of their sequence
  arc, keeping the same radial distance from the arc.

* New `gene_label_segment_linetype` argument controls the leader-line linetype.
  The default `"auto"` draws solid lines, except for labels that were moved to
  the other side of their arc (`gene_label_side`), which are drawn dashed. Any
  other valid ggplot2 linetype (e.g. `"dotted"` or a numeric dash pattern) is
  applied to all leader lines.

* Elbow leader lines no longer force fixed segment lengths: the stub scales
  with each label's text width and the horizontal space available between the
  gene and the label, so labels can be placed more flexibly without
  degenerate (near-zero) stubs.

* `geom_seq_label()` now documents and follows the intended `seq_label_radius`
  semantics: `1` sits on the arc, `> 1` places the label outside (away from
  the chord center) and `< 1` inside. Previously the multiplier was applied in
  the opposite direction (the default `1.15` put labels *inside* the chord).

* New `geom_seq_label()` options: `seq_label_orientation = "arc" | "horizontal"`
  (horizontal labels extend away from the chord center),
  `seq_label_hjust` / `seq_label_vjust` for per-sequence justification, and
  `check_overlap` to skip labels that would overlap.

* The default theme no longer draws grid lines (`panel.grid` is blank) and
  legend keys are transparent (they blend into the plot background instead of
  a fixed white rectangle).

## Performance
* Replaced the linear angle lookup in the layout mapping with a binary search
  (`findInterval`), speeding up layout computation for large plots.

## Dependency changes
* Declares `ggplot2 (>= 4.0.0)` and `R (>= 4.1.0)` to match the implementation
  (the package relies on ggplot2 4.x internals).

## Infrastructure
* Added a GitHub Actions `R CMD check` workflow (macOS, Windows, Linux).
* Removed the internal legacy `fill_ggnewscale_1` aesthetic name in favour of
  `fill_ribbon`.

## Bug fixes
* Tests no longer write to a hard-coded `/tmp` path: they use `tempfile()`, so
  the test suite passes on Windows and leaves no stray files behind for
  `R CMD check` (fixes the CRAN incoming-check failure).

* The Identity(%) colourbar no longer collapses into a thin/invisible line when
  the legend is placed at the top/bottom or the legend box is horizontal
  (`legend.box = "horizontal"`). It now follows the theme's legend position: a
  vertical bar filling the available height at the left/right, and a fixed-size
  horizontal bar at the top/bottom.

* Legend keys are transparent and do not inherit `panel.background` (ggplot2
  4.x lets unset legend keys follow the panel background, so the key fill is
  set explicitly to stay transparent).

* Sequence (and gene) labels no longer end up upside down when a global
  `rotation >= 90` is used: the readability flip is now re-applied after the
  layout rotation instead of only before it.

* The repulsion spring now pulls labels toward their own starting position
  rather than the leader-line anchor, which keeps labels moved with
  `gene_label_side = "outside"` on the outside while their leader line still
  starts at the gene.

* With `gene_label_side`, every label is kept on the requested side of its arc
  (previously only the labels moved by the side switch were re-checked, so a
  crowded repulsion layout could push other labels across the arc).

* The built-in `gene_data_example` annotations no longer contain URL-encoded
  `%2C` artifacts (e.g. "ribonucleotide reductase%2C large subunit" is now
  "ribonucleotide reductase large subunit").

## New features
* Each legend can now be positioned independently via the `legend_position`
  argument of `geom_seq()`, `geom_ribbon()` and `geom_gene()` (e.g.
  `geom_ribbon(legend_position = "bottom")`). Legends without an explicit
  position stay together at `theme(legend.position = ...)`.

* Parameter specification is now more flexible and human-friendly. Sequence
  parameters accept a single value, vectors, vectors/lists named by sequence
  ID, lists named by sequence order (`"1"`, `"2"`, ...) and unnamed lists;
  gene parameters additionally accept per-strand (`+`/`-`) specifications in
  any of those forms (e.g. `gene_label_rotation = c("+" = -15, "-" = -45)` or
  `list(c("+" = -15, "-" = -45), ...)`), including length-one lists that
  recycle.

# ggchord 0.5.0


## New features
* Added ribbon outline customization to `geom_ribbon()`:
  `ribbon_outline_color` (default `"black"`), `ribbon_outline_width`
  (default `0.05`) and `ribbon_outline_linetype` (default `1`, solid).

## Dependency changes
* Removed the `ggnewscale` dependency. The ribbon and gene fill scales are now
  kept independent via an internal renamed-fill aesthetic, so no external
  package is required for plots with both ribbon and gene layers.
* Removed the `RColorBrewer` dependency. The default Set1 categorical palette
  is now built into the package, so the rendered default colors are unchanged.

## Bug fixes
* Fixed the ribbon fill scale being overwritten by the gene fill scale when
  both `geom_ribbon()` and `geom_gene()` were present (previously produced
  wrong ribbon colors and a "Scale for fill is already present" message).
* Fixed `ribbon_alpha` rendering at the wrong opacity (e.g. `0.35` was drawn
  as ~0.55); the alpha value is now applied exactly as specified.
* Fixed `geom_axis(show_axis = FALSE)` failing with "object 'label' not found".
* Fixed `axis_label_orientation` rejecting mixed vectors such as
  `c("horizontal", 45, ...)`.
* Fixed warnings from `brewer.pal()` when fewer than 3 sequences or gene
  annotations are used (two-sequence plots now render cleanly).
* Fixed an error when `geom_gene()` was added before `geom_ribbon()`
  ("Continuous value supplied to a discrete scale").
* Fixed plots containing only `geom_axis()` where the axis path was
  misclassified as a sequence arc.
* Registered `+.ggchord` and `ggplot_build.ggchord()` as proper S3 methods and
  aligned the `ggplot_build.ggchord()` signature with the generic.

## Documentation
* Translated all code comments and user-facing messages to English.
* Added man pages for previously undocumented exported functions
  (`geom_seq()`, `geom_ribbon()`, `geom_gene()`, `geom_axis()`,
  `coord_chord()`, `+.ggchord()` and others).
* README: added per-column data preparation tables with example rows,
  rendered example plots under `examples/plots/`, generalized the package
  description beyond BLAST, and documented the ribbon outline parameters.
* Rewrote the package vignette for the layered `ggchord() + geom_*` API.

# ggchord 0.4.0

## Changes
* Parameter redistribution: layout parameters moved from `ggchord()` into the
  individual `geom_*` layers; `ggchord()` now only validates data and sets
  global style (`title`, `rotation`, `panel_margin`, `show_legend`, `debug`).
* Deferred computation: the coordinate layout is computed at `print()` time,
  collecting parameters from all layers during rendering.
* Custom `print.ggchord()` method: merge parameters, compute the layout,
  inject data into layers, then render.
* Added 15 unit tests.

# ggchord 0.3.0

* Layered API refactoring: split the monolithic function into
  `ggchord() + geom_seq() + geom_ribbon() + geom_gene() + geom_axis()`.
* Custom `+.ggchord` method that flattens layer lists automatically.
* Lightweight `coord_chord()` coordinate system.

# ggchord 0.2.0

* Enhanced arc and line mode optimization.
* Precise curvature and gap control.
* Enhanced color customization.

# ggchord 0.1.0

* Separate management of sequence, alignment, and gene data.
* Sequence orientation, custom order, gap and radius adjustment.
* Customizable axes; ribbons support 3 coloring schemes.

# ggchord 0.0.2

* Multi-sequence support; arc/line mode switching.

# ggchord 0.0.1

* Initial release: pairwise alignment chord diagram visualization.
