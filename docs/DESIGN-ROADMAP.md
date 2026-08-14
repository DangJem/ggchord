# ggchord version roadmap and technical design

This document records the version roadmap for ggchord and gives
concrete, executable API + implementation designs for the features that
are **not yet implemented**. Implemented features are marked
**\[done\]** with pointers to the code; the rest are designs to be
implemented in the next versions.

Priority order for the whole roadmap: **reliability → usability →
expressiveness → ecosystem**.

------------------------------------------------------------------------

## Version plan

| Version | Theme | Status |
|----|----|----|
| v0.7.0 | Reliability: validation, cleaning, tests, visual regression | \[done\] |
| v0.8.0 | Usability: import helpers, ribbon processing, sequence grouping | A/B/C \[done\], D \[designed\] |
| v0.9.0 | Expressiveness: ribbon direction mapping, highlights, features, crowded layouts | designed |
| v1.0.0 | Ecosystem: Plotly, export, themes, documentation, stable API | designed |

------------------------------------------------------------------------

## v0.7.0 — Reliability & maintenance \[done\]

- [`validate_ggchord_data()`](https://dangjem.github.io/ggchord/reference/validate_ggchord_data.md)
  — structured validation (`R/validate.R`).
- [`clean_ggchord_data()`](https://dangjem.github.io/ggchord/reference/clean_ggchord_data.md)
  — report-driven cleaning (`R/clean.R`).
- `ggchord(validate = "warn" | "error" | "none")` — single summary
  warning, cached report (`R/ggchord.R`).
- Test matrix and lightweight visual regression
  (`tests/testthat/test-*.R`, `helper-visual.R`).

## v0.8.0 — Usability & data convenience

### A. Import helpers \[done\] — `R/data-import.R`

[`read_blast()`](https://dangjem.github.io/ggchord/reference/read_blast.md),
[`read_gff3()`](https://dangjem.github.io/ggchord/reference/read_gff3.md),
[`read_fasta_lengths()`](https://dangjem.github.io/ggchord/reference/read_fasta_lengths.md).
No tidyverse dependency; optional columns preserved; clear diagnostics
for malformed files.

### B/C. Ribbon filtering, dedup, merge \[done\] — `R/ribbon-utils.R`

[`filter_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/filter_ggchord_ribbons.md),
[`deduplicate_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/deduplicate_ggchord_ribbons.md),
[`merge_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/merge_ggchord_ribbons.md).
All return `list(data, report)` and attach `source_rows` to `data`.

### D. Sequence grouping (`seq_group`) \[designed\]

Goal: group sequences (host/phage, chromosome sets, sample groups) with
visual separation: a larger inter-group gap, group labels and optional
group colors, without changing the layered API.

**Public API (draft):**

``` r

# Option 1: group column in seq_data (preferred)
seq_data$seq_group <- c("phage", "phage", "host", "host")

ggchord(seq_data, ribbon_data, gene_data) +
  geom_seq(
    seq_group = NULL,               # column name in seq_data or named vector
    seq_group_gap = 0.08,           # extra gap proportion BETWEEN groups
    seq_group_labels = TRUE,        # TRUE | FALSE | named character vector
    seq_group_label_radius = 1.35,
    seq_group_colors = NULL         # named vector by group
  )

# Option 2 (helper, returns a re-ordered seq_data with metadata):
seq_data2 <- ggchord_group_seq(seq_data, group = "seq_group",
                               order = c("host", "phage"))
```

**Layout changes
([`compute_chord_layout()`](https://dangjem.github.io/ggchord/reference/compute_chord_layout.md)):**

1.  After `seq_order` is resolved, split the ordered sequences into
    groups by the group vector (length 1 or n, or named by seq_id —
    reuse
    [`process_sequence_param()`](https://dangjem.github.io/ggchord/reference/process_sequence_param.md)).
2.  Angle allocation: instead of one `sum(seq_gap)`, add an extra gap of
    `seq_group_gap` between consecutive groups:
    `total_gap_prop = sum(seq_gap) + seq_group_gap * (n_groups - 1)`.
3.  New layout element `group_arcs`/`group_labels`: for each group, the
    label anchor is the midpoint angle between the first sequence start
    and the last sequence end, at radius
    `seqRadius * seq_group_label_radius`, rotated like `seq_labels_df`.
4.  Group coloring (optional): when `seq_group_colors` is set, add a
    thin “group band” polygon or color the group labels; the Seq ID
    legend is unchanged (no new aesthetic collisions — use a dedicated
    internal aesthetic like the ribbon `zfill`, or plain text colour).
5.  [`classify_ggchord_layers()`](https://dangjem.github.io/ggchord/reference/classify_ggchord_layers.md)
    /
    [`extract_ggchord_layer_data()`](https://dangjem.github.io/ggchord/reference/extract_ggchord_layer_data.md)
    gain a `group_label` case so `ggplotly()` and lazy data see the same
    geometry.

**Faces vs chord semantics:** faceting (`facet_wrap`/`facet_grid`)
breaks the single-center coordinate space, so it is **not** recommended
for chord diagrams. The supported alternatives are: (1) inter-group
gaps + group labels (above); (2) multiple chord panels arranged by the
user via `gridExtra` / `patchwork`; (3) a multi-ring layout (see
v0.9.0-D).

**Tests:** grouping changes gap allocation only (angles of the last
sequence of group k shift by the group gaps); fingerprints must stay
reproducible; group labels render; unknown group names error.

------------------------------------------------------------------------

## v0.9.0 — Expressiveness & complex layouts

### A. Ribbon direction & visual mapping \[designed\]

Extend
[`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md)
while keeping `"pident"`, `"query"`, `"subject"`, `"single"` behaviour
byte-for-byte compatible:

``` r

geom_ribbon(
  ribbon_color_scheme = "pident",      # unchanged
  ribbon_color_by = NULL,              # NEW: any numeric column, e.g. "bitscore"
  ribbon_color_limits = NULL,          # scale limits for continuous mapping
  ribbon_alpha_by = NULL,              # NEW: numeric column -> alpha
  ribbon_outline_by = NULL,            # NEW: discrete column -> outline colour
  ribbon_orientation = c("auto", "same", "reverse", "both"),
  ribbon_mix_colors = FALSE            # NEW: blend query+subject colors
)
```

**Design notes:**

- In
  [`compute_chord_layout()`](https://dangjem.github.io/ggchord/reference/compute_chord_layout.md),
  when `ribbon_color_by` is set, the polygon gets a continuous value
  column instead of `pident`;
  [`make_ggchord_scales()`](https://dangjem.github.io/ggchord/reference/make_ggchord_scales.md)
  adds a `scale_fill_stepsn()` with limits from `ribbon_color_limits`
  (default: data range) and the title from the column name. The
  Identity(%) legend is then not added — the ribbon legend title changes
  instead.
- Direction detection: a row is “reverse” when
  `(qstart > qend) != (sstart > send)` (or use the `sstrand` column when
  present). `"same"`/`"reverse"` draw only matching rows; `"both"` adds
  an outline/alpha split: two polygon layers per ribbon (same geometry,
  different outline/linetype/alpha), or one polygon with a `direction`
  column mapped to `linetype`/`alpha` — prefer the single-layer mapping
  to keep the legend count stable.
- Query/subject mixing (`ribbon_mix_colors = TRUE`): per-polygon fill =
  gradient between `seq_colors[q]` and `seq_colors[s]` (`colorRamp` at
  t=0.5), or a vertical two-tone gradient via
  [`grid::linearGradient()`](https://rdrr.io/r/grid/patterns.html)
  (needs a custom geom draw method; keep the identity-scale path for
  plotly).
- **Aesthetic isolation:** ribbon visual mappings must stay on the
  internal `zfill`/`zalpha`/`zlinetype` aesthetics (see
  `rename_fill_geom()`), so the gene legend is never polluted.
  [`ggchord_plotly_ggplot()`](https://dangjem.github.io/ggchord/reference/ggchord_plotly_ggplot.md)
  pre-maps the new columns in
  [`ggchord_plotly_ggplot()`](https://dangjem.github.io/ggchord/reference/ggchord_plotly_ggplot.md).
- Tests: each scheme produces distinct layout fingerprints; legend
  titles correct; gene legend unchanged; plotly conversion does not
  error.

### B. Highlight layers \[designed\]

``` r

geom_seq_region(
  data = NULL, mapping = NULL,
  regions = NULL,                    # data.frame(seq_id, start, end[, label, category, color])
  region_color = "#F59E0B",
  region_fill = "#F59E0B",
  region_alpha = 0.25,
  region_width = 0.08,
  region_offset = 0,
  region_side = c("inside", "outside", "auto"),
  show_legend = FALSE,
  legend_position = NULL,
  ...
)

geom_ribbon_highlight(
  data = NULL, mapping = NULL,
  ribbon_ids = NULL,                 # original row numbers or ribbon group IDs
  predicate = NULL,                  # a SAFE function(row) -> logical (no string eval)
  highlight_color = "#E11D48",
  highlight_alpha = 0.8,
  highlight_outline_color = NULL,
  highlight_outline_width = 0.3,
  highlight_expand = 0,              # radial expansion of the highlighted band
  show_legend = FALSE,
  ...
)
```

**Design notes:**

- `geom_seq_region` is a normal ggchord layer with
  `ggchord_type = "seq_region"`. In `compute_chord_geometry()`, collect
  `regions` (or `layer$data`) and map each region onto its sequence arc
  with `map_to_curve()` between the angles of `start`/`end`, drawing a
  band polygon at radius `seqRadius + region_offset ± region_width/2`
  (side-aware: `"inside"` → radius \> seqRadius, `"outside"` → radius \<
  seqRadius, `"auto"` → away from the chord center like gene labels).
  The band reuses the sequence’s reference path, so it respects
  `seq_order`, `seq_orientation`, `seq_radius`, `seq_curvature` for
  free.
- `geom_ribbon_highlight` does **not** recompute the layout: it reads
  the cached `layout$ribbon_polys`, selects the polygons whose `group`
  (or an added `source_row` column) matches `ribbon_ids` /
  `predicate(...)`, and draws them again with the highlight style
  *below/above* the ribbon layer (a second
  `ggchord_type = "ribbon_highlight"` layer whose data is derived, not
  re-layouted). `predicate` is a function taking a one-row data.frame
  and returning TRUE/FALSE; no `eval(parse(...))`.
- Both layers must register empty placeholder data, be listed in
  [`classify_ggchord_layers()`](https://dangjem.github.io/ggchord/reference/classify_ggchord_layers.md)
  and
  [`extract_ggchord_layer_data()`](https://dangjem.github.io/ggchord/reference/extract_ggchord_layer_data.md),
  and be added to
  [`ggchord_plotly_ggplot()`](https://dangjem.github.io/ggchord/reference/ggchord_plotly_ggplot.md)
  (pre-mapped colors) and the visual-regression set. Empty matches
  return the empty placeholder → safe no-op.

### C. Generic feature layer \[designed — keep geom_gene()\]

[`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
stays; `geom_feature()` is a generalization:

``` r

geom_feature(
  data = NULL, mapping = NULL,
  feature_type = NULL,               # column name, default "type"
  feature_shape = c("arrow", "block", "chevron", "lollipop"),
  feature_shape_by = NULL,           # column name -> shape
  feature_color_scheme = c("type", "category", "manual"),
  feature_width = 0.05,
  feature_offset = 0.1,
  show_legend = TRUE,
  ...
)
```

Input columns: `seq_id, start, end, strand, type, label, category`
(label → anno). Internally it reuses the gene-arrow path builder in
`layout.R` (parameterize the head/body polygon shape).
[`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
becomes a thin wrapper (`feature_type = "gene"`, arrow shape) so nothing
breaks. Tests:
[`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
output is byte-identical before/after the wrapper refactor.

### D. Crowded layouts \[designed\]

- **Auto ring placement:** `geom_seq(seq_ring = "auto")` groups
  sequences into nested rings when `n` is large; each ring gets its own
  `seqRadius` range. All automatic choices must be reproducible (seed),
  explainable (the chosen assignment is returned in
  `get_chord_layout()$rings`), overridable (`seq_ring = NULL` disables;
  explicit `seq_radius`/`seq_order` win).
- **Label avoidance:** cap `ggchord_repel_labels()` iterations, add
  `max_overlaps` sampling and a “hide labels beyond N” strategy for huge
  datasets (already partially present via `max_overlaps`).
- **Ribbon bundling/aggregation:** add
  `geom_ribbon(ribbon_reduce = c("none", "sample", "bundle", "density"))`.
  “sample” draws a stratified subsample; “bundle” merges nearby ribbons
  into a single band (reuse
  [`merge_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/merge_ggchord_ribbons.md)
  on the layout geometry); “density” tints the chord interior by
  alignment density. Quality-first (default) vs performance-first
  (`ribbon_reduce = "sample"`) policies documented.
- Performance baseline test: `ribbon_data` with 5k rows must build in a
  bounded time; layout caching (already on `p$ggchord$layout`) must be
  reused, and caches must never leak between plots (no global mutable
  state beyond the
  [`get_chord_layout()`](https://dangjem.github.io/ggchord/reference/get_chord_layout.md)
  inspection cache).

------------------------------------------------------------------------

## v1.0.0 — Ecosystem, interaction, export, documentation

### A. Plotly enrichment \[designed\]

Keep
[`ggplotly.ggchord()`](https://dangjem.github.io/ggchord/reference/ggplotly.ggchord.md)’s
strategy (convert geometry → pre-mapped standard ggplot →
[`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)),
then enrich the traces:

1.  **Hover text** (computed in
    [`ggchord_plotly_ggplot()`](https://dangjem.github.io/ggchord/reference/ggchord_plotly_ggplot.md)
    by attaching a `text` column to each data frame before conversion):
    - sequence arcs: `seq_id`, `seq_label`, `length`, `orientation`,
      display order;
    - ribbons: `qaccver`, `saccver`, `qstart`, `qend`, `sstart`, `send`,
      alignment `length`, `pident`, `evalue`/`bitscore` when present,
      ribbon group id + original row number, and the number of merged
      source rows (from
      [`merge_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/merge_ggchord_ribbons.md)’s
      `source_rows` attribute / report);
    - genes/features: `seq_id`, `start`, `end`, `strand`, `anno/label`,
      `type/category`, original row number.
2.  **Toggle:**
    `ggplotly.ggchord(p, hover = TRUE/FALSE, hover_template = ...)`.
    Avoid dumping repeated polygon coordinates (only show the first-row
    metadata via `text` + `hoverinfo = "text"`).
3.  **Click-to-highlight:** intercept plotly `plotly_click` events in JS
    is complex; the stable first step is a documented pattern: use
    `export_ggchord_layout()` (below) to find the clicked ribbon’s
    query/subject intervals, then re-plot with
    `geom_ribbon_highlight(ribbon_ids = ...)`. Ship a
    [`plotly::event_register`](https://rdrr.io/pkg/plotly/man/event_register.html)
    based example vignette; keep it out of the core test suite
    (`skip_if_not_installed("plotly")`).

### B. Layout / data export \[designed\]

``` r

export_ggchord_layout(
  plot,
  include = c("seq", "ribbon", "gene", "labels", "axis"),
  original_data = FALSE
)
# alias: ggchord_export_data()
```

Returns a tidy list (or a single `data.frame` with a `type` column) of
the computed geometry: sequence arcs, ribbon polygons (+ `source_row` /
merged counts), gene polygons (+ `source_row`), gene labels, axis ticks,
seq labels. It builds the plot if needed
([`prepare_ggchord_plot()`](https://dangjem.github.io/ggchord/reference/prepare_ggchord_plot.md)),
returns *copies* (never modifies `plot`), and documents the coordinate
system: Cartesian data units, origin = chord center, global `rotation`
applied, unit = 1 data unit per unit radius, `coord_fixed(ratio = 1)`.
This is the stable contract that power-user custom layers, SVG/JSON
exporters and interactive tools build on.

### C. Publication-ready output \[designed\]

``` r

theme_ggchord()             # current default, factored out (no visual change)
theme_ggchord_minimal()
theme_ggchord_dark()
theme_ggchord_paper()
theme_ggchord_nature()      # Nature-style, colorblind-friendly
scale_color_ggchord_seq()   # OKabe-Ito / Set1 alternatives
scale_fill_ggchord_identity()
```

Rules: default appearance unchanged unless the user calls a theme;
themes only touch theme elements (`theme()`), never scales; `scale_*`
helpers are plain `scale_*` wrappers so user-supplied `scale_*()` added
later still wins
([`attach_ggchord_scales()`](https://dangjem.github.io/ggchord/reference/attach_ggchord_scales.md)
already skips aesthetics with a user scale). Add SVG/PDF text-quality
notes (use `cairo_pdf`/`svglite`) in the vignette.

### D. Documentation & stable API

README/vignette additions: “Data validation and cleaning”, “Importing
FASTA / BLAST / GFF3 data”, “Ribbon filtering, deduplication and
merging”, “Highlights and features”, “Interactive Plotly output”,
“Exporting final layout data”, “Publication-ready themes and SVG/PDF
output”, and an FAQ covering dense plots, label overlap, wrong
orientation, too many ribbons, legend conflicts, Plotly-vs-static
differences, large-data performance and how to use exported layouts in
custom layers. Every documented parameter gets a minimal runnable
example; the full vignette walks from FASTA/BLAST/GFF3 to a publication
figure.

------------------------------------------------------------------------

## Acceptance & performance rules (all versions)

- `testthat::test_local(".", reporter = "summary")` must pass; every new
  feature needs success, error, boundary and compatibility tests.
- `R CMD check` must pass; optional `Suggests` (plotly etc.) skipped
  cleanly via `skip_if_not_installed()`.
- At least thousands of ribbons/genes should be considered; layout
  computation and filtering must be cached/reused and never recomputed
  needlessly; no cache may leak between plots.
- README examples must run with output written to temp files (never
  overwriting repo images); external command-line steps are skipped in
  tests.
- Public API, README examples and existing parameter behaviour stay
  backward compatible; new layers must work with lazy data,
  `ggplot_build()` and
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).
