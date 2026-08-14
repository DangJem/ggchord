# Changelog

## ggchord 0.7.0

### New features: structured data validation and cleaning

- New exported function
  [`validate_ggchord_data()`](https://dangjem.github.io/ggchord/reference/validate_ggchord_data.md)
  returns a structured `ggchord_validation` object: a `valid` flag,
  `errors` (severe problems), `warnings` (drawable but noteworthy
  issues), per-category `summary` counts, a `data_summary`
  (sequences/ribbons/genes, unknown IDs, out-of-range rows, …), the
  original row numbers of every problem (`invalid_rows`) and the
  automatically fixable issues (`cleanable`).
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods are
  provided; `strict = TRUE` stops on severe problems.

- New exported function
  [`clean_ggchord_data()`](https://dangjem.github.io/ggchord/reference/clean_ggchord_data.md)
  applies explicit, conservative policies (`unknown_id`, `out_of_range`,
  `reversed_interval`, `invalid_pident`, `empty_annotation`) and returns
  the cleaned tables plus a full report of every change (original row
  number, reason, original/new values, action). The input objects are
  never modified and nothing is dropped silently.

- [`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
  gains a `validate = c("warn", "error", "none")` argument. The default
  `"warn"` emits a single summary warning (never one warning per row)
  and caches the full report on the plot (`p$ggchord$validation`);
  `"error"` stops on severe problems; `"none"` keeps a fast path. Valid
  input renders exactly as before.

### New features: data import and ribbon preparation

- [`read_blast()`](https://dangjem.github.io/ggchord/reference/read_blast.md)
  parses BLAST `-outfmt 6/7` tabular output (12 or 17 columns,
  auto-detected) into `ribbon_data` format, preserving `evalue`,
  `bitscore`, `qcovs`, `qlen`, `slen`, `sstrand` and `stitle` when
  present.

- [`read_gff3()`](https://dangjem.github.io/ggchord/reference/read_gff3.md)
  parses GFF3 files into `gene_data` format, selecting `feature_types`
  (default `CDS`), extracting `anno` from attribute keys (`product`,
  `Name`, …), decoding percent-encoding, and mapping unstranded features
  to `+` (or dropping them).

- [`read_fasta_lengths()`](https://dangjem.github.io/ggchord/reference/read_fasta_lengths.md)
  reads FASTA headers and sequence lengths into `seq_data` format, with
  optional `header_delim` splitting.

- [`filter_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/filter_ggchord_ribbons.md)
  filters ribbons by sequence IDs, pident, length, E-value, bitscore,
  query/subject coverage, undirected sequence pairs and self-links, with
  optional sorting; missing columns produce clear errors.

- [`deduplicate_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/deduplicate_ggchord_ribbons.md)
  removes exact, coordinate-near or highly overlapping duplicate blocks
  (`by = "exact" | "coordinates" | "overlap"`) keeping the best pident,
  longest, or first representative.

- [`merge_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/merge_ggchord_ribbons.md)
  merges adjacent/overlapping blocks of the same sequence pair with
  length-weighted pident. Merging is deliberately conservative: blocks
  with inconsistent spans, large pident differences or different
  orientations are left unmerged.

- All ribbon utilities keep extra columns and the original column order,
  attach the original row numbers as the `source_rows` attribute, and
  return a report of what was removed/merged and why.

### Testing

- Added test files covering validation, cleaning, the `validate`
  integration, data import, ribbon filtering/deduplication/merging and a
  lightweight visual-regression suite (deterministic layout fingerprints
  plus an opt-in PNG md5 baseline behind `GGCHORD_VISUAL_REGRESSION=1`).

## ggchord 0.6.1

(No user-facing changes; internal bug fixes.)

## ggchord 0.6.0

### New features

- Plot objects are now fully self-contained: data and parameters are
  stored on the plot itself instead of in a package-wide environment.
  Multiple plots can be created and built independently in the same
  session, and plots survive
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html).

- The layout is now computed by
  [`ggplot_build()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)
  rather than by a custom [`print()`](https://rdrr.io/r/base/print.html)
  method. As a result [`print()`](https://rdrr.io/r/base/print.html),
  [`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html),
  [`ggplot_build()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)
  and other standard ggplot2 workflows
  (e.g. [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html))
  all work directly on ggchord plots, and rendering no longer modifies
  the user’s plot object.

- New layer
  [`geom_seq_label()`](https://dangjem.github.io/ggchord/reference/geom_seq_label.md):
  places sequence labels at the midpoint of each sequence arc with
  control over radial offset (`seq_label_radius`), rotation
  (`seq_label_rotation`) and font size (`seq_label_size`).

- New ribbon color scheme `"subject"`: colors ribbons by the subject
  sequence (`saccver`), complementing the existing `"query"` scheme.

- The layout accessor
  [`get_chord_layout()`](https://dangjem.github.io/ggchord/reference/get_chord_layout.md)
  is now exported, making the computed geometry available for custom
  layers and annotations.

- Themes, scales and other ggplot2 objects can now be added with `+`
  (e.g. `p + theme(legend.position = "bottom")`), and user-supplied
  colour/fill scales are respected instead of being overwritten.

- [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)
  now works on any ggchord plot, including plots that combine the ribbon
  and gene layers (previously this raised a scale error). A dedicated
  [`ggplotly.ggchord()`](https://dangjem.github.io/ggchord/reference/ggplotly.ggchord.md)
  method converts the computed geometry to a plotly-friendly plot and
  restores the Seq ID / Strand / Identity legends.

- [`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
  now warns about suspicious input data: reversed or out-of-range
  alignment/gene coordinates, `pident` outside \[0, 100\], and sequence
  IDs that are not present in `seq_data`.

- [`geom_gene_label_repel()`](https://dangjem.github.io/ggchord/reference/geom_gene_label_repel.md)
  gains `gene_label_side = "auto" | "inside" | "outside"`. With
  `"outside"`, labels that would sit inside the chord (where they can
  overlap the ribbons) are mirrored to the outside of their sequence
  arc, keeping the same radial distance from the arc.

- New `gene_label_segment_linetype` argument controls the leader-line
  linetype. The default `"auto"` draws solid lines, except for labels
  that were moved to the other side of their arc (`gene_label_side`),
  which are drawn dashed. Any other valid ggplot2 linetype
  (e.g. `"dotted"` or a numeric dash pattern) is applied to all leader
  lines.

- Elbow leader lines no longer force fixed segment lengths: the stub
  scales with each label’s text width and the horizontal space available
  between the gene and the label, so labels can be placed more flexibly
  without degenerate (near-zero) stubs.

- [`geom_seq_label()`](https://dangjem.github.io/ggchord/reference/geom_seq_label.md)
  now documents and follows the intended `seq_label_radius` semantics:
  `1` sits on the arc, `> 1` places the label outside (away from the
  chord center) and `< 1` inside. Previously the multiplier was applied
  in the opposite direction (the default `1.15` put labels *inside* the
  chord).

- New
  [`geom_seq_label()`](https://dangjem.github.io/ggchord/reference/geom_seq_label.md)
  options: `seq_label_orientation = "arc" | "horizontal"` (horizontal
  labels extend away from the chord center), `seq_label_hjust` /
  `seq_label_vjust` for per-sequence justification, and `check_overlap`
  to skip labels that would overlap.

- The default theme no longer draws grid lines (`panel.grid` is blank)
  and legend keys are transparent (they blend into the plot background
  instead of a fixed white rectangle).

### Performance

- Replaced the linear angle lookup in the layout mapping with a binary
  search (`findInterval`), speeding up layout computation for large
  plots.

### Dependency changes

- Declares `ggplot2 (>= 4.0.0)` and `R (>= 4.1.0)` to match the
  implementation (the package relies on ggplot2 4.x internals).

### Infrastructure

- Added a GitHub Actions `R CMD check` workflow (macOS, Windows, Linux).
- Removed the internal legacy `fill_ggnewscale_1` aesthetic name in
  favour of `fill_ribbon`.

### Bug fixes

- Tests no longer write to a hard-coded `/tmp` path: they use
  [`tempfile()`](https://rdrr.io/r/base/tempfile.html), so the test
  suite passes on Windows and leaves no stray files behind for
  `R CMD check` (fixes the CRAN incoming-check failure).

- The Identity(%) colourbar no longer collapses into a thin/invisible
  line when the legend is placed at the top/bottom or the legend box is
  horizontal (`legend.box = "horizontal"`). It now follows the theme’s
  legend position: a vertical bar filling the available height at the
  left/right, and a fixed-size horizontal bar at the top/bottom.

- Legend keys are transparent and do not inherit `panel.background`
  (ggplot2 4.x lets unset legend keys follow the panel background, so
  the key fill is set explicitly to stay transparent).

- [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)
  output now shows the Seq ID / Strand / Identity legends (the
  layout-level `showlegend` switch is enabled) and reproduces the
  [`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md)
  directional arrowheads as plotly annotations.

- Sequence (and gene) labels no longer end up upside down when a global
  `rotation >= 90` is used: the readability flip is now re-applied after
  the layout rotation instead of only before it.

- The repulsion spring now pulls labels toward their own starting
  position rather than the leader-line anchor, which keeps labels moved
  with `gene_label_side = "outside"` on the outside while their leader
  line still starts at the gene.

- With `gene_label_side`, every label is kept on the requested side of
  its arc (previously only the labels moved by the side switch were
  re-checked, so a crowded repulsion layout could push other labels
  across the arc).

- The built-in `gene_data_example` annotations no longer contain
  URL-encoded `%2C` artifacts (e.g. “ribonucleotide reductase%2C large
  subunit” is now “ribonucleotide reductase large subunit”).

### New features

- Each legend can now be positioned independently via the
  `legend_position` argument of
  [`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md),
  [`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md)
  and
  [`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
  (e.g. `geom_ribbon(legend_position = "bottom")`). Legends without an
  explicit position stay together at `theme(legend.position = ...)`.

- Parameter specification is now more flexible and human-friendly.
  Sequence parameters accept a single value, vectors, vectors/lists
  named by sequence ID, lists named by sequence order (`"1"`, `"2"`, …)
  and unnamed lists; gene parameters additionally accept per-strand
  (`+`/`-`) specifications in any of those forms
  (e.g. `gene_label_rotation = c("+" = -15, "-" = -45)` or
  `list(c("+" = -15, "-" = -45), ...)`), including length-one lists that
  recycle.

## ggchord 0.5.0

### New features

- Added ribbon outline customization to
  [`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md):
  `ribbon_outline_color` (default `"black"`), `ribbon_outline_width`
  (default `0.05`) and `ribbon_outline_linetype` (default `1`, solid).

### Dependency changes

- Removed the `ggnewscale` dependency. The ribbon and gene fill scales
  are now kept independent via an internal renamed-fill aesthetic, so no
  external package is required for plots with both ribbon and gene
  layers.
- Removed the `RColorBrewer` dependency. The default Set1 categorical
  palette is now built into the package, so the rendered default colors
  are unchanged.

### Bug fixes

- Fixed the ribbon fill scale being overwritten by the gene fill scale
  when both
  [`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md)
  and
  [`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
  were present (previously produced wrong ribbon colors and a “Scale for
  fill is already present” message).
- Fixed `ribbon_alpha` rendering at the wrong opacity (e.g. `0.35` was
  drawn as ~0.55); the alpha value is now applied exactly as specified.
- Fixed `geom_axis(show_axis = FALSE)` failing with “object ‘label’ not
  found”.
- Fixed `axis_label_orientation` rejecting mixed vectors such as
  `c("horizontal", 45, ...)`.
- Fixed warnings from `brewer.pal()` when fewer than 3 sequences or gene
  annotations are used (two-sequence plots now render cleanly).
- Fixed an error when
  [`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
  was added before
  [`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md)
  (“Continuous value supplied to a discrete scale”).
- Fixed plots containing only
  [`geom_axis()`](https://dangjem.github.io/ggchord/reference/geom_axis.md)
  where the axis path was misclassified as a sequence arc.
- Registered `+.ggchord` and `ggplot_build.ggchord()` as proper S3
  methods and aligned the `ggplot_build.ggchord()` signature with the
  generic.

### Documentation

- Translated all code comments and user-facing messages to English.
- Added man pages for previously undocumented exported functions
  ([`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md),
  [`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md),
  [`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md),
  [`geom_axis()`](https://dangjem.github.io/ggchord/reference/geom_axis.md),
  [`coord_chord()`](https://dangjem.github.io/ggchord/reference/coord_chord.md),
  `+.ggchord()` and others).
- README: added per-column data preparation tables with example rows,
  rendered example plots under `examples/plots/`, generalized the
  package description beyond BLAST, and documented the ribbon outline
  parameters.
- Rewrote the package vignette for the layered `ggchord() + geom_*` API.

## ggchord 0.4.0

### Changes

- Parameter redistribution: layout parameters moved from
  [`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
  into the individual `geom_*` layers;
  [`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
  now only validates data and sets global style (`title`, `rotation`,
  `panel_margin`, `show_legend`, `debug`).
- Deferred computation: the coordinate layout is computed at
  [`print()`](https://rdrr.io/r/base/print.html) time, collecting
  parameters from all layers during rendering.
- Custom `print.ggchord()` method: merge parameters, compute the layout,
  inject data into layers, then render.
- Added 15 unit tests.

## ggchord 0.3.0

- Layered API refactoring: split the monolithic function into
  `ggchord() + geom_seq() + geom_ribbon() + geom_gene() + geom_axis()`.
- Custom `+.ggchord` method that flattens layer lists automatically.
- Lightweight
  [`coord_chord()`](https://dangjem.github.io/ggchord/reference/coord_chord.md)
  coordinate system.

## ggchord 0.2.0

CRAN release: 2025-07-16

- Enhanced arc and line mode optimization.
- Precise curvature and gap control.
- Enhanced color customization.

## ggchord 0.1.0

- Separate management of sequence, alignment, and gene data.
- Sequence orientation, custom order, gap and radius adjustment.
- Customizable axes; ribbons support 3 coloring schemes.

## ggchord 0.0.2

- Multi-sequence support; arc/line mode switching.

## ggchord 0.0.1

- Initial release: pairwise alignment chord diagram visualization.
