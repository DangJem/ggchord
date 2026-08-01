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
