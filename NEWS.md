# ggchord 0.13.0

* Refined the pBluescript II SK(+) circular-map layout against the supplied
  visual reference. Restriction labels again retain their close polar contour
  and local ordered fans instead of being flattened onto shared left/right
  rails with feature callouts. Feature labels now exhaust internal feature
  gutters and deeper label tracks before becoming external, with short leaders
  for materially displaced internal labels. Compact, promoter and primer
  arrows use the same annular-arrow factory as other directional features,
  the same requested thickness, and their real genomic interval; the common
  short-arrow rule alone contracts glyphs that cannot contain a normal head.
  Short dotted/dashed segment boundaries are emitted as explicit broken
  geometry so graphics devices cannot silently render them solid. Plasmid axis
  ticks are longer, labels remain inside the backbone, and the automatic major
  tick target is now five.
* Enlarged circular feature-band separation into text-bearing corridors. Each
  corridor is subdivided into as many fixed-radius label tracks as its measured
  height permits. Feature labels select a radius first, then derive their
  position and tangent from that circle; labels placed on a deeper radius use a
  leader back to the owning outer feature instead of drifting around the map.
* Added `coord_circular()` as an independent coordinate contract for one
  circular sequence. It supports a closed circle (`gap = 0`), a degree-based
  opening, genomic-origin rotation, clockwise/counterclockwise direction,
  position scales, label-aware fitting, manual limits and clipping. Ribbon and
  link layers fail clearly because they do not belong to this coordinate.
* **Breaking Position change.** Before v0.13, `geom_gene()` implicitly used
  strand-separated placement with `gene_offset = 0.1`. Since v0.13,
  `geom_gene()` defaults to `position = "identity"`. To reproduce the former
  default, use `geom_gene(position = "strand")`. For plasmid-style shared-band
  placement, use `geom_gene(position = "plasmid")`.
* Added `position_strand()` and `position_plasmid()`. Signed offsets use the
  sequence-local outward normal: `position_strand(offset = 0.1)` exactly
  reproduces the former default (`+` strand outside, `-` strand inside), while
  named `+`/`-` values are already signed and are not flipped again.
  `position_plasmid(offset = -0.1)` puts both strands on the same inner band
  and leaves arrow direction controlled by strand.
* Feature Positions accept the standard ggplot2 identity spellings
  `position = "identity"`, `position = position_identity()`, and
  `position = PositionIdentity`. The latter two ggplot2 objects are re-exported
  for use after `library(ggchord)`.
* `position_feature_stack()` now accepts `base_position = "identity"`,
  `"strand"`, `"plasmid"`, or the corresponding Position object. Its no-argument
  preset uses wider 0.10-unit lanes for dense plasmid annotations. Layout export now
  records `position_name`, `base_offset`, `lane`, `lane_offset`, and
  `normal_offset` for gene and feature polygons.
* Explicit `gene_offset` and `feature_offset` remain temporarily accepted,
  translate through the former flexible accver/strand resolver, and emit a
  migration warning. Combining a legacy offset with a non-identity Position is
  an error.
* Refactored `geom_gene()` into a gene-specific convenience layer over the same
  neutral feature geometry pipeline as `geom_feature(shape = "arrow")`.
  Circular splitting, local tangent/normal calculations, arrow heads and
  short-feature fallback are implemented once. Gene data roles, gene fill
  scales, strand guide and key glyph remain independent public behavior.
* `geom_feature()` and `geom_gene()` now share `arrow_head_length`,
  `arrow_head_width`, `arrow_head_style = "shouldered" | "flush" |
  "triangle"`, and `short_feature = "auto" | "wedge" | "block"`.
  Features crossing the circular origin retain one source identity and split
  into drawable arc pieces without losing strand direction.
* Added `geom_feature_plasmid()` as a compact circular-map preset over the
  generic feature engine. Multi-segment records remain one biological feature
  during stacking and labelling; visible segment runs share one track and one
  label, preserve gaps and per-segment colours, and draw arrowheads only at the
  biological feature ends. Segment membership alone no longer draws a dashed
  divider; only source cleavage arrows or an explicit non-solid `line_style`
  create an internal boundary.
  Its default position assigns overlapping intervals to generic collision
  slots: explicit `display_priority`/`prioritized_display` comes first,
  structural spans stabilize the outer bands, and short local annotations
  follow genomic order so a non-overlapping cassette can continue on one lane.
  Feature type, name, and colour no longer choose a radial lane, and
  `find_common_features()` no longer emits semantic `preferred_lane` hints.
  The feature-shape vocabulary includes `compact_arrow`,
  `promoter_arrow`, `primer_arrow`, and `marker`. The three directional arrow
  variants share one outline and width algorithm. Default common-feature geometry now follows the
  feature's own nondirectional/forward/reverse/bidirectional value rather than
  its biological type. Extremely short directional intervals retain a
  restrained direction mark whose head and thickness contract with available
  display length; zero-length point features draw a radial tick rather than a
  minimum-width polygon. Small glyphs receive a lighter outline unless
  linewidth is explicit. Lane collision checks use the true rendered interval
  rather than a type-specific fabricated span. The plasmid preset now uses the visual
  reference's semantic feature colours and type-aware default shape hints while
  retaining user scale and explicit geometry priority.
* `geom_seq()` adds `seq_style = "single" | "double" | "band"`,
  `seq_backbone_gap`, and `seq_backbone_width` for circular backbones. Under
  `coord_circular()`, the implicit chord-direction arrow and redundant
  one-sequence legend are omitted; explicit `arrow` and `show.legend` values
  still take priority.
* Added `geom_feature_label()` and `geom_feature_label_repel()`. They retain a
  feature-specific public vocabulary while reusing the existing device-aware
  text measurement, wrapping, ellipsis, collision, leader routing and frame
  fitting engine. The default `label_layout = "feature"` compares measured text
  width with usable feature arc length, uses tangent text inside intervals that
  can contain it, and places compact-feature text nearby. Nearby collisions are
  nudged deterministically and only materially displaced labels receive a light
  leader; explicit radial/auto/callout layouts remain available. The feature
  mode now exports `feature_label_mode = "inside" | "adjacent" | "external"`,
  measures the final font family, face, size and line height, and checks
  oriented label boxes against other labels, rendered feature polygons, the
  centre reserve and the circular backbone. It searches feature-band gutters
  and deeper internal label tracks before unresolved labels become horizontal
  external callouts; setting
  `external = FALSE` hides them instead of pushing them indefinitely towards
  the centre. Label-aware coordinate fitting includes their measured boxes.
  Feature thickness responds to final label size within bounded limits, never
  to unbounded label length. Inside text automatically uses black or white from
  the final fill luminance unless its colour is explicitly supplied. External
  feature labels use rounded pastel callouts derived from their resolved fill.
  Restriction-site polar placement remains independent and is never rewritten
  into shared Cartesian side rails by the presence of a feature callout.
* Added `find_restriction_sites()` as a calculation-only API. It preserves one
  row per pattern match, stable pattern identity, multiple motifs per enzyme,
  overlapping/IUPAC/reverse-complement matches, circular-origin matches,
  negative and out-of-motif cleavage offsets, Type IIS, 1/2/4-cut records and
  unknown cleavage (`ncuts = 0`). Custom motifs may be a data frame or named
  character vector. Overlapping fixed anchors inside degenerate IUPAC motifs
  are now enumerated without consumption, fixing missed sites such as
  EcoO109I in pBluescript II SK(+). REBASE source motif casing is retained for
  auditability and lower-case reverse-cleavage definitions are oriented before
  calculating display positions.
* Added `filter_restriction_sites()` for display subsets (`all`, `unique`,
  `unique_dual`, `six_plus`, `unique_6plus`, and `commercial`) and intersecting
  enzyme, motif-length, site-count, window and commercial filters. Filtering
  never rewrites biological coordinates. `parent_set` additionally selects all
  enzymes, all commercial enzymes, or a stable nonredundant commercial subset;
  equivalence groups and preferred representatives come from the bundled
  REBASE `embossre.equ` source.
* `geom_restriction_site()` draws ticks, labels and independent radial leaders
  in one deterministic layer. Real site anchors are separated from final
  label positions. All enzyme-facing text edges follow one close circular
  contour in genomic order. Sparse callouts retain their natural angles, while
  dense lateral clusters switch deterministically to a compact ordered fan
  with uniform measured spacing rather than a rigid Cartesian text wall.
  Natural sparse callouts use a direct connector; displaced callouts
  use exactly two segments, an independent radial stub and fan connector.
  Enzymes sharing one cleavage coordinate are combined
  into one stable callout while all contributing rows remain in `source_rows`.
  Each connector ends on the enzyme-bearing left or right edge of the measured
  label box: right-travelling lines meet its left edge and left-travelling
  lines meet its right edge. Near-vertical leaders use the nearest corner on
  that same side rather than switching to the top/bottom edge. Mirrored text
  order keeps the enzyme at the connected side. The former `leader = "trunk"`
  spelling remains
  only as a compatibility alias and no longer creates a shared trunk.
  Dense-fan leaders retain individual site anchors, use a short radial root,
  and then fan out in genomic order from a narrow bundle. Default composite
  labels are rendered as two coordinated grobs. Only unique cutters use bold
  enzyme names; repeated cutters and coordinates use regular weight. Adjacent
  but distinct bp sites are never merged.
* Added `geom_seq_center_label()`, `theme_ggchord_plasmid()`,
  `scale_feature_fill_plasmid()`, and `scale_feature_shape_plasmid()`. Manual
  feature scales take priority over presets regardless of addition order. The
  plasmid theme places coordinate ticks and labels inside the backbone. The
  automatic circular-origin tick is longer and does not repeat a `0` label.
* Added reproducible one-row pUC19, pBR322 and pBluescript II SK(+) sequence
  fixtures generated byte-for-byte from unchanged FASTA files. Added a bundled
  normalized common-feature database with provenance metadata and dynamic
  exact DNA and six-frame protein matching. In `mode = "auto"`, records now
  follow their stored detection mode: `exactProteinMatch` uses exact protein
  followed by near-exact DNA, while other records use exact DNA. Audited
  annotations embedded in thirteen reference `.dna` files provide an exact
  known-sequence layer for pBR322, pUC19, pBluescript II SK(+), pSB1C3,
  pET-28a(+), pETDuet-1, pcDNA3.1(+), pTRE-Tight-BI,
  pSpCas9(BB)-2A-GFP (PX458), pDONR221, pCAMBIA1300, pEarleyGate 201, and
  pTRIPZ. The generated database currently contains 201 biological features and
  218 segments, plus seven true primer binding-site records from pSB1C3 and
  pETDuet-1 in a separate internal reference table rather than the Feature
  table. It records each source SHA-256 and record count. The former
  `plasmid_example_pUC19c` object remains as a pre-release compatibility alias
  of `plasmid_example_pUC19`. Complete REBASE 609
  pattern data are embedded for working-directory-independent restriction-site
  search and remain reproducible from `VERSION` plus the three `embossa_*.txt`
  source files and `misc/embossre.equ`. Other files in the extracted
  `examples/rebase/misc/` directory remain supplemental verification sources;
  runtime searches do not depend on that directory.
* Removed the v0.12 `seq_id` data/mapping compatibility entry. Rename it to
  `accver` before calling ggchord functions.
* Removed `ggchord::geom_ribbon()` after its v0.12 deprecation. Use
  `geom_link_ribbon()`; existing ribbon aesthetics, scales, themes, and layout
  export component names are unchanged.
* `geom_link_ribbon()` now defaults to `link_avoid = "uniform"`, using one
  obstacle-aware clearance across each endpoint. Local `"smooth"` avoidance
  and fixed `"none"` spacing remain explicit choices; an explicit
  `ribbon_gap` still takes priority.

# ggchord 0.12.0

* Added `geom_link_ribbon()` as the canonical interval link geometry. Ribbon
  aesthetics, scales, themes and exported layout components keep their names.
* `ggchord::geom_ribbon()` now emits a deprecation warning directing callers
  to `geom_link_ribbon()`; removal is planned for v0.13.0. Ribbon statistics
  use the canonical constructor without triggering the compatibility warning.

* Single-sequence tables now use `accver`. Legacy `seq_id` columns/mappings
  warn in v0.12 and will be removed in v0.13; supplying both is an error.
  Packaged example data have been regenerated from unchanged source files.
* Ribbon drawing requires only the six endpoint columns. `length` and
  `pident` are optional; requested filters/statistics still require their
  inputs. Missing identity uses a fixed fill without an Identity guide.
* Gene annotations and feature categories may be omitted. Coordinates,
  accessions, sequence lengths and feature strands remain required.
* Added point links (`geom_link_line()`), independent link scales and
  `legend.link.*`, and re-exported `grid::unit()` and `grid::arrow()`.
* Entity obstacle providers now use rendered shapes. Local clearance is opt-in
  through `link_avoid = "smooth"` or `"uniform"`; the default is `"none"`.
* `link_branch = "query"` or `"subject"` optionally shares equal endpoints
  with equal mapped styles. Trunks retain source membership in layout exports.
* Removed `geom_ribbon_highlight()`. Use an additional
  `geom_link_ribbon(data = subset, fill = ...)` layer instead.
* Fixed feature-shape legend mappings and increased gene-arrow key height.
* Shared-link fallback now preserves vertex-wise staged aesthetics instead of
  replacing each unshared path with its first vertex style.
* Shared ribbon branches now honor explicit edge control points, including
  layout rotation, while preserving the shared junction tangent.
* Split the constructor/build, automatic-scale, frame-fitting, sequence,
  ribbon, gene/region, axis, and label calculations into focused internal
  source modules without changing the public API or layout output.
* Clarified why direct printing in a small RStudio Plots pane can differ from
  export-size rendering, and recommends `view_ggchord()` for composition with
  an explicit matching `ggsave()` size for final files.
* Archived previous-release design records in `design-roadmap/DESIGN-HISTORY.md`.

# ggchord 0.11.0

* `geom_gene_label_repel()` now defaults to `gene_label_layout = "radial"`:
  horizontal names follow each sequence's offset contour, with normal leaders
  and short normal departures for displaced labels. `"auto"` combines side
  columns with radial placement elsewhere, replacing the removed `"aligned"`.
* Sequence curvature is continuous through zero and one and accepts signed
  values without silent clipping; opposite signs mirror the bow across the
  endpoint chord. Very large bows can make tracks intersect geometrically.
* Static previews use a wider default and fit their height to content when
  `height = NULL`; explicit output sizes remain respected.


* Radial packing explores a wider, ordered fan before adding an outer contour,
  and balances the extra clearance of opposite top/bottom sectors.
* Content-fitted previews now export the measured gtable directly, avoiding
  height-dependent relayout oscillation and excessive outer legend whitespace.
  Default plot margins retain a small physical safety inset.

## Preview and label fixes

* Nested text measurements close only their own graphics device and restore
  the previously active device. Alternating printing and previews no longer
  changes preview layout or closes the IDE device, leaving later plots hidden.
* Radial label contours now limit tangent extension past sequence endpoints.
  Crowded labels try an outer contour instead of drifting toward a neighbouring
  sequence, including mixed-radius plots with sequence names. The same bound
  applies to the radial portion of `"auto"`.
* Removed `"aligned"` requests now point directly to `"auto"`. Acceptance
  checks verify that auto actually combines side columns with radial labels.
* Static previews retain export warnings rather than silently suppressing
  font, device or missing-data problems. Content fitting warns when legends
  exceed the requested width; users can increase width or set guide rows.
  A regression checks real PNG
  dimensions and preservation of the caller's graphics device.
* `tools/validate-release.R` provides reproducible common-geometry, label,
  output, benchmark and executable-vignette checks, writing only temporary
  artifacts. Extreme parameter stress tests are deferred.
* Tutorial code uses the current API; retained figures are explicitly marked
  as historical. Obsolete pkgdown references and preview-size descriptions
  are corrected without rebuilding the published site or its assets.

## Advanced tracks and layout

* New `stat_ribbon_bundle()` and `stat_ribbon_density()` expose
  `after_stat(bundle_n)`, `after_stat(bundle_weight)`, and
  `after_stat(density)` while leaving `geom_ribbon()` conservative.

* New `focus_ggchord_data()` synchronizes sequence, gene/feature, and ribbon
  cropping across one or more loci. `position_feature_stack()` assigns the
  minimum deterministic radial lanes, while `seq_ring` and
  `scale_seq_ring_manual()` provide explicit multi-ring radii.


* `geom_gene_label_repel()` adds `gene_label_segment_overlap = "fade" |
  "clip" | "show"` and `gene_label_segment_overlap_alpha`. Covered leader
  portions now default to a faint continuous line instead of a hard gap;
  clipping and full display remain explicit per-layer choices.

* `geom_gene_label_repel()` gains `gene_label_fit` and
  `gene_label_max_lines`. Initially crowded annotations can be wrapped,
  shortened with an ellipsis, handled automatically, or left unchanged;
  explicit `gene_label_wrap` remains the fixed-width override.

* `geom_gene_label(gene_label_side = "auto")` now aligns horizontal text from
  its displacement relative to the local sequence curve rather than the plot
  origin, fixing inward labels in several quadrants. Single-strand gene plots
  also generate a Strand guide with the correct number and direction of keys.

* The packaged gene example is rebuilt deterministically from the unchanged
  source table with eight well-spaced, visibly varied features per sequence
  and balanced display strands. `source_strand` retains the source annotation
  direction, and the complete rule is recorded in
  `data-raw/generate_example_data.R`.

* The packaged ribbon example is rebuilt from the unchanged BLAST files with
  a 300-base minimum and at most three spatially distributed alignments per
  dense sequence pair. Sparse pairs remain intact, reducing the fixture from
  31 to 13 ribbons without fabricating alignments.

## Standard ggplot2 layer grammar

* Every public geom and stat now returns one standard `LayerInstance`.
  `geom_gene_label_repel()` uses a composite Geom/gTree implementation instead
  of returning child-layer lists. The custom `+.ggchord` list-flattening path
  has been removed.

* Public layers consistently expose `mapping`, `data`, `position`,
  `show.legend`, `inherit.aes`, and `...`. A layer `data` argument may be a
  data frame or a function receiving that layer's default input table.

* `prepare_ggchord_plot()` is now the single orchestration path for geometry,
  scales and coordinates. `ggplot_build.ggchord()` delegates to it before the
  standard ggplot2 build, keeping printing, `ggsave()` and `view_ggchord()`
  consistent.

* `coord_chord()` now constructs a real `CoordChord` ggproto. Global
  most-recent-layout state has been removed; `get_chord_layout()` requires an
  explicit ggchord plot.

* Role aesthetics use tidy evaluation once per build and preserve source
  columns. Geometry roles reject `after_stat()`/`after_scale()`, while visual
  mappings continue to support staged aesthetics. User scales always win;
  missing visual scales are inferred as continuous or discrete, and
  incompatible mappings fail clearly.

## Themes and namespace

* `theme_ggchord()`, `theme_ggchord_minimal()`,
  `theme_ggchord_publication()`, and `theme_ggchord_dark()` are now complete,
  theme-like appearance interfaces with dotted arguments. `text` supplies the
  parent typography and child sizes inherit through `rel()`; the redundant
  `base_size`, `base_family`, underscore element arguments, and
  `theme_ggchord_elements()` have been removed.

* Sequence axes are displayed automatically and configured with `axis.*`
  theme arguments. `geom_axis()` has been removed. Axis gaps, major/minor tick
  lengths and text offsets use physical `grid::unit()` values; breaks, labels,
  limits and transforms remain in `scale_seq_position_continuous()`.

* Sequence, ribbon, gene, feature and region guides can independently inherit
  or override position, direction, background, key dimensions, text, title,
  margin and spacing through `legend.<role>.*`. Ribbon colourbar ticks and its
  axis line are independently themeable. Explicit scale/`guides()` settings
  still take precedence, and feature guides no longer inherit gene placement.

* The default Identity colourbar is placed at the left edge of the panel;
  compact sequence, strand, feature and region guides use the right edge.
  `view_ggchord()` defaults to an 11-inch-wide content-fitted preview. Explicit
  common or role guide positions and explicit preview dimensions still win.

* British colour spelling is canonical. `seq_color`, `ribbon_color`, `colors`,
  `scale_*_color_*()` and the new `guide_ggchord_colorbar()` are complete US
  aliases; supplying both spellings in one call now fails explicitly.

* Layout-only operations no longer open an implicit `Rplots.pdf` device while
  converting physical units. This also fixes a following `ggsave()` or
  `view_ggchord()` occasionally closing the wrong graphics device.

* Selected ggplot2 helpers such as `theme()`, `ggtitle()`, `labs()` and
  `guides()` remain exact re-exports. ggchord does not wrap them or expose
  unrelated Cartesian geoms, facets, coordinates, scales, or complete themes.

## Breaking API cleanup

v0.11.0 removes the previous soft-compatibility layer. Removed arguments fail
immediately with a migration message instead of warning or being ignored.

| Removed entry | Replacement |
| --- | --- |
| `ggchord(title/rotation/panel_margin/show_legend)` | `labs()` / `coord_chord()` / `theme()` |
| `show_legend` | `show.legend` |
| geom legend position/key arguments | `guides()` + `guide_ggchord_*()` |
| geom colours, order, limits, breaks and names | role-specific `scale_*()` |
| ribbon `*_by`, direction and visual-value arguments | `aes(ribbon_* = ...)` + scales |
| `ribbon_alpha` and ribbon outline names | `alpha`, `colour`, `linewidth`, `linetype` |
| feature type/category/label column strings | `aes(feature_type/feature_fill/feature_label = ...)` |
| `geom_seq(seq_labels/seq_colors)` | `geom_seq_label(labels = ...)` / sequence scale |
| `geom_axis()` and its style/layout arguments | automatic axis + `theme_ggchord(axis.*)`; breaks/labels use the position scale |
| `geom_seq_region(regions/region_* style)` | `data`, `fill`, `colour`, `alpha` |

Vignette code examples are migrated to the current API and checked separately
from document rendering. Existing figures and the public site are retained
for a separate documentation rebuild after this release.

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
  rotation and explicit rings. They share fixed-obstacle avoidance,
  cross-sequence collision handling, order-preserving leader routing and
  device-aware clipping. Rotated labels use oriented-rectangle collision and
  clipping, preventing spurious whitespace and oversized breaks in leaders.

* The layout ideas adapt orderly multi-sequence callout figures and the
  external/inside feature-label approaches offered by SnapGene and Geneious.
  ggchord uses generic API names and integrates the ideas into its own geometry
  system. See the
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
