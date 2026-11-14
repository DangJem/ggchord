Compare arrow outlines on the same observed CDS. Choose shape before size; neither changes the interval, strand or biological function.

## Prepare a real annotation and local window {#input}

Use the bundled integrase in `gene_data_example`: OQ646790.1 at 43,765–46,140 bp. Restore the original positive strand from `source_strand`. Focus on 43,000–47,000 bp for visibility and display source coordinates on the axis.

{{code:comparegenedata}}

`gene_window` retains clipping and source metadata. Your CDS needs a complete sequence length, matching `accver`, one-based inclusive `start/end`, and an evidence-backed strand. Use `geom_feature()` for unknown strands. This page runs independently in a new R session and requires the optional dplyr package.

## Select the real neighbourhood by explicit ID {#locus-selection}

Combine source accession/endpoints into a workflow target ID, not a claimed database gene ID. Find_ggchord_loci selects exact IDs and extends 100 bp upstream/downstream in biological direction. Negative-strand upstream points to larger source coordinates; unknown/bidirectional strands require symmetric extensions. Pass returned windows to focus_ggchord_data and inspect full-length boundary clipping flags.

{{code:comparefindloci}}
{{console:comparefindloci}}

## Three arrow constructions {#arrow-heads}

All figures use 9 × 4 inches and identical input, window and layers. To expose shoulders, all three explicitly use `arrow_head_width = 1.6` (default 1). Only `arrow_head_style` changes.

{{code:comparegeneshouldered}}
{{figure:comparegeneshouldered}}

`shouldered` is the default construction. A wider head creates shoulders where it joins the shaft, emphasizing direction.

{{code:comparegeneflush}}
{{figure:comparegeneflush}}

`flush` uses a flush head-to-shaft outline. Compare the head junction while keeping the same endpoints and shaft.

{{code:comparegenetriangle}}
{{figure:comparegenetriangle}}

`triangle` draws the entire interval as a triangular wedge rather than adding a head to a rectangular shaft. Shape differences do not imply different functions or measured values.

## When a short feature cannot fit its head {#short-features}

Keep the same observed integrase and window, requesting head length 1.5 in all three plots so the interval is short relative to its head. Change only short_feature: auto chooses a compressed directional glyph, wedge a wedge, block an undirected block. Do not lengthen source coordinates to fit a head; a block does not make the biological strand unknown.

{{code:compareshortauto}}
{{figure:compareshortauto}}
{{code:compareshortwedge}}
{{figure:compareshortwedge}}
{{code:compareshortblock}}
{{figure:compareshortblock}}

## Dimensions, short features and presets {#parameters}

| Parameter | `geom_gene()` default | Values and choices |
| --- | --- | --- |
| `width` | `NULL`, resolving to 0.05 | Shaft width; may vary by sequence or strand |
| `arrow_head_length` | `0.04` | Head length in the local frame, not bp |
| `arrow_head_width` | `1` | Head width relative to the shaft; 1.6 in all three figures |
| `arrow_head_style` | `"shouldered"` | shouldered / flush / triangle |
| `short_feature` | `"auto"` | Selects wedge or block when the head cannot fit; explicit wedge/block fixes the fallback |
| `position` | `"identity"` | position_strand separates strands; position_feature_stack separates overlaps |

Do not hide short intervals by enlarging heads. Check source length before choosing a fallback. Generic `geom_feature()` accepts fixed `shape` or mapped `aes(feature_shape = ...)` with a corresponding scale; see [styling](styling.html).

`geom_feature_plasmid()` is a preset with different width, head dimensions and Position defaults. Changing several defaults while comparing it to gene is not a single-parameter experiment. Segmented annotation boundaries mark segments, not different genes.

## Modes, input adaptation and next steps {#adapt}

Gene arrows support bp inputs in chord, collinear and single circular modes. Use generic features for aa intervals. Focusing and display reversal do not rewrite biological strand; see [compatibility](compatibility.html).

First change the common head width from 1.6 back to 1 and inspect shoulder visibility. Then change only `width` and compare shaft versus backbone width. Keep export dimensions fixed rather than editing source intervals to improve appearance.

Continue with [feature stacking](annotations.html#stack), [labels](labels.html), or add this CDS to the [four-genome example](collinear.html). Full interfaces: [`geom_gene()`](reference/geom_gene.html), [`geom_feature()`](reference/geom_feature.html).
