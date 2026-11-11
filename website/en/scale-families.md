Choose the meaning of data before roles and scales. Use the same observed HSPs, annotations and GC windows for continuous values, categories, identity, opacity, widths, line types and themes, grouping related helpers on one page.

## Independent input and scale responsibilities {#input}

Independently load ggchord and optional dplyr. Complete lengths, HSP identity and query/subject IDs come from bundled source records. Named vectors assign display values to existing categories, not new samples or measurements. Map inside aes, fix styles in geoms, and let explicit scales override automatic defaults. Explicit fixed aesthetics override matching mappings; omit fixed alpha when mapping opacity.

{{code:comparescaledata}}


## Continuous, stepped and direct colour values {#numeric}

Both map observed pident to ribbon_fill with fixed 0–100 domain, palette and dimensions. Gradientn interpolates continuously; stepsn uses explicit breaks. Matching colours do not establish new homology. Values place palette anchors, not source pident; breaks control ticks/steps and limits the domain. Check NA and out-of-bounds policy rather than treating missing values as zero.

{{code:comparescalegradient}}
{{figure:comparescalegradient}}
{{code:comparescalesteps}}
{{figure:comparescalesteps}}


## Identity preserves supplied colour values {#identity}

First derive display colours from observed identity using an explicit 90% threshold, then preserve those strings with scale_ribbon_fill_identity. This is declared display classification, not a way to recover continuous identity. Use a gradient for a continuous legend. Identity normally hides its guide; provide deliberate breaks, labels and guide when one is needed.

{{code:comparescaleidentity}}


## Category fill, opacity, outline and line type {#categories}

The same qaccver receives the same display values. Name manual values by actual categories to avoid order-dependent mismatches. This combines independently trained roles rather than claiming a single-parameter visual experiment. Use `scale_ribbon_alpha_continuous(range=...)` for numeric opacity and `scale_ribbon_alpha_manual` for categories. Colour and linetype change outlines, not alignment direction.

{{code:comparescalemanual}}


## Independent scales for point links {#point-links}

Explicitly calculate midpoints of observed HSPs; point links do not display full spans. Manual helpers map link_colour, link_alpha, link_linewidth and link_linetype by category; gradientn, alpha_continuous and linewidth_continuous map numeric values. Construct category and numeric combinations separately.

Automatic numeric pident mappings use the 0–100 domain and six Identity ticks: 0, 50, 80, 90, 95 and 100. Explicit scales keep their own settings. The example sets the same domain, ticks and title for all three channels so their legend keys can combine.

{{code:comparescalelinkdata}}
{{code:comparescalelinkmanual}}
{{code:comparescalelinkcontinuous}}
{{figure:comparescalelinkcontinuous}}


## Feature, gene and region roles {#features}

Map feature_fill and feature_shape to the same type to merge their guide silhouettes; request a separate shape guide explicitly. Fixed shape overrides its mapping. Mapping all types to arrow demonstrates mapping without inventing functional shape semantics. Plasmid fill/shape presets use their documented category rules.

{{code:comparescalefeature}}
{{code:comparescalepreset}}


## Gene strand, sequence rings and regions {#rings}

Restore biological strands from source_strand. Region, gene and ribbon fill roles are independent. `seq_ring` requires explicit positive radii, never inferred track counts; it conflicts with seq_radius and collinear mode. Sequence radius mapping is distinct from quantitative position_track bands.

{{code:comparescalegenering}}
{{code:comparescaleregion}}


## Named track domains and colours {#tracks}

Compute windows from actual pUC19 DNA. `track_value` controls height; track_fill/track_colour train separate colour scales, with IDs matching position_track. The 0.5 threshold is a display choice and full fractions remain. Tiles do not average data; NA heights do not become zero. Declare scales for each ID independently.

{{code:trackdata}}
{{code:comparescaletrack}}
{{code:comparescaletrackgradient}}


## Guides, spelling aliases and staged mapping {#guides}

`guide_ggchord_legend` manages discrete role guides and guide_ggchord_colourbar continuous colour bars. `show.legend` controls participation, guides chooses guides, and theme controls position/direction/appearance. Colour/color aliases are equivalent: guide_ggchord_colorbar, scale_seq_color_manual, scale_ribbon_color_manual, scale_link_color_manual and scale_link_color_gradientn correspond to colour versions. Do not add both synonymous scales expecting cumulative control.

`after_stat` reads computed variables such as bundle_n in the next lesson; after_scale computes supported visual attributes after scaling. They come from ggplot2 and have distinct staged semantics even when their R objects compare equal. All ggplot2/grid reexports follow upstream help, with ownership explicit in the parameter atlas.



## Themes style non-data appearance {#themes}

`theme_ggchord`, minimal, publication and dark are alternative non-data style starting points; plasmid is tailored to single circular maps. Only the theme differs below; sequence data and mappings remain. Axis gaps, tick lengths and text offsets use documented physical units. Theme fields are not accver-keyed geometry lists. Hiding text does not remove sequences.

{{code:comparethemes}}
{{figure:comparethememinimal}}
{{figure:comparethemedark}}


## Exercises and selection rules {#adapt}

Switch pident from gradient to explicit steps, restore the gradient, then change only the alpha output range. Check that guides still represent the same data units. Give additional layers independent data/mapping and retain explicit user scales. See the [parameter atlas](parameter-atlas.html) for complete defaults, types and forwarding, and [link statistics](ribbon-statistics.html) for aggregation.
