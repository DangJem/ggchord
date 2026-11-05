Choose a lesson by task. A means chord, B collinear, C single-sequence circular. Check [compatibility boundaries](compatibility.html) for conditional entries. Complete the basics, then select advanced topics as needed. Browser find and site search index these tasks.

## Tasks and parameters {#tasks}

| Task | Function | Key parameters | Mode | Lesson and level |
| --- | --- | --- | --- | --- |
| Change backbone style | [`geom_seq()`](reference/geom_seq.html) | `seq_style, seq_backbone_gap, seq_backbone_width` | A / B / C | [Basics](plasmid-tuning.html#style) |
| Set radii and curvature | [`geom_seq()`](reference/geom_seq.html) | `seq_radius, seq_curvature` | A | [Advanced](chord-geometry.html) |
| Order and orient sequences | [`geom_seq()`](reference/geom_seq.html) | `seq_order, seq_orientation` | A / B | [Advanced](chord-geometry.html) |
| Route links around features | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `link_ctrl_point, link_gap, link_avoid` | A / B | [Advanced](chord-routing.html) |
| Connect points with arrows | [`geom_link_line()`](reference/geom_link_line.html) | `link_type, qpos/spos, arrow` | A / B | [Advanced](chord-branches.html#point-links) |
| Change gene arrow shapes | [`geom_gene()`](reference/geom_gene.html) | `arrow_head_style, short_feature` | A / B / C (bp) | [Advanced](feature-geometry.html#arrow-heads) |
| Highlight an interval | [`geom_seq_region()`](reference/geom_seq_region.html) | `region_side, region_width, region_offset` | A / B / C | [Advanced](chord-routing.html) |
| Colour by identity | [`scale_ribbon_fill_stepsn()`](reference/scale_ribbon_fill_stepsn.html) | `ribbon_fill` | A / B | [Basics](adjustments.html) |
| Use categories and legends | [`scale_seq_colour_manual()`](reference/scale_seq_colour_manual.html) | `values, guide` | A / B / C | [Advanced](styling.html) |
| Filter, deduplicate and merge HSPs | [`filter_ggchord_ribbons()`](reference/filter_ggchord_ribbons.html) | `min_pident, min_length` | A / B | [Advanced](ribbons.html) |
| Summarize link density | [`stat_ribbon_density()`](reference/stat_ribbon_density.html) | `bins, weight` | A / B | [Advanced](ribbons.html) |
| Bundle existing links | [`stat_ribbon_bundle()`](reference/stat_ribbon_bundle.html) | `min_bundle, group_by` | A / B | [Advanced](ribbons.html) |
| Align sequences left, center or right | [`coord_collinear()`](reference/coord_collinear.html) | `align` | B | [Basics](collinear-layout.html#alignment) |
| Anchor to common loci | [`coord_collinear()`](reference/coord_collinear.html) | `anchor` | B | [Advanced](collinear-layout.html#anchor) |
| Choose display directions automatically | [`coord_collinear()`](reference/coord_collinear.html) | `orientation, min_support` | B | [Advanced](collinear-layout.html#orientation) |
| Connect nonadjacent rows | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `adjacent_only, pairs, nonadjacent_linewidth` | B | [Advanced](collinear-routing.html) |
| Show local windows and breaks | [`segment_ggchord_data()`](reference/segment_ggchord_data.html) | `loci; coord_collinear: segments, segment_gap` | B | [Advanced](views.html#segments) |
| Draw a fixed-length scale bar | [`geom_seq_scale()`](reference/geom_seq_scale.html) | `length, location, offset` | B | [Advanced](collinear-tuning.html#scale) |
| Rotate a circular view | [`coord_circular()`](reference/coord_circular.html) | `rotation, direction, gap, fit` | C | [Basics](plasmid-tuning.html) |
| Stack feature annotations | [`position_feature_stack()`](reference/position_feature_stack.html) | `spacing, side, base_position` | A / B / C | [Advanced](annotations.html) |
| Arrange feature labels automatically | [`geom_feature_label_repel()`](reference/geom_feature_label_repel.html) | `label_layout, label_fit, label_min_size` | A / B / C (conditional) | [Advanced](labels.html) |
| Find primer bindings | [`find_primer_bindings()`](reference/find_primer_bindings.html) | `set, circular` | DNA | [Advanced](recognition.html) |
| Draw primers and labels | [`geom_primer()`](reference/geom_primer.html) | `position` | C (bp) | [Advanced](plasmids.html) |
| Search and filter restriction sites | [`find_restriction_sites()`](reference/find_restriction_sites.html) | `methylation, enzymes` | DNA | [Advanced](recognition.html#restriction) |
| Adjust restriction label leaders | [`geom_restriction_site()`](reference/geom_restriction_site.html) | `leader, label_style, label_side` | C (bp) | [Advanced](plasmids.html) |
| Draw GC bars and lines | [`geom_track_bar()`](reference/geom_track_tile.html) | `track_value, baseline` | A / B / C | [Advanced](tracks.html) |
| Draw a tile colour band | [`geom_track_tile()`](reference/geom_track_tile.html) | `track_fill, position` | A / B / C | [Advanced](tracks.html#tile) |
| Set track width and order | [`position_track()`](reference/position_track.html) | `track, side, order, width, gap` | A / B / C | [Advanced](tracks.html#sharing) |
| Import FASTA, GFF3 and alignments | [`read_fasta_lengths()`](reference/read_fasta_lengths.html) | `sequence_unit` | general | [Basics](formats.html) |
| Plot native protein comparisons | [`ggchord()`](reference/ggchord.html) | `seq_data: sequence_unit` | A / B (aa) | [Advanced](protein-alignments.html) |
| Preview and export publication figures | [`view_ggchord()`](reference/view_ggchord.html) | `width, height, units` | general | [Basics](output.html) |
| Diagnose and export layout | [`export_ggchord_layout()`](reference/export_ggchord_layout.html) | `plot` | general | [Advanced](faq.html) |
| Share a trunk at equal endpoints | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `link_branch, link_branch_fraction` | A / B | [Advanced](chord-branches.html#branches) |

## Three independent learning paths {#paths}

Chord: backbones → links and genes → geometry → routing → labels → four-genome case. Collinear: backbones → input and links → spacing → alignment and orientation → cross-row links and windows → four-genome case. Circular: backbone → features and labels → style → annotation layout → primers and restriction sites → tracks → complete pUC19 → optional reference gallery.

[Mappings, scales and guides](styling.html) is the shared home for appearance concepts; [output](output.html) covers preview, dimensions and saving. Reference pages provide full signatures; task links locate teaching and do not certify equal depth for every parameter.
