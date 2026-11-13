## Choose a plotting module {#modules}

Choose the figure you need; the four modules are independent learning routes. Confirm installation and complete a first plot, then add data, styling and output skills as needed.

{{taskcards}}

## Shared foundations {#foundations}

- [Download and install](installation.html)
- [Table fields and coordinate conventions](data.html)
- [Plotting concepts](modes.html)
- [Mappings, scales, guides and themes](styling.html)
- [Preview and export at publication size](output.html)

<!-- generated-tasks -->
Find a lesson by goal; check [compatibility](compatibility.html). Historical A/B/C labels below mean Chord/Collinear/Circular and do not describe the complete module structure.

## Tasks and parameters {#tasks}

| Task | Function | Key parameters | Mode | Lesson and level |
| --- | --- | --- | --- | --- |
| Change backbone style | [`geom_seq()`](reference/geom_seq.html) | `seq_style, seq_backbone_gap, seq_backbone_width` | A / B / C | [Basics](circular-layout.html#style) |
| Set radii and curvature | [`geom_seq()`](reference/geom_seq.html) | `seq_radius, seq_curvature` | A | [Advanced](chord-geometry.html) |
| Order and orient sequences | [`geom_seq()`](reference/geom_seq.html) | `seq_order, seq_orientation` | A / B | [Advanced](chord-geometry.html) |
| Route links around features | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `link_ctrl_point, link_gap, link_avoid` | A / B | [Advanced](multi-sequence-links.html#control) |
| Connect points with arrows | [`geom_link_line()`](reference/geom_link_line.html) | `link_type, qpos/spos, arrow` | A / B | [Advanced](multi-sequence-links.html#point-links) |
| Change gene arrow shapes | [`geom_gene()`](reference/geom_gene.html) | `arrow_head_style, short_feature` | A / B / C (bp) | [Advanced](feature-geometry.html#arrow-heads) |
| Highlight an interval | [`geom_seq_region()`](reference/geom_seq_region.html) | `region_side, region_width, region_offset` | A / B / C | [Advanced](multi-sequence-links.html#control) |
| Colour by identity | [`scale_ribbon_fill_stepsn()`](reference/scale_ribbon_fill_stepsn.html) | `ribbon_fill` | A / B | [Basics](styling.html#roles) |
| Use categories and legends | [`scale_seq_colour_manual()`](reference/scale_seq_colour_manual.html) | `values, guide` | A / B / C | [Advanced](styling.html) |
| Filter, deduplicate and merge HSPs | [`filter_ggchord_ribbons()`](reference/filter_ggchord_ribbons.html) | `min_pident, min_length` | A / B | [Advanced](ribbons.html) |
| Summarize link density | [`stat_ribbon_density()`](reference/stat_ribbon_density.html) | `bins, weight` | A / B | [Advanced](ribbons.html#statistics-density) |
| Bundle existing links | [`stat_ribbon_bundle()`](reference/stat_ribbon_bundle.html) | `min_bundle, group_by` | A / B | [Advanced](ribbons.html#statistics-bundle) |
| Align sequences left, center or right | [`coord_collinear()`](reference/coord_collinear.html) | `align` | B | [Basics](collinear-layout.html#alignment) |
| Anchor to common loci | [`coord_collinear()`](reference/coord_collinear.html) | `anchor` | B | [Advanced](collinear-layout.html#anchor) |
| Choose display directions automatically | [`coord_collinear()`](reference/coord_collinear.html) | `orientation, min_support` | B | [Advanced](collinear-layout.html#orientation) |
| Connect nonadjacent rows | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `adjacent_only, pairs, nonadjacent_linewidth` | B | [Advanced](multi-sequence-links.html#routes) |
| Show local windows and breaks | [`segment_ggchord_data()`](reference/segment_ggchord_data.html) | `windows; coord_collinear: segments, segment_gap` | B | [Advanced](views.html#segments) |
| Draw a fixed-length scale bar | [`geom_seq_scale()`](reference/geom_seq_scale.html) | `length, location, offset` | B | [Advanced](collinear-layout.html#gap) |
| Rotate a circular view | [`coord_circular()`](reference/coord_circular.html) | `rotation, direction, gap, fit` | C | [Basics](circular-layout.html#style) |
| Stack feature annotations | [`position_feature_stack()`](reference/position_feature_stack.html) | `spacing, side, base_position` | A / B / C | [Advanced](annotations.html) |
| Arrange feature labels automatically | [`geom_feature_label_repel()`](reference/geom_feature_label_repel.html) | `label_layout, label_fit, label_min_size` | A / B / C (conditional) | [Advanced](labels.html) |
| Find primer bindings | [`find_primer_bindings()`](reference/find_primer_bindings.html) | `set, circular` | DNA | [Advanced](primers.html#input) |
| Draw primers and labels | [`geom_primer()`](reference/geom_primer.html) | `position` | C (bp) | [Advanced](primers.html#primer-arcs) |
| Search and filter restriction sites | [`find_restriction_sites()`](reference/find_restriction_sites.html) | `methylation, enzymes` | DNA | [Advanced](recognition.html#restriction) |
| Adjust restriction label leaders | [`geom_restriction_site()`](reference/geom_restriction_site.html) | `leader, label_style, label_side` | C (bp) | [Advanced](restriction-sites.html#restriction-labels) |
| Draw GC bars and lines | [`geom_track_col()`](reference/geom_track.html) | `track_value, baseline` | A / B / C | [Advanced](tracks.html) |
| Draw a tile colour band | [`geom_track_tile()`](reference/geom_track_tile.html) | `track_fill, position` | A / B / C | [Advanced](tracks.html) |
| Set track width and order | [`position_track()`](reference/position_track.html) | `track, side, order, width, gap` | A / B / C | [Advanced](tracks.html) |
| Import FASTA, GFF3 and alignments | [`read_fasta_lengths()`](reference/read_fasta_lengths.html) | `sequence_unit` | general | [Basics](formats.html) |
| Plot native protein comparisons | [`ggchord()`](reference/ggchord.html) | `seq_data: sequence_unit` | A / B (aa) | [Advanced](protein-alignments.html) |
| Preview and export publication figures | [`view_ggchord()`](reference/view_ggchord.html) | `width, height, units` | general | [Basics](output.html) |
| Diagnose and export layout | [`export_ggchord_layout()`](reference/export_ggchord_layout.html) | `plot` | general | [Advanced](faq.html) |
| Share a trunk at equal endpoints | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `link_branch, link_branch_fraction` | A / B | [Advanced](multi-sequence-links.html#branches) |
| Style centre names and lengths | [`geom_seq_center_label()`](reference/geom_seq_center_label.html) | `name, length, separator, name_style, length_style` | C | [Advanced](circular-layout.html#center) |
| Set fixed label direction and side | [`geom_feature_label()`](reference/geom_feature_label.html) | `label_orientation, label_side, label_rotation` | A / B / C | [Advanced](circular-labels.html#fixed-labels) |
| Fit or abbreviate automatic labels | [`geom_feature_label_repel()`](reference/geom_feature_label_repel.html) | `label_fit, label_max_lines, external` | A / B / C | [Advanced](circular-labels.html#fitting) |
| Choose continuous and binned role scales | [`scale_ribbon_fill_gradientn()`](reference/scale_ribbon_fill_stepsn.html) | `colours, values, limits, guide` | A / B | [Advanced](scale-families.html) |
| Inspect windows and built-in catalogs | [`diagnose_ggchord()`](reference/diagnose_ggchord.html) | `plot, width, height; catalog: aliases` | general | [Advanced](workflows.html#tools) |
| Locate every API lesson | [`ggchord()`](reference/ggchord.html) | `all exported APIs` | general | [Advanced](api-curriculum.html) |
| Find parameter types and full defaults | [`ggchord()`](reference/ggchord.html) | `all public formals` | general | [Advanced](parameter-atlas.html) |
| Count source-position observations | [`geom_track()`](reference/geom_track.html) | `geom, stat, track, position` | A / B / C / Sankey | [Advanced](tracks-statistics.html#tracks) |
| Route repeated sequence instances | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `between, allow_repeated` | A / B | [Advanced](instances.html#instances) |
| Draw mirrored interval relationships | [`coord_sankey()`](reference/coord_sankey.html) | `column_gap, self_links` | Sankey | [Advanced](sankey.html#sankey) |
| Show HSPs on independent axes | [`coord_dotplot()`](reference/coord_dotplot.html) | `x_seqs, y_seqs, x_orientation, y_orientation` | Dotplot | [Advanced](dotplot-start.html#dotplot) |
| Retain four alignment endpoints | [`geom_alignment_segment()`](reference/geom_alignment_segment.html) | `mapping, data, allow_repeated` | Dotplot | [Advanced](dotplot-start.html#dotplot) |
| Draw cyclic radial synteny | [`coord_radial_synteny()`](reference/coord_radial_synteny.html) | `start_angle, clockwise, inner_radius` | Radial | [Advanced](radial-synteny.html#radial) |
| Compose source-anchored views | [`ggchord_fusion()`](reference/ggchord_fusion.html) | `center, center_size, layout` | Fusion | [Advanced](fusion-routing.html#intervals) |
| Connect real intervals across panels | [`fusion_ribbon()`](reference/fusion_ribbon.html) | `data, qpanel, spanel, gradient, tip` | Fusion | [Advanced](fusion-routing.html#intervals) |
| Audit composition sources and physical anchors | [`export_ggchord_fusion()`](reference/export_ggchord_fusion.html) | `width, height` | Fusion | [Advanced](fusion-routing.html#export) |
| Import paired BED intervals | [`read_bedpe()`](reference/read_bedpe.html) | `file, files, source_file, unmapped` | Input | [Advanced](import-intervals.html#bedpe) |
| Show actual characters and gaps from an alignment | [`geom_residue()`](reference/geom_residue.html) | `mapping, overflow, size` | Cartesian aligned excerpts | [Advanced](aligned-sequences.html#aligned-characters) |
| Look up source positions and alignment columns | [`map_alignment_position()`](reference/map_alignment_position.html) | `data, accver, position, inverse` | DNA / RNA / AA | [Advanced](aligned-sequences.html#aligned-source) |
| Connect a genomic region to its independently rebuilt detail | [`magnify_ggchord()`](reference/magnify_ggchord.html) | `region, detail_plot, side, connector` | Circular / Chord + Collinear | [Advanced](fusion-details.html#magnify) |

## Independent learning paths {#paths}

Chord: backbones → links and genes → geometry → routing → labels → four-genome case. Collinear: backbones → input and links → spacing → alignment and orientation → cross-row links and windows → four-genome case. Circular: backbone → features and labels → style → annotation layout → primers and restriction sites → tracks → complete pUC19 → optional reference gallery.

[Mappings, scales and guides](styling.html) is the shared home for appearance concepts; [output](output.html) covers preview, dimensions and saving. Reference pages provide full signatures; task links locate teaching and do not certify equal depth for every parameter.
