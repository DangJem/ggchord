Find purpose, primary section, modes, complete cases and full parameters by API. References define exact defaults, formats, units and constraints; maintenance levels appear at the end. Upstream reexports and spelling aliases are explained by family rather than duplicate sections.

<a id="apis"></a>

## API usage index {#apis}

| API/reference | Purpose | Usage section | Mode/input boundary | Parameter atlas | Case/workflow |
| --- | --- | --- | --- | --- | --- |
| [`adapt_ggchord`](../../man/adapt_ggchord.Rd) | Opt into bounded theme scaling and joint guide layout at the target size | [Usage](styling.md#adaptive) | fixed/scale/adaptive; explicit settings remain authoritative; no global collision guarantee | [↗](parameter-atlas.md#api-adapt_ggchord) | [↗](styling.md) |
| [`aes`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-aes) | [↗](chord.md) |
| [`after_scale`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-after_scale) | [↗](chord.md) |
| [`after_stat`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-after_stat) | [↗](chord.md) |
| [`annotate`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-annotate) | [↗](chord.md) |
| [`arrow`](../../man/grid-helpers.Rd) | Physical units and direction arrows | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-arrow) | [↗](chord.md) |
| [`as_ggchord_alignment`](../../man/as_ggchord_alignment.Rd) | Adapt existing aligned strings with exact nullable source positions | [Usage](aligned-sequences.md#aligned-source) | Existing DNA/RNA/AA columns; Cartesian local cells; explicit source mapping | [↗](parameter-atlas.md#api-as_ggchord_alignment) | [↗](formats.md) |
| [`as_ggchord_feature`](../../man/as_ggchord_feature.Rd) | Convert intervals to generic ggchord features | [Usage](import-sequences.md#adapters) | table/S4 adapter; coordinates and biological direction retained | [↗](parameter-atlas.md#api-as_ggchord_feature) | [↗](formats.md) |
| [`as_ggchord_link`](../../man/as_ggchord_link.Rd) | Convert paired genomic intervals to ggchord links | [Usage](import-sequences.md#adapters) | table/GInteractions adapter; no invented alignment statistics | [↗](parameter-atlas.md#api-as_ggchord_link) | [↗](formats.md) |
| [`as_ggchord_seq`](../../man/as_ggchord_seq.Rd) | Convert sequence information to ggchord input | [Usage](import-sequences.md#adapters) | table/S4 adapter; complete lengths and explicit IDs/units | [↗](parameter-atlas.md#api-as_ggchord_seq) | [↗](formats.md) |
| [`build_ggchord_cds`](../../man/build_ggchord_cds.Rd) | Build an explicit spliced CDS coordinate map | [Usage](projection.md#boundaries) | known genomic CDS, phase and genetic code; no guessed correspondence | [↗](parameter-atlas.md#api-build_ggchord_cds) | [↗](homology.md) |
| [`bundle_ggchord_ribbons`](../../man/bundle_ggchord_ribbons.Rd) | Bundle dense alignment ribbons explicitly | [Usage](ribbons.md#statistics-bundle) | explicit normalized-bin aggregation; original observations retained separately | [↗](parameter-atlas.md#api-bundle_ggchord_ribbons) | [↗](chord.md) |
| [`clean_ggchord_data`](../../man/clean_ggchord_data.Rd) | Clean ggchord input data with explicit, report-driven policies | [Usage](data.md#repairs) | explicit repair policies with original tables retained | [↗](parameter-atlas.md#api-clean_ggchord_data) | [↗](real-data.md) |
| [`coord_chord`](../../man/coord_chord.Rd) | Chord diagram coordinate system | [Usage](chord-geometry.md#shape) | chord | [↗](parameter-atlas.md#api-coord_chord) | [↗](chord.md) |
| [`coord_circular`](../../man/coord_circular.Rd) | Single-sequence circular coordinate system | [Usage](circular-layout.md#coordinates) | one sequence; bp or aa; no ribbon/link layers | [↗](parameter-atlas.md#api-coord_circular) | [↗](circular-start.md) |
| [`coord_collinear`](../../man/coord_collinear.Rd) | Linear multi-sequence collinearity coordinates | [Usage](collinear-layout.md#parameters) | collinear; bp or aa; circular-only local controls rejected | [↗](parameter-atlas.md#api-coord_collinear) | [↗](collinear.md) |
| [`coord_dotplot`](../../man/coord_dotplot.Rd) | Two-axis HSP projection | [Usage](dotplot-start.md#dotplot) | Identity HSP segments on independent native axes | [↗](parameter-atlas.md#api-coord_dotplot) | [↗](dotplot-start.md) |
| [`coord_radial_synteny`](../../man/coord_radial_synteny.Rd) | Project linear sequences onto radial spokes | [Usage](radial-synteny.md#radial) | Proportional/equal spokes, cyclic actual HSPs and shared feature/Track frames | [↗](parameter-atlas.md#api-coord_radial_synteny) | [↗](radial-synteny.md) |
| [`coord_sankey`](../../man/coord_sankey.Rd) | Mirrored query-left/subject-right interval relationships | [Usage](sankey.md#sankey) | Mirrored two-column relationships | [↗](parameter-atlas.md#api-coord_sankey) | [↗](multi-sequence.md) |
| [`deduplicate_ggchord_ribbons`](../../man/deduplicate_ggchord_ribbons.Rd) | Deduplicate alignment ribbons | [Usage](ribbons.md#transformations) | select existing representatives with source membership | [↗](parameter-atlas.md#api-deduplicate_ggchord_ribbons) | [↗](chord.md) |
| [`diagnose_ggchord`](../../man/diagnose_ggchord.Rd) | Inspect one target-size render without changing the plot | [Usage](diagnostics.md#target-size) | one target-size device; advisory, not universal quality certification | [↗](parameter-atlas.md#api-diagnose_ggchord) | [↗](chord.md) |
| [`element_blank`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-element_blank) | [↗](chord.md) |
| [`element_line`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-element_line) | [↗](chord.md) |
| [`element_rect`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-element_rect) | [↗](chord.md) |
| [`element_text`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-element_text) | [↗](chord.md) |
| [`expansion`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-expansion) | [↗](chord.md) |
| [`export_ggchord_fusion`](../../man/export_ggchord_fusion.Rd) | Export measured composition fronts and sources | [Usage](fusion-routing.md#export) | Registered genomic plots, nested Fusion and registered character details; explicit panel/instance routing; Dotplot excluded | [↗](parameter-atlas.md#api-export_ggchord_fusion) | [↗](fusion.md) |
| [`export_ggchord_layout`](../../man/export_ggchord_layout.Rd) | Export computed ggchord geometry | [Usage](diagnostics.md#geometry) | layer-aware public tables for all supported plots | [↗](parameter-atlas.md#api-export_ggchord_layout) | [↗](chord.md) |
| [`filter_ggchord_ribbons`](../../man/filter_ggchord_ribbons.Rd) | Filter alignment ribbons before plotting | [Usage](ribbons.md#input) | observed alignment filtering, not recomputation | [↗](parameter-atlas.md#api-filter_ggchord_ribbons) | [↗](chord.md) |
| [`filter_restriction_sites`](../../man/filter_restriction_sites.Rd) | Filter restriction sites for display | [Usage](restriction-sites.md#parameters) | display set selection preserves complete-target counts | [↗](parameter-atlas.md#api-filter_restriction_sites) | [↗](plasmids.md) |
| [`find_common_features`](../../man/find_common_features.Rd) | Find common biological features in DNA sequences | [Usage](recognition.md#features) | DNA/reference feature recognition with method/provenance, not wet-lab validation | [↗](parameter-atlas.md#api-find_common_features) | [↗](plasmids.md) |
| [`find_ggchord_loci`](../../man/find_ggchord_loci.Rd) | Find genomic loci around explicitly identified features | [Usage](feature-geometry.md#locus-selection) | genomic workflow with explicit target IDs and source strands | [↗](parameter-atlas.md#api-find_ggchord_loci) | [↗](views.md) |
| [`find_primer_bindings`](../../man/find_primer_bindings.Rd) | Find sequencing-primer binding sites | [Usage](primers.md#input) | actual DNA and explicit oligonucleotides/catalogue; matches not experiments | [↗](parameter-atlas.md#api-find_primer_bindings) | [↗](plasmids.md) |
| [`find_restriction_sites`](../../man/find_restriction_sites.Rd) | Find restriction pattern matches and cleavage coordinates | [Usage](restriction-sites.md#input) | actual DNA and explicit methylation assumptions; full site results | [↗](parameter-atlas.md#api-find_restriction_sites) | [↗](plasmids.md) |
| [`focus_ggchord_data`](../../man/focus_ggchord_data.Rd) | Focus ggchord data on one or more loci | [Usage](views.md#focus) | explicit source windows with reversible origins and clipping reports | [↗](parameter-atlas.md#api-focus_ggchord_data) | [↗](homology.md) |
| [`focus_ggchord_track`](../../man/focus_ggchord_track.Rd) | Focus quantitative intervals or samples without bridging omitted windows | [Usage](tracks-statistics.md#track-focus) | point/interval track data; shares declared loci with focus/segment | [↗](parameter-atlas.md#api-focus_ggchord_track) | [↗](tracks.md) |
| [`fusion_panel`](../../man/fusion_panel.Rd) | Register a panel at a physical side | [Usage](fusion-layout.md#advanced-placement) | Registered genomic/character panels; eight sides, explicit boxes and nested paths | [↗](parameter-atlas.md#api-fusion_panel) | [↗](fusion.md) |
| [`fusion_ribbon`](../../man/fusion_ribbon.Rd) | Connect full source intervals | [Usage](fusion-routing.md#intervals) | Registered genomic plots, nested Fusion and registered character details; explicit panel/instance routing; Dotplot excluded | [↗](parameter-atlas.md#api-fusion_ribbon) | [↗](fusion.md) |
| [`geom_alignment_segment`](../../man/geom_alignment_segment.Rd) | Draw complete HSP endpoints | [Usage](dotplot-start.md#dotplot) | Identity HSP segments on independent native axes | [↗](parameter-atlas.md#api-geom_alignment_segment) | [↗](dotplot-start.md) |
| [`geom_feature`](../../man/geom_feature.Rd) | Draw generic genomic features | [Usage](circular-features.md#feature-shapes) | all three; generic annotations including unknown strand and aa; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_feature) | [↗](chord.md) |
| [`geom_feature_label`](../../man/geom_feature_label.Rd) | Label generic genomic features | [Usage](circular-labels.md#fixed-labels) | all three; linear uses vertical/horizontal rows for radial/tangent; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_feature_label) | [↗](chord.md) |
| [`geom_gene_label`](../../man/geom_gene_label.Rd) | Add a gene label layer | [Usage](labels.md#fixed) | bp gene labels; circular/linear orientation differs; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_gene_label) | [↗](chord.md) |
| [`geom_gene_label_repel`](../../man/geom_gene_label_repel.Rd) | Add an automatically arranged gene label layer | [Usage](labels.md#comparison) | bp gene labels; linear shares row layout; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_gene_label_repel) | [↗](chord.md) |
| [`geom_link_line`](../../man/geom_link_line.Rd) | Draw point-to-point sequence links | [Usage](multi-sequence-links.md#point-links) | chord/collinear; explicit qpos/spos; circular rejected; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_link_line) | [↗](chord.md) |
| [`geom_link_ribbon`](../../man/geom_link_ribbon.Rd) | Add an alignment ribbon layer | [Usage](multi-sequence-links.md#branches) | chord/collinear; circular rejected; linear adjacency follows displayed rows; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_link_ribbon) | [↗](chord.md) |
| [`geom_primer`](../../man/geom_primer.Rd) | Draw primer binding ranges on a circular sequence | [Usage](primers.md#primer-arcs) | bp single circular recommended; single chord composition available; linear rejected | [↗](parameter-atlas.md#api-geom_primer) | [↗](plasmids.md) |
| [`geom_primer_label_repel`](../../man/geom_primer_label_repel.Rd) | Arrange primer labels and leaders | [Usage](primers.md#primer-labels) | bp single circular recommended; linear rejected | [↗](parameter-atlas.md#api-geom_primer_label_repel) | [↗](plasmids.md) |
| [`geom_residue`](../../man/geom_residue.Rd) | Draw local characters with measured physical space and explicit tile/error policy | [Usage](aligned-sequences.md#aligned-characters) | Existing DNA/RNA/AA columns; Cartesian local cells; explicit source mapping | [↗](parameter-atlas.md#api-geom_residue) | [↗](formats.md) |
| [`geom_restriction_site`](../../man/geom_restriction_site.Rd) | Draw restriction sites, labels and deterministic leaders | [Usage](restriction-sites.md#restriction-labels) | bp single circular recommended; linear rejected | [↗](parameter-atlas.md#api-geom_restriction_site) | [↗](plasmids.md) |
| [`geom_seq`](../../man/geom_seq.Rd) | Add a sequence arc layer | [Usage](circular-layout.md#backbones) | all three; linear rejects radius/curvature/ring; circular owns gap/direction; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_seq) | [↗](chord.md) |
| [`geom_seq_break`](../../man/geom_seq_break.Rd) | Mark cropped ends of linear sequences | [Usage](views.md#segments) | collinear only, focused/segmented truncated-end metadata required; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_seq_break) | [↗](collinear.md) |
| [`geom_seq_center_label`](../../man/geom_seq_center_label.Rd) | Add a circular-sequence centre label | [Usage](circular-layout.md#center) | one sequence; circular recommended; linear rejected | [↗](parameter-atlas.md#api-geom_seq_center_label) | [↗](chord.md) |
| [`geom_seq_label`](../../man/geom_seq_label.Rd) | Add a sequence label layer | [Usage](multi-sequence-labels.md#sequence-parameters) | all three; linear rejects explicit seq_label_radius; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_seq_label) | [↗](chord.md) |
| [`geom_seq_region`](../../man/geom_seq_region.Rd) | Highlight regions along sequence arcs | [Usage](annotations.md#region-parameters) | all three, valid source intervals and units; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_seq_region) | [↗](chord.md) |
| [`geom_seq_scale`](../../man/geom_seq_scale.Rd) | Add a shared sequence-position scale bar to a linear display | [Usage](collinear-layout.md#gap) | collinear only; native bp/aa ruler | [↗](parameter-atlas.md#api-geom_seq_scale) | [↗](collinear.md) |
| [`geom_track`](../../man/geom_track.Rd) | Source-space statistical tracks | [Usage](tracks-statistics.md#tracks) | Chord/Circular/Collinear/Sankey within the first support matrix; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_track) | [↗](multi-sequence.md) |
| [`geom_track_bar`](../../man/geom_track.Rd) | Count observations at source positions | [Usage](tracks-statistics.md#tracks) | Count observations at source positions in compatible backbone views; zero baseline; binning is a separate explicit Stat | [↗](parameter-atlas.md#api-geom_track_bar) | [↗](multi-sequence.md) |
| [`geom_track_col`](../../man/geom_track.Rd) | Display supplied numerical heights | [Usage](tracks-statistics.md#tracks) | Chord/Circular/Collinear/Sankey; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_track_col) | [↗](multi-sequence.md) |
| [`geom_track_line`](../../man/geom_track_tile.Rd) | Draw explicit quantitative tracks along sequence backbones | [Usage](tracks.md#multiple) | all three; increasing positions per group; closed requires complete circular gap=0; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_track_line) | [↗](tracks.md) |
| [`geom_track_tile`](../../man/geom_track_tile.Rd) | Draw explicit quantitative tracks along sequence backbones | [Usage](tracks.md#tile) | all three; default or named track; overlapping input tiles keep drawing order; Radial per new-view support matrix | [↗](parameter-atlas.md#api-geom_track_tile) | [↗](tracks.md) |
| [`get_chord_layout`](../../man/get_chord_layout.Rd) | Get the chord layout from a ggchord plot | [Usage](diagnostics.md#geometry) | complete documented layout for all supported plots | [↗](parameter-atlas.md#api-get_chord_layout) | [↗](chord.md) |
| [`ggchord`](../../man/ggchord.Rd) | Create sequence plots from tables and independent layers | [Usage](modes.md#composition) | constructor; one declared unit per plot | [↗](parameter-atlas.md#api-ggchord) | [↗](real-data.md) |
| [`ggchord_fusion`](../../man/ggchord_fusion.Rd) | Compose independent genomic panels | [Usage](fusion-routing.md#intervals) | Registered genomic/character panels; eight sides, explicit boxes and nested paths | [↗](parameter-atlas.md#api-ggchord_fusion) | [↗](fusion.md) |
| [`ggsave`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-ggsave) | [↗](chord.md) |
| [`ggtitle`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-ggtitle) | [↗](chord.md) |
| [`guide_ggchord_colorbar`](../../man/guide_ggchord_legend.Rd) | Guides for ggchord role aesthetics; spelling alias of guide_ggchord_colourbar | [Usage](scale-families.md#guides) | matching discrete/continuous role guide; all supported plots | [↗](parameter-atlas.md#api-guide_ggchord_colorbar) | [↗](chord.md) |
| [`guide_ggchord_colourbar`](../../man/guide_ggchord_legend.Rd) | Guides for ggchord role aesthetics | [Usage](scale-families.md#guides) | matching discrete/continuous role guide; all supported plots | [↗](parameter-atlas.md#api-guide_ggchord_colourbar) | [↗](chord.md) |
| [`guide_ggchord_legend`](../../man/guide_ggchord_legend.Rd) | Guides for ggchord role aesthetics | [Usage](scale-families.md#guides) | matching discrete/continuous role guide; all supported plots | [↗](parameter-atlas.md#api-guide_ggchord_legend) | [↗](chord.md) |
| [`guide_seq_axis`](../../man/guide_seq_axis.Rd) | Select sequence axes and control tick label overlap | [Usage](scale-families.md#guides) | Backbone views with sequence axes; select focused axes by window ID | [↗](parameter-atlas.md#api-guide_seq_axis) | [↗](tracks.md) |
| [`guide_track_axis`](../../man/guide_track_axis.Rd) | Set tick angle, overlap policy and local theme on a quantitative track axis | [Usage](scale-families.md#guides) | Views using a local quantitative track value scale | [↗](parameter-atlas.md#api-guide_track_axis) | [↗](tracks.md) |
| [`guides`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-guides) | [↗](chord.md) |
| [`labs`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-labs) | [↗](chord.md) |
| [`last_plot`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-last_plot) | [↗](chord.md) |
| [`layout_sequences`](../../man/layout_sequences.Rd) | Set backbone order, direction, spacing and relative display lengths | [Usage](multi-sequence-layout.md#order) | Multi-sequence chord and linear layouts; instance IDs can override source directions; circular and dotplot axes own their directions | [↗](parameter-atlas.md#api-layout_sequences) | [↗](chord.md) |
| [`link_ggchord_clusters`](../../man/link_ggchord_clusters.Rd) | Convert explicit homology-cluster membership to feature links | [Usage](homology.md#links) | explicit biological membership; matching names alone are insufficient | [↗](parameter-atlas.md#api-link_ggchord_clusters) | [↗](homology.md) |
| [`magnify_ggchord`](../../man/magnify_ggchord.Rd) | Connect a real source interval to an explicit independent local detail plot | [Usage](fusion-details.md#magnify) | Registered genomic plots, nested Fusion and registered character details; explicit panel/instance routing; Dotplot excluded | [↗](parameter-atlas.md#api-magnify_ggchord) | [↗](fusion.md) |
| [`map_alignment_position`](../../man/map_alignment_position.Rd) | Look up exact alignment columns and source positions in a selected row | [Usage](aligned-sequences.md#aligned-characters) | Existing DNA/RNA/AA columns; Cartesian local cells; explicit source mapping | [↗](parameter-atlas.md#api-map_alignment_position) | [↗](formats.md) |
| [`map_ggchord_position`](../../man/map_ggchord_position.Rd) | Map source positions and display coordinates | [Usage](radial-synteny.md#export) | Native/source positions; explicit Dotplot axis, or both x/y for radial inverse | [↗](parameter-atlas.md#api-map_ggchord_position) | [↗](dotplot-start.md) |
| [`margin`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-margin) | [↗](chord.md) |
| [`merge_ggchord_ribbons`](../../man/merge_ggchord_ribbons.Rd) | Merge adjacent or overlapping alignment blocks of the same sequence pair | [Usage](ribbons.md#transformations) | compatible display envelopes with estimation/provenance flags | [↗](parameter-atlas.md#api-merge_ggchord_ribbons) | [↗](chord.md) |
| [`optimize_ggchord_layout`](../../man/optimize_ggchord_layout.Rd) | Optimize sequence order and orientation for ribbon readability | [Usage](diagnostics.md#target-size) | explicit input sequence order/orientation heuristic, not biological inference | [↗](parameter-atlas.md#api-optimize_ggchord_layout) | [↗](chord.md) |
| [`position_feature_stack`](../../man/position_feature_stack.Rd) | Stack overlapping genes or features on local lanes | [Usage](circular-features.md#stacking) | feature lanes; base Position before overlap allocation | [↗](parameter-atlas.md#api-position_feature_stack) | [↗](chord.md) |
| [`position_identity`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-position_identity) | [↗](chord.md) |
| [`position_offset`](../../man/position_offset.Rd) | Place both feature strands on one shared band | [Usage](circular-features.md#stacking) | shared local offset for compatible sequence coordinates | [↗](parameter-atlas.md#api-position_offset) | [↗](circular-features.md) |
| [`position_strand`](../../man/position_strand.Rd) | Place features on opposite sides according to strand | [Usage](parameters.md#strand) | local signed placement on valid feature/gene layers | [↗](parameter-atlas.md#api-position_strand) | [↗](parameters.md) |
| [`PositionIdentity`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-PositionIdentity) | [↗](chord.md) |
| [`prepare_ggchord_data`](../../man/prepare_ggchord_data.Rd) | Prepare named ggchord tables with an explicit audit trail | [Usage](data.md#validation) | workflow preparation; defaults differ from standalone cleaning | [↗](parameter-atlas.md#api-prepare_ggchord_data) | [↗](real-data.md) |
| [`primer_catalog`](../../man/primer_catalog.Rd) | Inspect the bundled sequencing-primer catalogue | [Usage](recognition.md#catalog) | catalogue metadata, not a binding search | [↗](parameter-atlas.md#api-primer_catalog) | [↗](plasmids.md) |
| [`project_ggchord_features`](../../man/project_ggchord_features.Rd) | Project protein domains or genomic features through a spliced CDS map | [Usage](projection.md#cds) | supported amino-acid intervals through explicit CDS segment maps | [↗](parameter-atlas.md#api-project_ggchord_features) | [↗](homology.md) |
| [`project_ggchord_links`](../../man/project_ggchord_links.Rd) | Project explicitly ungapped protein links to genomic coding segments | [Usage](projection.md#links) | known paired protein/CDS correspondence; gaps require justified handling | [↗](parameter-atlas.md#api-project_ggchord_links) | [↗](homology.md) |
| [`read_aligned_fasta`](../../man/read_aligned_fasta.Rd) | Read existing DNA/RNA/AA alignment characters without computing a new alignment | [Usage](aligned-sequences.md#aligned-characters) | Existing DNA/RNA/AA columns; Cartesian local cells; explicit source mapping | [↗](parameter-atlas.md#api-read_aligned_fasta) | [↗](formats.md) |
| [`read_bed`](../../man/read_bed.Rd) | Read BED intervals as ggchord features | [Usage](import-intervals.md#zero-based) | zero-based half-open BED converted once | [↗](parameter-atlas.md#api-read_bed) | [↗](formats.md) |
| [`read_bedgraph`](../../man/read_bedgraph.Rd) | Import continuous interval signals without inferring units | [Usage](import-signals.md#signals) | quantitative Track input; BigWig optional on supported platforms | [↗](parameter-atlas.md#api-read_bedgraph) | [↗](formats.md) |
| [`read_bedpe`](../../man/read_bedpe.Rd) | Read paired genomic intervals without inventing alignment evidence | [Usage](import-intervals.md#bedpe) | Canonical q/s intervals, paired-anchor direction and provenance | [↗](parameter-atlas.md#api-read_bedpe) | [↗](formats.md) |
| [`read_bigwig`](../../man/read_bigwig.Rd) | Import continuous interval signals without inferring units | [Usage](import-signals.md#signals) | quantitative Track input; BigWig optional on supported platforms | [↗](parameter-atlas.md#api-read_bigwig) | [↗](formats.md) |
| [`read_blast`](../../man/read_blast.Rd) | Read one or more BLAST tabular output files into ribbon_data format | [Usage](alignments.md#standard) | local observed BLAST records; retain q/s coordinate meaning | [↗](parameter-atlas.md#api-read_blast) | [↗](real-data.md) |
| [`read_chain`](../../man/read_chain.Rd) | Read existing genomic alignment blocks and source coordinates | [Usage](import-relations.md#genome-alignments) | optional Biopython; MAF cells and chain intervals | [↗](parameter-atlas.md#api-read_chain) | [↗](formats.md) |
| [`read_delta`](../../man/read_delta.Rd) | Read existing delta alignments through an optional native decoder | [Usage](import-relations.md#delta) | optional show-coords; bp coordinates and declared alignment metric unit | [↗](parameter-atlas.md#api-read_delta) | [↗](formats.md) |
| [`read_fasta_lengths`](../../man/read_fasta_lengths.Rd) | Read one or more FASTA files into seq_data format | [Usage](aligned-sequences.md#characters) | local FASTA; declared bp/aa unit and complete lengths | [↗](parameter-atlas.md#api-read_fasta_lengths) | [↗](real-data.md) |
| [`read_feature_relations`](../../man/read_feature_relations.Rd) | Read precomputed feature pairs and memberships | [Usage](import-relations.md#feature-relations) | Exact feature IDs; contiguous supplied anchors; no full lengths inferred | [↗](parameter-atlas.md#api-read_feature_relations) | [↗](formats.md) |
| [`read_genbank`](../../man/read_genbank.Rd) | Import local GenBank or GBFF records using an optional mature parser | [Usage](import-sequences.md#genbank) | local GenBank; explicit installed Biopython parser; provenance retained | [↗](parameter-atlas.md#api-read_genbank) | [↗](formats.md) |
| [`read_gff3`](../../man/read_gff3.Rd) | Read one or more GFF3 files into gene_data format | [Usage](real-data.md#read) | local one-based inclusive genomic annotation | [↗](parameter-atlas.md#api-read_gff3) | [↗](real-data.md) |
| [`read_maf`](../../man/read_maf.Rd) | Read existing genomic alignment blocks and source coordinates | [Usage](import-relations.md#genome-alignments) | optional Biopython; MAF cells and chain intervals | [↗](parameter-atlas.md#api-read_maf) | [↗](formats.md) |
| [`read_paf`](../../man/read_paf.Rd) | Read PAF mappings as ggchord links | [Usage](import-intervals.md#zero-based) | PAF with complete source lengths and strand-aware endpoint conversion | [↗](parameter-atlas.md#api-read_paf) | [↗](formats.md) |
| [`read_signature_matches`](../../man/read_signature_matches.Rd) | Read existing protein signature locations and source metadata | [Usage](import-relations.md#protein-signatures) | AA features; source-native scores; JSON needs optional jsonlite | [↗](parameter-atlas.md#api-read_signature_matches) | [↗](formats.md) |
| [`read_vcf`](../../man/read_vcf.Rd) | Import variant records and explicitly decoded structural breakpoints | [Usage](import-relations.md#variants) | VCF; BCF with explicit pysam Python; optional structural decoder | [↗](parameter-atlas.md#api-read_vcf) | [↗](formats.md) |
| [`register_ggchord_alignment`](../../man/register_ggchord_alignment.Rd) | Register exact character/source positions for independent Fusion details | [Usage](fusion-details.md#character-detail) | Explicit source/column map; independent Cartesian detail; complete observed codons | [↗](parameter-atlas.md#api-register_ggchord_alignment) | [↗](fusion.md) |
| [`rel`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-rel) | [↗](chord.md) |
| [`scale_feature_fill_manual`](../../man/scale_gene_fill_manual.Rd) | Gene, feature and region fill scales | [Usage](scale-families.md#features) | all three; matching generic feature role/preset | [↗](parameter-atlas.md#api-scale_feature_fill_manual) | [↗](chord.md) |
| [`scale_feature_fill_plasmid`](../../man/scale_feature_fill_plasmid.Rd) | Plasmid-map feature presets | [Usage](scale-families.md#features) | all three; matching generic feature role/preset | [↗](parameter-atlas.md#api-scale_feature_fill_plasmid) | [↗](chord.md) |
| [`scale_feature_shape_manual`](../../man/scale_feature_shape_manual.Rd) | Feature shape scale | [Usage](scale-families.md#features) | all three; matching generic feature role/preset | [↗](parameter-atlas.md#api-scale_feature_shape_manual) | [↗](chord.md) |
| [`scale_feature_shape_plasmid`](../../man/scale_feature_fill_plasmid.Rd) | Plasmid-map feature presets | [Usage](scale-families.md#features) | all three; matching generic feature role/preset | [↗](parameter-atlas.md#api-scale_feature_shape_plasmid) | [↗](chord.md) |
| [`scale_fill_residue`](../../man/scale_fill_residue.Rd) | Choose chemical-property colors for nucleotide or protein residues | [Usage](aligned-sequences.md#aligned-characters) | Existing DNA/RNA/AA columns; Cartesian local cells; explicit source mapping | [↗](parameter-atlas.md#api-scale_fill_residue) | [↗](formats.md) |
| [`scale_gene_fill_manual`](../../man/scale_gene_fill_manual.Rd) | Gene, feature and region fill scales | [Usage](scale-families.md#rings) | bp gene layers; compatible backbone coordinates | [↗](parameter-atlas.md#api-scale_gene_fill_manual) | [↗](chord.md) |
| [`scale_link_alpha_continuous`](../../man/scale_link.Rd) | Scales for point links | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_alpha_continuous) | [↗](chord.md) |
| [`scale_link_alpha_manual`](../../man/scale_link.Rd) | Scales for point links | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_alpha_manual) | [↗](chord.md) |
| [`scale_link_color_gradientn`](../../man/scale_link.Rd) | Scales for point links; spelling alias of scale_link_colour_gradientn | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_color_gradientn) | [↗](chord.md) |
| [`scale_link_color_manual`](../../man/scale_link.Rd) | Scales for point links; spelling alias of scale_link_colour_manual | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_color_manual) | [↗](chord.md) |
| [`scale_link_colour_gradientn`](../../man/scale_link.Rd) | Scales for point links | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_colour_gradientn) | [↗](chord.md) |
| [`scale_link_colour_manual`](../../man/scale_link.Rd) | Scales for point links | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_colour_manual) | [↗](chord.md) |
| [`scale_link_linetype_manual`](../../man/scale_link.Rd) | Scales for point links | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_linetype_manual) | [↗](chord.md) |
| [`scale_link_linewidth_continuous`](../../man/scale_link.Rd) | Scales for point links | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_linewidth_continuous) | [↗](chord.md) |
| [`scale_link_linewidth_manual`](../../man/scale_link.Rd) | Scales for point links | [Usage](scale-families.md#point-links) | chord/collinear matching point-link roles | [↗](parameter-atlas.md#api-scale_link_linewidth_manual) | [↗](chord.md) |
| [`scale_region_fill_manual`](../../man/scale_gene_fill_manual.Rd) | Gene, feature and region fill scales | [Usage](scale-families.md#rings) | all three; matching region_fill role | [↗](parameter-atlas.md#api-scale_region_fill_manual) | [↗](chord.md) |
| [`scale_ribbon_alpha_continuous`](../../man/scale_ribbon_alpha_continuous.Rd) | Ribbon alpha, outline and linetype scales | [Usage](scale-families.md#categories) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_alpha_continuous) | [↗](chord.md) |
| [`scale_ribbon_alpha_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | Ribbon alpha, outline and linetype scales | [Usage](scale-families.md#categories) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_alpha_manual) | [↗](chord.md) |
| [`scale_ribbon_color_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | Ribbon alpha, outline and linetype scales; spelling alias of scale_ribbon_colour_manual | [Usage](scale-families.md#categories) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_color_manual) | [↗](chord.md) |
| [`scale_ribbon_colour_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | Ribbon alpha, outline and linetype scales | [Usage](scale-families.md#categories) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_colour_manual) | [↗](chord.md) |
| [`scale_ribbon_fill_gradientn`](../../man/scale_ribbon_fill_stepsn.Rd) | Ribbon fill scales | [Usage](scale-families.md#numeric) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_fill_gradientn) | [↗](chord.md) |
| [`scale_ribbon_fill_identity`](../../man/scale_ribbon_fill_stepsn.Rd) | Ribbon fill scales | [Usage](scale-families.md#identity) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_fill_identity) | [↗](chord.md) |
| [`scale_ribbon_fill_manual`](../../man/scale_ribbon_fill_stepsn.Rd) | Ribbon fill scales | [Usage](scale-families.md#categories) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_fill_manual) | [↗](chord.md) |
| [`scale_ribbon_fill_stepsn`](../../man/scale_ribbon_fill_stepsn.Rd) | Ribbon fill scales | [Usage](scale-families.md#numeric) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_fill_stepsn) | [↗](chord.md) |
| [`scale_ribbon_linetype_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | Ribbon alpha, outline and linetype scales | [Usage](scale-families.md#categories) | chord/collinear matching ribbon roles | [↗](parameter-atlas.md#api-scale_ribbon_linetype_manual) | [↗](chord.md) |
| [`scale_seq_color_manual`](../../man/scale_seq_colour_manual.Rd) | Sequence colour scales; spelling alias of scale_seq_colour_manual | [Usage](scale-families.md#categories) | all three; matching sequence colour role | [↗](parameter-atlas.md#api-scale_seq_color_manual) | [↗](chord.md) |
| [`scale_seq_colour_manual`](../../man/scale_seq_colour_manual.Rd) | Sequence colour scales | [Usage](scale-families.md#categories) | all three; matching sequence colour role | [↗](parameter-atlas.md#api-scale_seq_colour_manual) | [↗](chord.md) |
| [`scale_seq_position_continuous`](../../man/scale_seq_position_continuous.Rd) | Genomic sequence position scale | [Usage](views.md#focus) | all three; position axis, separate from the collinear ruler | [↗](parameter-atlas.md#api-scale_seq_position_continuous) | [↗](chord.md) |
| [`scale_seq_ring_manual`](../../man/scale_seq_ring_manual.Rd) | Explicit sequence ring scale | [Usage](scale-families.md#rings) | chord/circular; positive declared radii; linear rejected | [↗](parameter-atlas.md#api-scale_seq_ring_manual) | [↗](chord.md) |
| [`scale_track_colour_gradientn`](../../man/scale_track_fill_gradientn.Rd) | Set a track&#x27;s fill or colour mapping | [Usage](scale-families.md#tracks) | all three; scale ID must match track_band | [↗](parameter-atlas.md#api-scale_track_colour_gradientn) | [↗](tracks.md) |
| [`scale_track_colour_manual`](../../man/scale_track_fill_gradientn.Rd) | Set a track&#x27;s fill or colour mapping | [Usage](scale-families.md#tracks) | all three; scale ID must match track_band | [↗](parameter-atlas.md#api-scale_track_colour_manual) | [↗](tracks.md) |
| [`scale_track_fill_gradientn`](../../man/scale_track_fill_gradientn.Rd) | Set a track&#x27;s fill or colour mapping | [Usage](scale-families.md#tracks) | all three; scale ID must match track_band | [↗](parameter-atlas.md#api-scale_track_fill_gradientn) | [↗](tracks.md) |
| [`scale_track_fill_manual`](../../man/scale_track_fill_gradientn.Rd) | Set a track&#x27;s fill or colour mapping | [Usage](scale-families.md#tracks) | all three; scale ID must match track_band | [↗](parameter-atlas.md#api-scale_track_fill_manual) | [↗](tracks.md) |
| [`scale_track_value_continuous`](../../man/scale_track_value_continuous.Rd) | Set a named track&#x27;s numerical height scale | [Usage](scale-families.md#tracks) | all three; scale ID must match track_band | [↗](parameter-atlas.md#api-scale_track_value_continuous) | [↗](tracks.md) |
| [`segment_ggchord_data`](../../man/segment_ggchord_data.Rd) | Prepare independent windows for a segmented linear display | [Usage](views.md#segments) | independent windows and declared local display reversal | [↗](parameter-atlas.md#api-segment_ggchord_data) | [↗](views.md) |
| [`stat_coverage`](../../man/stat_coverage.Rd) | Compute union coverage and mean depth for inclusive intervals in each window | [Usage](tracks-statistics.md#interval-coverage) | Source intervals with reference lengths and units; statistics within each observed group | [↗](parameter-atlas.md#api-stat_coverage) | [↗](tracks.md) |
| [`stat_ribbon_bundle`](../../man/stat_ribbon_bundle.Rd) | Bundle ribbons as a ggchord stat | [Usage](ribbons.md#statistics-bundle) | chord/collinear ribbon layer; explicit aggregation | [↗](parameter-atlas.md#api-stat_ribbon_bundle) | [↗](chord.md) |
| [`stat_ribbon_density`](../../man/stat_ribbon_density.Rd) | Compute local ribbon density as a ggchord stat | [Usage](ribbons.md#statistics-density) | chord/collinear ribbon layer; one ribbon per source observation | [↗](parameter-atlas.md#api-stat_ribbon_density) | [↗](chord.md) |
| [`theme`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-theme) | [↗](chord.md) |
| [`theme_ggchord`](../../man/theme_ggchord.Rd) | ggchord themes | [Usage](scale-families.md#themes) | non-data style; plasmid preset intended circular | [↗](parameter-atlas.md#api-theme_ggchord) | [↗](chord.md) |
| [`theme_ggchord_dark`](../../man/theme_ggchord.Rd) | ggchord themes | [Usage](scale-families.md#themes) | non-data style; plasmid preset intended circular | [↗](parameter-atlas.md#api-theme_ggchord_dark) | [↗](chord.md) |
| [`theme_ggchord_minimal`](../../man/theme_ggchord.Rd) | ggchord themes | [Usage](scale-families.md#themes) | non-data style; plasmid preset intended circular | [↗](parameter-atlas.md#api-theme_ggchord_minimal) | [↗](chord.md) |
| [`theme_ggchord_plasmid`](../../man/theme_ggchord.Rd) | ggchord themes | [Usage](scale-families.md#themes) | non-data style; plasmid preset intended circular | [↗](parameter-atlas.md#api-theme_ggchord_plasmid) | [↗](chord.md) |
| [`theme_ggchord_publication`](../../man/theme_ggchord.Rd) | ggchord themes | [Usage](scale-families.md#themes) | non-data style; plasmid preset intended circular | [↗](parameter-atlas.md#api-theme_ggchord_publication) | [↗](chord.md) |
| [`track_band`](../../man/track_band.Rd) | Declare a quantitative band along sequence backbones | [Usage](tracks-statistics.md#independent-bands) | all three; one placement and numerical domain per track ID | [↗](parameter-atlas.md#api-track_band) | [↗](tracks.md) |
| [`translate_ggchord_cds`](../../man/translate_ggchord_cds.Rd) | Translate real CDS strings with codon provenance using optional Biostrings | [Usage](aligned-sequences.md#cds-translation) | Explicit source/column map; independent Cartesian detail; complete observed codons | [↗](parameter-atlas.md#api-translate_ggchord_cds) | [↗](formats.md) |
| [`unit`](../../man/grid-helpers.Rd) | Physical units and direction arrows | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-unit) | [↗](chord.md) |
| [`validate_ggchord_data`](../../man/validate_ggchord_data.Rd) | Validate ggchord input data before plotting | [Usage](data.md#validation) | validation reports, no silent input repair | [↗](parameter-atlas.md#api-validate_ggchord_data) | [↗](real-data.md) |
| [`validate_ggchord_hierarchy`](../../man/validate_ggchord_hierarchy.Rd) | Validate explicit genome and contig relationships | [Usage](data.md#hierarchy) | explicit genome/contig membership; never inferred from ID prefixes | [↗](parameter-atlas.md#api-validate_ggchord_hierarchy) | [↗](real-data.md) |
| [`view_ggchord`](../../man/view_ggchord.Rd) | Preview a ggchord plot at its intended export size | [Usage](output.md#size) | target-size preview of all supported plots; temporary path output | [↗](parameter-atlas.md#api-view_ggchord) | [↗](chord.md) |
| [`waiver`](../../man/ggplot2-helpers.Rd) | Selected ggplot2 helpers for ggchord plots | [upstream](parameter-atlas.md#upstream) | upstream; applies only to matching supported ggplot roles | [↗](parameter-atlas.md#api-waiver) | [↗](chord.md) |

<a id="reproduce"></a>

## How to reproduce {#reproduce}

Start at the preparation code in the primary section. Inspect fields in [basic tables](data.md) and conversions in [formats](formats.md). Shared example snippets execute in independent environments and important comparisons are rendered by ggchord. Explicit synthetic geometry demonstrations are not biological evidence. Retain inputs, source and reproduction notes; see [installation](installation.md) for setup.

<details class="maintenance-record"><summary>Maintenance coverage record</summary>

L2 denotes minimal runnable use and input/output explanation; L3 adds key-parameter contrasts. Neither certifies every combination.

| API | Level |
| --- | --- |
| `adapt_ggchord` | L2 |
| `aes` | L1 |
| `after_scale` | L1 |
| `after_stat` | L1 |
| `annotate` | L1 |
| `arrow` | L1 |
| `as_ggchord_alignment` | L2 |
| `as_ggchord_feature` | L2 |
| `as_ggchord_link` | L2 |
| `as_ggchord_seq` | L2 |
| `build_ggchord_cds` | L2 |
| `bundle_ggchord_ribbons` | L3 |
| `clean_ggchord_data` | L2 |
| `coord_chord` | L3 |
| `coord_circular` | L3 |
| `coord_collinear` | L3 |
| `coord_dotplot` | L2 |
| `coord_radial_synteny` | L2 |
| `coord_sankey` | L2 |
| `deduplicate_ggchord_ribbons` | L2 |
| `diagnose_ggchord` | L2 |
| `element_blank` | L1 |
| `element_line` | L1 |
| `element_rect` | L1 |
| `element_text` | L1 |
| `expansion` | L1 |
| `export_ggchord_fusion` | L2 |
| `export_ggchord_layout` | L2 |
| `filter_ggchord_ribbons` | L2 |
| `filter_restriction_sites` | L2 |
| `find_common_features` | L2 |
| `find_ggchord_loci` | L2 |
| `find_primer_bindings` | L2 |
| `find_restriction_sites` | L2 |
| `focus_ggchord_data` | L2 |
| `focus_ggchord_track` | L2 |
| `fusion_panel` | L2 |
| `fusion_ribbon` | L2 |
| `geom_alignment_segment` | L2 |
| `geom_feature` | L3 |
| `geom_feature_label` | L3 |
| `geom_gene_label` | L3 |
| `geom_gene_label_repel` | L3 |
| `geom_link_line` | L3 |
| `geom_link_ribbon` | L3 |
| `geom_primer` | L3 |
| `geom_primer_label_repel` | L3 |
| `geom_residue` | L2 |
| `geom_restriction_site` | L3 |
| `geom_seq` | L3 |
| `geom_seq_break` | L2 |
| `geom_seq_center_label` | L3 |
| `geom_seq_label` | L3 |
| `geom_seq_region` | L3 |
| `geom_seq_scale` | L3 |
| `geom_track` | L2 |
| `geom_track_bar` | L3 |
| `geom_track_col` | L2 |
| `geom_track_line` | L3 |
| `geom_track_tile` | L3 |
| `get_chord_layout` | L2 |
| `ggchord` | L2 |
| `ggchord_fusion` | L2 |
| `ggsave` | L1 |
| `ggtitle` | L1 |
| `guide_ggchord_colorbar` | L2 |
| `guide_ggchord_colourbar` | L2 |
| `guide_ggchord_legend` | L2 |
| `guide_seq_axis` | L2 |
| `guide_track_axis` | L2 |
| `guides` | L1 |
| `labs` | L1 |
| `last_plot` | L1 |
| `layout_sequences` | L2 |
| `link_ggchord_clusters` | L2 |
| `magnify_ggchord` | L2 |
| `map_alignment_position` | L2 |
| `map_ggchord_position` | L2 |
| `margin` | L1 |
| `merge_ggchord_ribbons` | L2 |
| `optimize_ggchord_layout` | L2 |
| `position_feature_stack` | L3 |
| `position_identity` | L1 |
| `position_offset` | L2 |
| `position_strand` | L2 |
| `PositionIdentity` | L1 |
| `prepare_ggchord_data` | L2 |
| `primer_catalog` | L2 |
| `project_ggchord_features` | L2 |
| `project_ggchord_links` | L2 |
| `read_aligned_fasta` | L2 |
| `read_bed` | L2 |
| `read_bedgraph` | L2 |
| `read_bedpe` | L2 |
| `read_bigwig` | L2 |
| `read_blast` | L2 |
| `read_chain` | L2 |
| `read_delta` | L2 |
| `read_fasta_lengths` | L2 |
| `read_feature_relations` | L2 |
| `read_genbank` | L2 |
| `read_gff3` | L2 |
| `read_maf` | L2 |
| `read_paf` | L2 |
| `read_signature_matches` | L2 |
| `read_vcf` | L2 |
| `register_ggchord_alignment` | L2 |
| `rel` | L1 |
| `scale_feature_fill_manual` | L2 |
| `scale_feature_fill_plasmid` | L2 |
| `scale_feature_shape_manual` | L2 |
| `scale_feature_shape_plasmid` | L2 |
| `scale_fill_residue` | L2 |
| `scale_gene_fill_manual` | L2 |
| `scale_link_alpha_continuous` | L2 |
| `scale_link_alpha_manual` | L2 |
| `scale_link_color_gradientn` | L2 |
| `scale_link_color_manual` | L2 |
| `scale_link_colour_gradientn` | L2 |
| `scale_link_colour_manual` | L2 |
| `scale_link_linetype_manual` | L2 |
| `scale_link_linewidth_continuous` | L2 |
| `scale_link_linewidth_manual` | L2 |
| `scale_region_fill_manual` | L2 |
| `scale_ribbon_alpha_continuous` | L2 |
| `scale_ribbon_alpha_manual` | L2 |
| `scale_ribbon_color_manual` | L2 |
| `scale_ribbon_colour_manual` | L2 |
| `scale_ribbon_fill_gradientn` | L2 |
| `scale_ribbon_fill_identity` | L2 |
| `scale_ribbon_fill_manual` | L2 |
| `scale_ribbon_fill_stepsn` | L2 |
| `scale_ribbon_linetype_manual` | L2 |
| `scale_seq_color_manual` | L2 |
| `scale_seq_colour_manual` | L2 |
| `scale_seq_position_continuous` | L2 |
| `scale_seq_ring_manual` | L2 |
| `scale_track_colour_gradientn` | L2 |
| `scale_track_colour_manual` | L2 |
| `scale_track_fill_gradientn` | L2 |
| `scale_track_fill_manual` | L2 |
| `scale_track_value_continuous` | L2 |
| `segment_ggchord_data` | L2 |
| `stat_coverage` | L2 |
| `stat_ribbon_bundle` | L3 |
| `stat_ribbon_density` | L3 |
| `theme` | L1 |
| `theme_ggchord` | L2 |
| `theme_ggchord_dark` | L2 |
| `theme_ggchord_minimal` | L2 |
| `theme_ggchord_plasmid` | L2 |
| `theme_ggchord_publication` | L2 |
| `track_band` | L3 |
| `translate_ggchord_cds` | L2 |
| `unit` | L1 |
| `validate_ggchord_data` | L2 |
| `validate_ggchord_hierarchy` | L2 |
| `view_ggchord` | L2 |
| `waiver` | L1 |

</details>
