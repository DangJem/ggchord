# Package index

## Plot construction

Build and inspect ggchord plots.

- [`ggchord()`](ggchord.md)
  : ggchord: layered multi-sequence alignment chord diagrams for ggplot2
- [`coord_chord()`](coord_chord.md)
  : Chord diagram coordinate system
- [`get_chord_layout()`](get_chord_layout.md)
  : Get the chord layout from the package environment
- [`` `+`( ``*`<ggchord>`*`)`](plus-.ggchord.md)
  : Combine a ggchord plot with ggplot2 objects

## Methods and datasets

S3 methods and example data shipped with ggchord.

- [`as.data.frame(`*`<ggchord_validation>`*`)`](as.data.frame.ggchord_validation.md)
  : Coerce a validation result to a flat data.frame
- [`print(`*`<ggchord_clean>`*`)`](print.ggchord_clean.md)
  : Print a cleaned ggchord data result
- [`seq_data_example`](seq_data_example.md)
  : Example sequence data
- [`ribbon_data_example`](ribbon_data_example.md)
  : Example alignment data
- [`gene_data_example`](gene_data_example.md)
  : Example gene annotation data

## Layers

Stack geometry, annotation, axis and label layers.

- [`geom_seq()`](geom_seq.md)
  : Add a sequence arc layer
- [`geom_ribbon()`](geom_ribbon.md)
  : Add an alignment ribbon layer
- [`geom_gene()`](geom_gene.md)
  : Add a gene arrow layer
- [`geom_gene_label()`](geom_gene_label.md)
  : Add a gene label layer
- [`geom_gene_label_repel()`](geom_gene_label_repel.md)
  : Add an automatically arranged gene label layer
- [`geom_axis()`](geom_axis.md)
  : Add an axis layer
- [`geom_seq_label()`](geom_seq_label.md)
  : Add a sequence label layer
- [`geom_seq_region()`](geom_seq_region.md)
  : Highlight regions along sequence arcs
- [`geom_ribbon_highlight()`](geom_ribbon_highlight.md)
  : Highlight selected alignment ribbons
- [`geom_feature()`](geom_feature.md)
  : Draw generic genomic features

## Scales

Control sequence, ribbon, feature and genomic-position aesthetics.

- [`scale_seq_colour_manual()`](scale_seq_colour_manual.md)
  [`scale_seq_color_manual()`](scale_seq_colour_manual.md)
  : Sequence colour scales
- [`scale_ribbon_fill_stepsn()`](scale_ribbon_fill_stepsn.md)
  [`scale_ribbon_fill_gradientn()`](scale_ribbon_fill_stepsn.md)
  [`scale_ribbon_fill_manual()`](scale_ribbon_fill_stepsn.md)
  [`scale_ribbon_fill_identity()`](scale_ribbon_fill_stepsn.md)
  : Ribbon fill scales
- [`scale_ribbon_alpha_continuous()`](scale_ribbon_alpha_continuous.md)
  [`scale_ribbon_alpha_manual()`](scale_ribbon_alpha_continuous.md)
  [`scale_ribbon_colour_manual()`](scale_ribbon_alpha_continuous.md)
  [`scale_ribbon_color_manual()`](scale_ribbon_alpha_continuous.md)
  [`scale_ribbon_linetype_manual()`](scale_ribbon_alpha_continuous.md)
  : Ribbon alpha, outline and linetype scales
- [`scale_gene_fill_manual()`](scale_gene_fill_manual.md)
  [`scale_feature_fill_manual()`](scale_gene_fill_manual.md)
  [`scale_region_fill_manual()`](scale_gene_fill_manual.md)
  : Gene, feature and region fill scales
- [`scale_feature_shape_manual()`](scale_feature_shape_manual.md)
  : Feature shape scale
- [`scale_seq_position_continuous()`](scale_seq_position_continuous.md)
  : Genomic sequence position scale

## Themes and guides

Apply publication-ready themes and legend defaults.

- [`theme_ggchord()`](theme_ggchord.md)
  [`theme_ggchord_minimal()`](theme_ggchord.md)
  [`theme_ggchord_dark()`](theme_ggchord.md)
  [`theme_ggchord_publication()`](theme_ggchord.md)
  : ggchord themes
- [`guide_ggchord_legend()`](guide_ggchord_legend.md)
  [`guide_ggchord_colourbar()`](guide_ggchord_legend.md)
  : Guides for ggchord role aesthetics

## Data import

Read FASTA, BLAST and GFF3 files in R.

- [`read_fasta_lengths()`](read_fasta_lengths.md)
  : Read one or more FASTA files into seq_data format
- [`read_blast()`](read_blast.md)
  : Read one or more BLAST tabular output files into ribbon_data format
- [`read_gff3()`](read_gff3.md)
  : Read one or more GFF3 files into gene_data format

## Validation and cleaning

Check and repair ggchord input tables before plotting.

- [`validate_ggchord_data()`](validate_ggchord_data.md)
  : Validate ggchord input data before plotting
- [`clean_ggchord_data()`](clean_ggchord_data.md)
  : Clean ggchord input data with explicit, report-driven policies

## Ribbon preparation

Filter, deduplicate and merge alignment blocks.

- [`filter_ggchord_ribbons()`](filter_ggchord_ribbons.md)
  : Filter alignment ribbons before plotting
- [`deduplicate_ggchord_ribbons()`](deduplicate_ggchord_ribbons.md)
  : Deduplicate alignment ribbons
- [`merge_ggchord_ribbons()`](merge_ggchord_ribbons.md)
  : Merge adjacent or overlapping alignment blocks of the same sequence
  pair
