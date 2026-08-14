# Package index

## Plot construction

Build, convert and inspect ggchord plots.

- [`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md)
  : ggchord: layered multi-sequence alignment chord diagrams for ggplot2
- [`coord_chord()`](https://dangjem.github.io/ggchord/reference/coord_chord.md)
  : Chord diagram coordinate system
- [`get_chord_layout()`](https://dangjem.github.io/ggchord/reference/get_chord_layout.md)
  : Get the chord layout from the package environment
- [`` `+`( ``*`<ggchord>`*`)`](https://dangjem.github.io/ggchord/reference/plus-.ggchord.md)
  : Combine a ggchord plot with ggplot2 objects
- [`ggplotly(`*`<ggchord>`*`)`](https://dangjem.github.io/ggchord/reference/ggplotly.ggchord.md)
  : Convert a ggchord plot to a plotly object

## Methods and datasets

S3 methods and example data shipped with ggchord.

- [`as.data.frame(`*`<ggchord_validation>`*`)`](https://dangjem.github.io/ggchord/reference/as.data.frame.ggchord_validation.md)
  : Coerce a validation result to a flat data.frame
- [`print(`*`<ggchord_clean>`*`)`](https://dangjem.github.io/ggchord/reference/print.ggchord_clean.md)
  : Print a cleaned ggchord data result
- [`seq_data_example`](https://dangjem.github.io/ggchord/reference/seq_data_example.md)
  : Example sequence data
- [`ribbon_data_example`](https://dangjem.github.io/ggchord/reference/ribbon_data_example.md)
  : Example alignment data
- [`gene_data_example`](https://dangjem.github.io/ggchord/reference/gene_data_example.md)
  : Example gene annotation data

## Layers

Stack geometry, annotation, axis and label layers.

- [`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md)
  : Add a sequence arc layer
- [`geom_ribbon()`](https://dangjem.github.io/ggchord/reference/geom_ribbon.md)
  : Add an alignment ribbon layer
- [`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
  : Add a gene arrow layer
- [`geom_gene_label()`](https://dangjem.github.io/ggchord/reference/geom_gene_label.md)
  : Add a gene label layer
- [`geom_gene_label_repel()`](https://dangjem.github.io/ggchord/reference/geom_gene_label_repel.md)
  : Add a repelled gene label layer (ggrepel-style)
- [`geom_axis()`](https://dangjem.github.io/ggchord/reference/geom_axis.md)
  : Add an axis layer
- [`geom_seq_label()`](https://dangjem.github.io/ggchord/reference/geom_seq_label.md)
  : Add a sequence label layer
- [`geom_seq_region()`](https://dangjem.github.io/ggchord/reference/geom_seq_region.md)
  : Highlight regions along sequence arcs
- [`geom_ribbon_highlight()`](https://dangjem.github.io/ggchord/reference/geom_ribbon_highlight.md)
  : Highlight selected alignment ribbons
- [`geom_feature()`](https://dangjem.github.io/ggchord/reference/geom_feature.md)
  : Draw generic genomic features

## Data import

Read FASTA, BLAST and GFF3 files in R.

- [`read_fasta_lengths()`](https://dangjem.github.io/ggchord/reference/read_fasta_lengths.md)
  : Read one or more FASTA files into seq_data format
- [`read_blast()`](https://dangjem.github.io/ggchord/reference/read_blast.md)
  : Read one or more BLAST tabular output files into ribbon_data format
- [`read_gff3()`](https://dangjem.github.io/ggchord/reference/read_gff3.md)
  : Read one or more GFF3 files into gene_data format

## Validation and cleaning

Check and repair ggchord input tables before plotting.

- [`validate_ggchord_data()`](https://dangjem.github.io/ggchord/reference/validate_ggchord_data.md)
  : Validate ggchord input data before plotting
- [`clean_ggchord_data()`](https://dangjem.github.io/ggchord/reference/clean_ggchord_data.md)
  : Clean ggchord input data with explicit, report-driven policies

## Ribbon preparation

Filter, deduplicate and merge alignment blocks.

- [`filter_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/filter_ggchord_ribbons.md)
  : Filter alignment ribbons before plotting
- [`deduplicate_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/deduplicate_ggchord_ribbons.md)
  : Deduplicate alignment ribbons
- [`merge_ggchord_ribbons()`](https://dangjem.github.io/ggchord/reference/merge_ggchord_ribbons.md)
  : Merge adjacent or overlapping alignment blocks of the same sequence
  pair
