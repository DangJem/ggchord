## Future directions

ggchord will continue to evolve around **a shared model of genomic relationships, flexible ggplot2-style layers, and publication-ready static figures**. The following are planned or exploratory directions; **new functions and arguments mentioned here are not necessarily available yet**.

- **More genomic views**: Explore `coord_dotplot()`, spoke-based `coord_radial_synteny()`, and mirrored `coord_sankey()`, reusing sequence and interval-alignment tables to show the same relationships in different layouts.
- **Unified features and annotations**: Move toward `feature_data` as a common input for genes and other annotations. Extend `geom_feature()` through shapes, mappings, and positioning to cover variants, motifs, repeats, and protein domains, while retaining convenient layers such as `geom_gene()`.
- **Tracks with a familiar ggplot2 grammar**: Explore `geom_track(geom = ..., stat = ...)` for lines, points, areas, tiles, columns, and statistical transformations. For example, `geom = "bar", stat = "count"` would count observations, whereas `geom = "col", stat = "identity"` would plot supplied values; track positioning, overlays, and scales would be handled separately.
- **More expressive sequence relationships**: Improve ribbon endpoint colours, gradients, and directional arrowheads. Add optional equal-length sequence display and smarter first/last-position labels and intermediate ticks, without changing underlying bp/aa coordinates.
- **Magnified regions and sequence-level detail**: Explore coordinate-aware insets that connect whole-genome views to local features, nucleotide letters, or amino-acid residues, with independent layers and layouts in the detail panel.
- **Automatic layout and multi-view composition**: Improve dense annotations, leader lines, legends, and adaptation to different figure sizes. Explore figure compositions linking overview and detail panels through exact genomic coordinates.
- **A stable, extensible data foundation**: Keep ordinary `data.frame` inputs simple while standardizing sequences, annotations, alignments, and display instances internally. Continue to focus on reproducible PNG, PDF, and SVG output.
